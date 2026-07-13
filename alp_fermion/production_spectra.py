"""SensCalc 2501 LHC light-meson spectra for the BC10 production audit.

SensCalc stores the physical parent multiplicity separately from each
normalized ``(theta, energy)`` shape.  This module preserves that separation:
the sampler returns parent four-vectors, while callers apply the per-collision
yield and the relevant ALP branching ratio as production weights.

The source checkout must first pass
``tools/export_senscalc_2501_production.py --check-only``.  The values in
``LIGHT_MESON_PARENTS`` are the LHC entries from the pinned
``codes/experiments.nb`` file.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .paths import ANALYSIS_DIR


@dataclass(frozen=True)
class ParentSpectrumSpec:
    """One SensCalc LHC parent shape and its multiplicity per pp collision."""

    process_name: str
    mass_gev: float
    yield_per_collision: float
    filename: str


LIGHT_MESON_PARENTS = {
    "Eta": ParentSpectrumSpec(
        "Eta", 0.547862, 3.64, "DoubleDistr_LHC_Eta.txt"
    ),
    "EtaPr": ParentSpectrumSpec(
        "EtaPr", 0.95778, 0.46, "DoubleDistr_LHC_EtaPr.txt"
    ),
    "Omega": ParentSpectrumSpec(
        "Omega", 0.78266, 4.42, "DoubleDistr_LHC_Omega.txt"
    ),
    "RhoCh": ParentSpectrumSpec(
        "RhoCh", 0.77511, 8.68181, "DoubleDistr_LHC_RhoCh.txt"
    ),
    "KS": ParentSpectrumSpec(
        "KS", 0.497611, 3.1, "DoubleDistr_LHC_KS.txt"
    ),
}


def _node_bin_edges(values: np.ndarray) -> np.ndarray:
    if values.ndim != 1 or len(values) < 2:
        raise ValueError("spectrum sampling requires at least two grid nodes")
    if np.any(np.diff(values) <= 0):
        raise ValueError("spectrum grid nodes must be strictly increasing")
    edges = np.empty(len(values) + 1, dtype=float)
    edges[0] = values[0]
    edges[-1] = values[-1]
    edges[1:-1] = 0.5 * (values[:-1] + values[1:])
    return edges


def _refined_log_nodes(
    values: np.ndarray,
    subdivisions: int,
    extra_nodes: tuple[float, ...] = (),
) -> np.ndarray:
    if subdivisions < 1:
        raise ValueError("subdivisions must be positive")
    extra = [value for value in extra_nodes if values[0] < value < values[-1]]
    base = np.unique(np.concatenate((values, np.asarray(extra, dtype=float))))
    if np.any(base <= 0.0):
        raise ValueError("log-interpolated spectrum nodes must be positive")
    pieces = [
        np.geomspace(base[i], base[i + 1], subdivisions + 1)[:-1]
        for i in range(len(base) - 1)
    ]
    return np.concatenate((*pieces, base[-1:]))


@dataclass(frozen=True)
class TabulatedSpectrum:
    """A validated Cartesian density grid in polar angle and energy."""

    theta: np.ndarray
    energy: np.ndarray
    density: np.ndarray
    source: Path

    @classmethod
    def load(cls, path: Path) -> "TabulatedSpectrum":
        path = Path(path)
        data = np.loadtxt(path, dtype=float)
        if data.ndim != 2 or data.shape[1] != 3:
            raise ValueError(f"expected three spectrum columns in {path}")
        if not np.all(np.isfinite(data)):
            raise ValueError(f"non-finite spectrum value in {path}")

        theta = np.unique(data[:, 0])
        energy = np.unique(data[:, 1])
        if len(data) != len(theta) * len(energy):
            raise ValueError(f"spectrum is not a Cartesian grid: {path}")
        theta_mesh, energy_mesh = np.meshgrid(theta, energy, indexing="ij")
        if not (
            np.allclose(data[:, 0], theta_mesh.ravel(), rtol=0.0, atol=1e-12)
            and np.allclose(
                data[:, 1], energy_mesh.ravel(), rtol=0.0, atol=1e-12
            )
        ):
            raise ValueError(f"spectrum rows are not in theta-energy order: {path}")

        density = data[:, 2].reshape(len(theta), len(energy))
        if np.any(density < 0):
            raise ValueError(f"negative spectrum density in {path}")
        return cls(theta, energy, density, path.resolve())

    @property
    def theta_edges(self) -> np.ndarray:
        return _node_bin_edges(self.theta)

    @property
    def energy_edges(self) -> np.ndarray:
        return _node_bin_edges(self.energy)

    def cell_weights(self) -> np.ndarray:
        theta_width = np.diff(self.theta_edges)
        energy_width = np.diff(self.energy_edges)
        return self.density * theta_width[:, None] * energy_width[None, :]

    def _smooth_density(
        self, theta_nodes: np.ndarray, energy_nodes: np.ndarray
    ) -> np.ndarray:
        log_density = np.log(np.where(self.density > 0.0, self.density, 1e-90))
        interpolator = RegularGridInterpolator(
            (np.log(self.theta), np.log(self.energy)),
            log_density,
            method="linear",
            bounds_error=False,
            fill_value=-np.inf,
        )
        theta_mesh, energy_mesh = np.meshgrid(
            theta_nodes, energy_nodes, indexing="ij"
        )
        points = np.column_stack((
            np.log(theta_mesh.ravel()),
            np.log(energy_mesh.ravel()),
        ))
        return np.exp(interpolator(points)).reshape(theta_mesh.shape)

    def _smooth_integral(
        self, theta_nodes: np.ndarray, energy_nodes: np.ndarray
    ) -> float:
        density = self._smooth_density(theta_nodes, energy_nodes)
        energy_integral = np.trapz(density, energy_nodes, axis=1)
        return float(np.trapz(energy_integral, theta_nodes))

    def integral(self, subdivisions: int = 4) -> float:
        """Integrate the same log-linear interpolant used by SensCalc."""
        theta_nodes = _refined_log_nodes(self.theta, subdivisions)
        energy_nodes = _refined_log_nodes(self.energy, subdivisions)
        return self._smooth_integral(theta_nodes, energy_nodes)

    def polar_fraction(
        self, theta_min: float, theta_max: float, subdivisions: int = 4
    ) -> float:
        """Integrate the SensCalc log-linear shape over a polar interval."""
        if not 0.0 <= theta_min < theta_max <= math.pi:
            raise ValueError("polar interval must lie within [0, pi]")
        clipped_min = max(theta_min, float(self.theta[0]))
        clipped_max = min(theta_max, float(self.theta[-1]))
        if clipped_min >= clipped_max:
            return 0.0

        selected_base = self.theta[
            (self.theta > clipped_min) & (self.theta < clipped_max)
        ]
        selected_theta = _refined_log_nodes(
            np.concatenate(([clipped_min], selected_base, [clipped_max])),
            subdivisions,
        )
        energy_nodes = _refined_log_nodes(self.energy, subdivisions)
        selected = self._smooth_integral(selected_theta, energy_nodes)
        norm = self.integral(subdivisions)
        if not math.isfinite(norm) or norm <= 0.0:
            raise ValueError(f"spectrum has no finite positive weight: {self.source}")
        return selected / norm

    def sample(
        self,
        n_events: int,
        parent_mass: float,
        rng: np.random.Generator | None = None,
        subdivisions: int = 4,
    ) -> dict[str, np.ndarray]:
        """Sample on-shell parents from the SensCalc log-linear grid shape."""
        if n_events <= 0:
            raise ValueError("n_events must be positive")
        if parent_mass <= 0.0:
            raise ValueError("parent_mass must be positive")
        if rng is None:
            rng = np.random.default_rng()

        if parent_mass >= self.energy[-1]:
            raise ValueError(
                f"parent mass outside the spectrum energy range: {self.source}"
            )
        theta_nodes = _refined_log_nodes(self.theta, subdivisions)
        energy_nodes = _refined_log_nodes(
            self.energy, subdivisions, extra_nodes=(parent_mass,)
        )
        density = self._smooth_density(theta_nodes, energy_nodes)
        density[:, energy_nodes < parent_mass] = 0.0

        theta_edges = _node_bin_edges(theta_nodes)
        energy_edges = _node_bin_edges(energy_nodes)
        energy_lower = np.maximum(energy_edges[:-1], parent_mass)
        energy_upper = np.maximum(energy_edges[1:], parent_mass)
        energy_width = np.maximum(energy_upper - energy_lower, 0.0)
        theta_width = np.diff(theta_edges)
        weights = (density * theta_width[:, None] * energy_width[None, :]).ravel()
        norm = float(weights.sum())
        if not math.isfinite(norm) or norm <= 0.0:
            raise ValueError(f"spectrum has no finite positive weight: {self.source}")
        cdf = np.cumsum(weights)
        cdf /= cdf[-1]

        flat_index = np.searchsorted(cdf, rng.random(n_events), side="right")
        flat_index = np.clip(flat_index, 0, weights.size - 1)
        energy_count = len(energy_nodes)
        theta_index = flat_index // energy_count
        energy_index = flat_index % energy_count

        theta = rng.uniform(
            theta_edges[theta_index], theta_edges[theta_index + 1]
        )
        energy = rng.uniform(
            energy_lower[energy_index], energy_upper[energy_index]
        )
        if np.any(energy < parent_mass - 1e-12):
            raise ValueError(
                f"spectrum energy below parent mass {parent_mass}: {self.source}"
            )

        momentum = np.sqrt(np.maximum(energy**2 - parent_mass**2, 0.0))
        phi = rng.uniform(0.0, 2.0 * math.pi, n_events)
        transverse = momentum * np.sin(theta)
        return {
            "E": energy,
            "px": transverse * np.cos(phi),
            "py": transverse * np.sin(phi),
            "pz": momentum * np.cos(theta),
            "theta": theta,
            "phi": phi,
        }


def spectrum_path(senscalc_root: Path, spec: ParentSpectrumSpec) -> Path:
    return (
        Path(senscalc_root)
        / "spectra"
        / "SM particles"
        / spec.filename
    )


def transverse_theta_interval(eta_max: float) -> tuple[float, float]:
    """Return the polar interval equivalent to ``|eta| < eta_max``."""
    if eta_max <= 0.0:
        raise ValueError("eta_max must be positive")
    theta_min = 2.0 * math.atan(math.exp(-eta_max))
    return theta_min, math.pi - theta_min


def summarize_parent_spectra(
    senscalc_root: Path, eta_max: float = 0.5
) -> list[dict[str, float | str]]:
    theta_min, theta_max = transverse_theta_interval(eta_max)
    rows = []
    for spec in LIGHT_MESON_PARENTS.values():
        spectrum = TabulatedSpectrum.load(spectrum_path(senscalc_root, spec))
        fraction = spectrum.polar_fraction(theta_min, theta_max)
        rows.append({
            "parent": spec.process_name,
            "mass_GeV": spec.mass_gev,
            "yield_per_collision": spec.yield_per_collision,
            "shape_integral": spectrum.integral(),
            "parent_fraction_abs_eta_lt": fraction,
            "parent_yield_abs_eta_lt_per_collision": (
                spec.yield_per_collision * fraction
            ),
        })
    return rows


def write_summary(
    rows: list[dict[str, float | str]], output_path: Path, eta_max: float
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "parent",
        "mass_GeV",
        "yield_per_collision",
        "shape_integral",
        f"parent_fraction_abs_eta_lt_{eta_max:g}",
        f"parent_yield_abs_eta_lt_{eta_max:g}_per_collision",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                fieldnames[0]: row["parent"],
                fieldnames[1]: row["mass_GeV"],
                fieldnames[2]: row["yield_per_collision"],
                fieldnames[3]: row["shape_integral"],
                fieldnames[4]: row["parent_fraction_abs_eta_lt"],
                fieldnames[5]: row["parent_yield_abs_eta_lt_per_collision"],
            })


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("senscalc_root", type=Path)
    parser.add_argument("--eta-max", type=float, default=0.5)
    parser.add_argument(
        "--output",
        type=Path,
        default=ANALYSIS_DIR / "senscalc_2501_parent_spectra.csv",
    )
    args = parser.parse_args(argv)

    from .tools.export_senscalc_2501_production import verify_sources

    verify_sources(args.senscalc_root)
    rows = summarize_parent_spectra(args.senscalc_root, args.eta_max)
    write_summary(rows, args.output, args.eta_max)
    for row in rows:
        print(
            f"{row['parent']:6s}  yield={row['yield_per_collision']:.6g}  "
            f"parent |eta|<{args.eta_max:g} fraction="
            f"{row['parent_fraction_abs_eta_lt']:.6f}"
        )
    print(f"Wrote {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
