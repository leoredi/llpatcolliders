"""Generate BC1 Drell--Yan spectra with native Pythia 8 + NNPDF4.0.

This is the production generator for the low-mass region where MadEvent's
``p p > zp`` 2->1 phase-space mapping returns NaNs. It uses Pythia's pure
``f fbar -> Z'`` process with vector couplings mapped to
``epsilon * e * Q_f`` and the same pinned NNPDF40 NLO central member used by
the HNL production chain. The showered A' pT/rapidity histogram and inclusive
cross section at epsilon^2=1 are written in the format consumed by
``dark_photon.production``.

The NNPDF set is frozen below Q=1.65 GeV, so this driver refuses lower masses.
Sub-1.65-GeV production must be covered by resolved hadronic channels rather
than an extrapolated perturbative DY calculation.

Usage:
    python -m dark_photon.generator.dy_pythia --masses 1.65 2 3 5 10
"""
from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

import numpy as np

from dark_photon.mass_grid import DY_MASS_GRID

_ROOT = Path(__file__).resolve().parents[1]
_SOURCE = Path(__file__).with_name("dp_dy_pythia.cc")
_BUILD_DIR = _ROOT / "tmp" / "pythia_dy"
_BINARY = _BUILD_DIR / "dp_dy_pythia"
_OUTPUT_DIR = _ROOT / "data" / "spectra" / "dy"

EPS_REF = 1.0e-2
PDF_QMIN = 1.65
PHOTON_LIKE_MAX = 50.0
PT_MIN_RANGE, PT_BIN_WIDTH = 30.0, 0.2
Y_MIN_RANGE, Y_BIN_WIDTH = 8.0, 0.1


def _find_workspace(rel: str) -> Path:
    for parent in Path(__file__).resolve().parents:
        candidate = parent / rel
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"workspace dependency not found: {rel}")


def _dependencies():
    pythia = _find_workspace("shared/vendored/pythia8315")
    hnl_env = _find_workspace("environments/conda/hnl")
    return pythia, hnl_env


def _build_binary() -> Path:
    pythia, hnl_env = _dependencies()
    if _BINARY.exists() and _BINARY.stat().st_mtime >= _SOURCE.stat().st_mtime:
        return _BINARY
    _BUILD_DIR.mkdir(parents=True, exist_ok=True)
    cmd = [
        "c++", str(_SOURCE), "-o", str(_BINARY), "-std=c++17", "-O2",
        f"-I{pythia / 'include'}", f"-I{hnl_env / 'include'}",
        f"-L{pythia / 'lib'}", f"-Wl,-rpath,{pythia / 'lib'}", "-lpythia8",
        f"-L{hnl_env / 'lib'}", f"-Wl,-rpath,{hnl_env / 'lib'}", "-lLHAPDF",
        "-ldl",
    ]
    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(f"Pythia DY build failed:\n{result.stdout}{result.stderr}")
    return _BINARY


def _mass_label(mass: float) -> str:
    return f"{mass:.3f}".replace(".", "p")


def _parse_metadata(stderr: str) -> dict[str, float]:
    keys = {
        "N_TRIED", "N_ACCEPTED", "N_ZP", "MASS_GEV", "EPSILON",
        "WIDTH_GEV", "SIGMA_MB", "SIGMA_ERR_MB",
    }
    metadata = {}
    for line in stderr.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0] in keys:
            metadata[parts[0]] = float(parts[1])
    missing = keys - metadata.keys()
    if missing:
        raise RuntimeError(f"Pythia DY output missing metadata: {sorted(missing)}")
    return metadata


def make_dy_spectrum(mass: float, n_events: int = 20_000, seed: int = 42) -> Path:
    if mass < PDF_QMIN:
        raise ValueError(
            f"m_A={mass:g} GeV is below the NNPDF4.0 validity floor "
            f"({PDF_QMIN:g} GeV); do not extrapolate perturbative DY"
        )
    if mass > PHOTON_LIKE_MAX:
        raise ValueError(
            f"m_A={mass:g} GeV is above the validated photon-like DY range "
            f"({PHOTON_LIKE_MAX:g} GeV); use the full electroweak HAHM model"
        )
    binary = _build_binary()
    pythia, hnl_env = _dependencies()
    env = os.environ.copy()
    env["PYTHIA8DATA"] = str(pythia / "share" / "Pythia8" / "xmldoc")
    env["LHAPDF_DATA_PATH"] = str(hnl_env / "share" / "LHAPDF")
    for key in ("DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH"):
        paths = [str(pythia / "lib"), str(hnl_env / "lib")]
        if env.get(key):
            paths.append(env[key])
        env[key] = os.pathsep.join(paths)

    result = subprocess.run(
        [str(binary), f"{mass:.12g}", str(n_events), str(EPS_REF), str(seed)],
        text=True, capture_output=True, env=env,
    )
    if result.returncode:
        raise RuntimeError(f"Pythia DY generation failed:\n{result.stderr[-4000:]}")
    metadata = _parse_metadata(result.stderr)
    rows = []
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) == 3 and parts[0] == "32":
            rows.append((float(parts[1]), float(parts[2])))
    if not rows or len(rows) != int(metadata["N_ZP"]):
        raise RuntimeError(
            f"Pythia reported {int(metadata['N_ZP'])} A' but parsed {len(rows)}"
        )
    points = np.asarray(rows)
    # Preserve the original 0.2-GeV / 0.1 rapidity resolution, but extend each
    # histogram far enough to retain every generated event. A fixed 30-GeV pT
    # ceiling discarded the high-pT tail most relevant to a transverse detector.
    pt_max = max(PT_MIN_RANGE,
                 np.ceil(points[:, 0].max() / 10.0) * 10.0 + PT_BIN_WIDTH)
    y_abs = max(Y_MIN_RANGE,
                np.ceil(np.abs(points[:, 1]).max()) + Y_BIN_WIDTH)
    pt_bins = int(np.ceil(pt_max / PT_BIN_WIDTH))
    y_bins_half = int(np.ceil(y_abs / Y_BIN_WIDTH))
    pt_edges = np.linspace(0.0, pt_bins * PT_BIN_WIDTH, pt_bins + 1)
    y_edges = np.linspace(-y_bins_half * Y_BIN_WIDTH,
                          y_bins_half * Y_BIN_WIDTH, 2 * y_bins_half + 1)
    hist, _, _ = np.histogram2d(points[:, 0], points[:, 1],
                                bins=[pt_edges, y_edges])
    sigma_mg_pb = metadata["SIGMA_MB"] * 1.0e9
    sigma_dy_pb = sigma_mg_pb / EPS_REF**2
    _OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output = _OUTPUT_DIR / f"dy_{_mass_label(mass)}.npz"
    np.savez_compressed(
        output,
        hist=hist,
        pt_edges=pt_edges,
        y_edges=y_edges,
        sigma_dy_pb=np.float64(sigma_dy_pb),
        n_sample=np.int64(len(points)),
        sigma_generator_pb=np.float64(sigma_mg_pb),
        sigma_generator_err_pb=np.float64(metadata["SIGMA_ERR_MB"] * 1.0e9),
        eps_ref=np.float64(EPS_REF),
        mass_GeV=np.float64(mass),
        width_GeV=np.float64(metadata["WIDTH_GEV"]),
        generator=np.str_("Pythia 8.315 NewGaugeBoson pure-Zprime"),
        pdf=np.str_("NNPDF40_nlo_as_01180 member 0"),
        seed=np.int64(seed),
        histogram_coverage=np.float64(hist.sum() / len(points)),
    )
    print(
        f"m_A={mass:.3f} GeV  n={len(points)}  "
        f"sigma(eps={EPS_REF:g})={sigma_mg_pb:.4g} pb  "
        f"sigma(eps^2=1)={sigma_dy_pb:.4g} pb  -> {output.name}",
        flush=True,
    )
    return output


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--masses", type=float, nargs="+", default=None,
        help="mass points in GeV (default: canonical DY grid)",
    )
    parser.add_argument("--n-events", type=int, default=20_000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args(argv)
    masses = args.masses if args.masses is not None else DY_MASS_GRID
    for mass in masses:
        # Stable under subset generation and future grid refinement.
        mass_seed = args.seed + int(round(mass * 1000.0))
        make_dy_spectrum(mass, args.n_events, mass_seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
