#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
ENV = ROOT / "env"
FONLL = ROOT / "src" / "fonll"
LINUX = FONLL / "Linux"
OUTDIR = ROOT / "output"
RUNDIR = ROOT / "run"

PDFS = {
    "cteq66": ("CTEQ6.6", 10550, "cteq66"),
    "nlo": ("NNPDF40_nlo_as_01180", 331700),
    "nlo_as_01170": ("NNPDF40_nlo_as_01170", 333900, "nnpdf40_nlo_as_01170"),
    "nlo_as_01190": ("NNPDF40_nlo_as_01190", 334100, "nnpdf40_nlo_as_01190"),
    "nnlo": ("NNPDF40_nnlo_as_01180", 331100),
}

QUARKS = {
    "bottom": {
        "mass": 4.75,
        "frag_mode": 2,
        "frag_param": 24.2,
        "frag_param_name": "alpha",
        "fragmentation": "Kartvelishvili, normalized (1-z)*z^alpha; public FONLL N=5 central default for B hadrons",
        "final_state_detail": "B hadron",
    },
    "charm": {
        "mass": 1.50,
        "frag_mode": 5,
        "frag_param": 0.1,
        "frag_param_name": "r",
        "fragmentation": "BCFY pseudoscalar plus calibrated D* feeddown; public FONLL D0 convention",
        "final_state_detail": "D0 with public FONLL D* feeddown convention",
        "feeddown": "public_d0",
    },
}

# Match the legacy production tables so changing the PDF set does not silently
# alter the sampler resolution used by downstream detector-acceptance studies.
Y_VALUES = [round(-3.0 + 6.0 * i / 99.0, 10) for i in range(100)]
PT_VALUES = [round(50.0 * i / 99.0, 10) for i in range(100)]

# The internal quark grid is wider than the requested meson grid because
# fragmentation samples quark pT above the observed meson pT.
FONLL_GRID_PTMAX = 200.0
FONLL_GRID_NPT = 80
FRAG_FRAME = 1
FRAG_MESON_MASS = -1.0
MATCHING_C = 5.0

# ---------------------------------------------------------------------------
# Production-uncertainty variation axes.
#
# Standard 7-point scale variation, expressed as (muR/mu0, muF/mu0) with
# mu0 = sqrt(m^2 + pT^2). The two antipodal extremes (2, 0.5) and (0.5, 2)
# are excluded by convention. FONLL's grid driver reads the two scale factors
# in the order (ffact=muF, fren=muR); see misc1/fonllgrid.f:120
# `read(*,*) ffact,fren` and main/fonll0.f (ffact -> xmu_fact, fren -> xmu_ren).
SCALE_POINTS_7 = [
    (1.0, 1.0),
    (2.0, 2.0),
    (0.5, 0.5),
    (2.0, 1.0),
    (1.0, 2.0),
    (0.5, 1.0),
    (1.0, 0.5),
]

# Heavy-quark-mass variations (GeV), FONLL benchmark conventions
# (Cacciari et al.): m_b = 4.75 +/- 0.25, m_c = 1.5 -0.2/+0.2.
MASS_VARIATIONS = {
    "bottom": [4.50, 5.00],
    "charm": [1.30, 1.70],
}

# Set by the campaign driver: gzip each per-command log on success so a full
# replica campaign does not leave gigabytes of uncompressed FONLL output.
_COMPRESS_LOGS = False


@dataclass(frozen=True)
class Variation:
    """One point in the (scale, PDF member, heavy-quark mass) variation space.

    ``muR``/``muF`` are the renormalization/factorization scale factors relative
    to mu0; ``pdf_member`` is the offset added to the set's base LHAPDF id;
    ``mass`` is the heavy-quark mass actually used (already resolved per quark).
    """

    kind: str  # "central" | "scale" | "pdf" | "mass"
    tag: str  # filesystem-safe label embedded in filenames and the manifest
    muR: float
    muF: float
    pdf_member: int
    mass: float


def _fmt_float_tag(value: float) -> str:
    """Filesystem-safe float, e.g. 2.0 -> '2p0', 0.5 -> '0p5', 4.75 -> '4p75'."""
    return f"{value:g}".replace(".", "p").replace("-", "m")


def enumerate_variations(
    quark: str,
    pdf_members: Iterable[int],
    include_scale: bool,
    include_mass: bool,
) -> list[Variation]:
    """Build the ordered variation list for one quark.

    The ``central`` point (member 0, central scale, central mass) is emitted
    once and reused downstream as the scale (1,1) point, the PDF member-0
    replica, and the mass-central point; the other axes only add their
    non-central points so no grid is computed twice.
    """
    default_mass = float(QUARKS[quark]["mass"])
    variations = [Variation("central", "central", 1.0, 1.0, 0, default_mass)]
    if include_scale:
        for muR, muF in SCALE_POINTS_7:
            if (muR, muF) == (1.0, 1.0):
                continue
            tag = f"scale_muR{_fmt_float_tag(muR)}_muF{_fmt_float_tag(muF)}"
            variations.append(Variation("scale", tag, muR, muF, 0, default_mass))
    for member in pdf_members:
        if member == 0:
            continue
        tag = f"pdf_{member:04d}"
        variations.append(Variation("pdf", tag, 1.0, 1.0, int(member), default_mass))
    if include_mass:
        for m in MASS_VARIATIONS[quark]:
            tag = f"mass_{_fmt_float_tag(m)}"
            variations.append(Variation("mass", tag, 1.0, 1.0, 0, float(m)))
    return variations


def point_key(pt: float, y: float) -> tuple[float, float]:
    # fragmfonll perturbs dense-grid coordinates at its printed precision.
    # Five decimals remain far below the node spacing and preserve all
    # distinct charm feeddown requests.
    return (round(pt, 5), round(y, 5))


def daughter_energy_fraction(parent_mass: float, daughter_mass: float, spectator_mass: float) -> float:
    return (
        parent_mass * parent_mass
        + daughter_mass * daughter_mass
        - spectator_mass * spectator_mass
    ) / (2.0 * parent_mass * parent_mass)


# Public FONLL v1.3.2 CTEQ6.6 D0 references at 14 TeV, y=0,
# queried from the public form on 2026-05-29. These calibrate only the
# charm D* feeddown weight; the generated central grid still uses the
# requested LHAPDF set.
PUBLIC_CHARM_D0_CTEQ66_Y0 = {
    1.0: 2.5497e8,
    5.0: 3.5221e7,
    8.0: 6.8479e6,
    15.0: 5.0952e5,
    22.0: 9.0635e4,
    29.0: 2.4923e4,
    36.0: 8.8653e3,
    43.0: 3.7354e3,
    50.0: 1.7742e3,
}

# Masses in GeV and branching fractions used for the collinear two-body
# D* feeddown approximation. The final normalization is fixed by the
# public-FONLL CTEQ6.6 calibration above.
CHARM_D0_FEEDDOWN_CHANNELS = [
    {
        "name": "D*0 -> D0 pi0",
        "branching_fraction": 0.647,
        "parent_mass": 2.00685,
        "daughter_mass": 1.86484,
        "spectator_mass": 0.134977,
    },
    {
        "name": "D*0 -> D0 gamma",
        "branching_fraction": 0.353,
        "parent_mass": 2.00685,
        "daughter_mass": 1.86484,
        "spectator_mass": 0.0,
    },
    {
        "name": "D*+ -> D0 pi+",
        "branching_fraction": 0.677,
        "parent_mass": 2.01026,
        "daughter_mass": 1.86484,
        "spectator_mass": 0.139570,
    },
]


def charm_d0_feeddown_channels() -> list[dict[str, float | str]]:
    channels = []
    for channel in CHARM_D0_FEEDDOWN_CHANNELS:
        enriched = dict(channel)
        enriched["z_collinear"] = daughter_energy_fraction(
            float(channel["parent_mass"]),
            float(channel["daughter_mass"]),
            float(channel["spectator_mass"]),
        )
        channels.append(enriched)
    return channels


def run_command(
    cmd: list[str],
    cwd: Path,
    stdin: str,
    log_path: Path,
    env: dict[str, str],
) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        proc = subprocess.run(
            cmd,
            cwd=cwd,
            input=stdin,
            text=True,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
            check=False,
        )
    if proc.returncode != 0:
        raise RuntimeError(f"command failed with exit {proc.returncode}: {' '.join(cmd)}; see {log_path}")
    if _COMPRESS_LOGS:
        # Keep failing logs uncompressed (the RuntimeError above points at the
        # plain .log); only successful, bulky logs are gzipped.
        gz_path = log_path.with_suffix(log_path.suffix + ".gz")
        with log_path.open("rb") as src, gzip.open(gz_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        log_path.unlink()


def fonll_env() -> dict[str, str]:
    env = os.environ.copy()
    env["PATH"] = f"{ENV / 'bin'}:{env.get('PATH', '')}"
    env["LHAPDF_DATA_PATH"] = str(ENV / "share" / "LHAPDF")
    return env


def grid_input(
    prefix: str,
    pdf_id: int,
    mass: float,
    ffact: float = 1.0,
    fren: float = 1.0,
    y_values: Iterable[float] = Y_VALUES,
    ptmax: float = FONLL_GRID_PTMAX,
    npt: int = FONLL_GRID_NPT,
) -> str:
    lines = [
        prefix,
        f" 1 7000. 0 0 {pdf_id}",
        f" 1 7000. 0 0 {pdf_id}",
        f" {mass:.8g}",
        " -1.",
        # FONLL grid driver order is (ffact=muF, fren=muR); see fonllgrid.f:120.
        f" {ffact:.8g} {fren:.8g}",
    ]
    lines.extend(f" {y:.8g}" for y in y_values)
    lines.extend(
        [
            " 10000",
            " 0",
            f" {ptmax:.8g} {npt}",
            " 1",
        ]
    )
    return "\n".join(lines) + "\n"


def _contiguous_chunks(values: list[float], n_chunks: int) -> list[list[float]]:
    n_chunks = min(n_chunks, len(values))
    size, remainder = divmod(len(values), n_chunks)
    chunks = []
    start = 0
    for index in range(n_chunks):
        stop = start + size + (1 if index < remainder else 0)
        chunks.append(values[start:stop])
        start = stop
    return chunks


def build_quark_grid(
    run_dir: Path,
    grid_file: Path,
    prefix: str,
    pdf_id: int,
    mass: float,
    pdf_key: str,
    quark: str,
    env: dict[str, str],
    grid_workers: int,
    ffact: float = 1.0,
    fren: float = 1.0,
    log_label: str = "",
) -> None:
    suffix = f"_{log_label}" if log_label else ""
    if grid_workers == 1:
        run_command(
            [str(LINUX / "fonllgridlha")],
            run_dir,
            grid_input(prefix, pdf_id, mass, ffact=ffact, fren=fren),
            ROOT / "logs" / f"fonllgrid_{pdf_key}_{quark}{suffix}.log",
            env,
        )
        return

    chunk_root = run_dir / "grid_chunks"
    if chunk_root.exists():
        shutil.rmtree(chunk_root)
    chunk_root.mkdir()
    chunks = _contiguous_chunks(Y_VALUES, grid_workers)

    def run_chunk(index: int, y_values: list[float]) -> Path:
        chunk_dir = chunk_root / f"{index:02d}"
        chunk_dir.mkdir()
        run_command(
            [str(LINUX / "fonllgridlha")],
            chunk_dir,
            grid_input(prefix, pdf_id, mass, ffact=ffact, fren=fren, y_values=y_values),
            ROOT / "logs" / f"fonllgrid_{pdf_key}_{quark}{suffix}_chunk{index:02d}.log",
            env,
        )
        chunk_file = chunk_dir / f"{prefix}.out"
        if not chunk_file.exists():
            raise RuntimeError(f"expected chunk grid file not found: {chunk_file}")
        return chunk_file

    with ThreadPoolExecutor(max_workers=len(chunks)) as executor:
        chunk_files = list(executor.map(run_chunk, range(len(chunks)), chunks))

    combined_lines = []
    header = None
    for chunk_file in chunk_files:
        lines = chunk_file.read_text().splitlines()
        if not lines:
            raise RuntimeError(f"empty chunk grid file: {chunk_file}")
        if header is None:
            header = lines[0]
            combined_lines.append(header)
        elif lines[0] != header:
            raise RuntimeError(f"inconsistent raw-grid header in {chunk_file}")
        combined_lines.extend(lines[1:])

    expected_lines = 1 + len(Y_VALUES) * FONLL_GRID_NPT
    if len(combined_lines) != expected_lines:
        raise RuntimeError(
            f"merged raw grid has {len(combined_lines)} lines, expected {expected_lines}"
        )
    grid_file.write_text("\n".join(combined_lines) + "\n")


def observed_points() -> list[tuple[float, float]]:
    return [
        (pt, y)
        for pt in PT_VALUES
        if pt > 0.0
        for y in Y_VALUES
    ]


def charm_feeddown_parent_points(points: Iterable[tuple[float, float]]) -> list[tuple[float, float]]:
    parent_points = set()
    for pt, y in points:
        for channel in charm_d0_feeddown_channels():
            z = float(channel["z_collinear"])
            parent_points.add(point_key(pt / z, y))
    return sorted(parent_points)


def frag_input_points(
    grid_file_name: str,
    frag_mode: int,
    frag_param: float,
    points: Iterable[tuple[float, float]],
) -> str:
    lines = [
        "0",
        "0",
        f"{MATCHING_C:.8g}",
        grid_file_name,
        "",
        "0",
        str(frag_mode),
        str(FRAG_FRAME),
        f"{FRAG_MESON_MASS:.8g}",
        f"{frag_param:.8g}",
        "9",
    ]
    for pt, y in points:
        lines.append(f"{pt:.10g} {y:.8g}")
    lines.append("-1 0")
    return "\n".join(lines) + "\n"


def frag_input(grid_file_name: str, frag_mode: int, frag_param: float) -> str:
    return frag_input_points(grid_file_name, frag_mode, frag_param, observed_points())


def parse_frag(path: Path) -> dict[tuple[float, float], float]:
    values: dict[tuple[float, float], float] = {}
    for line in path.read_text().splitlines()[1:]:
        if not line.strip():
            continue
        cols = line.split()
        if len(cols) < 3:
            continue
        pt = float(cols[0].replace("D", "E"))
        y = float(cols[1].replace("D", "E"))
        ds = float(cols[2].replace("D", "E"))
        values[point_key(pt, y)] = ds
    return values


def trapz2(values: dict[tuple[float, float], float]) -> float:
    total = 0.0
    for ipt, pt in enumerate(PT_VALUES):
        wpt = 0.5 if ipt == 0 or ipt == len(PT_VALUES) - 1 else 1.0
        for iy, y in enumerate(Y_VALUES):
            wy = 0.5 if iy == 0 or iy == len(Y_VALUES) - 1 else 1.0
            total += wpt * wy * values[point_key(pt, y)]
    return total * (PT_VALUES[1] - PT_VALUES[0]) * (Y_VALUES[1] - Y_VALUES[0])


def write_final(
    out_path: Path,
    pdf_name: str,
    pdf_id: int,
    quark: str,
    quark_config: dict[str, object],
    values: dict[tuple[float, float], float],
    grid_file: Path,
    extra_headers: dict[str, object] | None = None,
    mass: float | None = None,
    ffact: float = 1.0,
    fren: float = 1.0,
) -> dict[str, object]:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    mass = float(quark_config["mass"]) if mass is None else float(mass)
    for y in Y_VALUES:
        values[point_key(0.0, y)] = 0.0

    rows = []
    for pt in PT_VALUES:
        for y in Y_VALUES:
            key = point_key(pt, y)
            if key not in values:
                raise RuntimeError(f"missing fragmented point pT={pt}, y={y}")
            val = values[key]
            if not math.isfinite(val):
                raise RuntimeError(f"non-finite value at pT={pt}, y={y}: {val}")
            rows.append((pt, y, val))

    integral = trapz2(values)
    with out_path.open("w") as out:
        out.write("# FONLL heavy-flavor meson grid\n")
        out.write("# columns: pT y dsigma/dpT/dy\n")
        out.write("# units: GeV 1 pb/GeV\n")
        out.write("# collision: pp\n")
        out.write("# ebeam1_GeV: 7000\n")
        out.write("# ebeam2_GeV: 7000\n")
        out.write("# sqrt_s_GeV: 14000\n")
        out.write(f"# quark: {quark}\n")
        out.write("# final_state: meson\n")
        out.write(f"# final_state_detail: {quark_config['final_state_detail']}\n")
        out.write(f"# heavy_quark_mass_GeV: {mass:.8g}\n")
        out.write(
            f"# scale: mu0=sqrt(m^2+pT^2); ffact(muF)={ffact:.8g}, fren(muR)={fren:.8g}\n"
        )
        out.write(f"# pdf: {pdf_name}\n")
        out.write(f"# lhapdf_id: {pdf_id}\n")
        out.write("# perturbative_order: FONLL\n")
        out.write(f"# fragmentation_function: {quark_config['fragmentation']}\n")
        out.write(f"# fragmentation_mode: {quark_config['frag_mode']}\n")
        out.write(
            f"# fragmentation_{quark_config['frag_param_name']}: "
            f"{quark_config['frag_param']:.8g}\n"
        )
        if extra_headers:
            for key, value in extra_headers.items():
                out.write(f"# {key}: {value}\n")
        out.write("# fragmentation_frame: y=0\n")
        out.write("# fragmentation_fraction: 1\n")
        out.write("# meson_mass: heavy quark mass\n")
        out.write(f"# internal_quark_grid: {Path(grid_file).name}\n")
        out.write(f"# output_pT_values_GeV: 0..50 step {PT_VALUES[1] - PT_VALUES[0]:.8g}\n")
        out.write(f"# output_y_values: -3..3 step {Y_VALUES[1] - Y_VALUES[0]:.8g}\n")
        out.write("# note_pT0: pT=0 rows are set to 0; fragmfonll evaluates differential points for pT>0\n")
        out.write(f"# trapezoid_integral_pb_y-3to3_pt0to50: {integral:.12e}\n")
        out.write("# pT y dsigma/dpT/dy\n")
        for pt, y, val in rows:
            out.write(f"{pt:.8g} {y:.8g} {val:.12e}\n")

    vals = [r[2] for r in rows]
    return {
        "path": str(out_path),
        "rows": len(rows),
        "pt_min": min(r[0] for r in rows),
        "pt_max": max(r[0] for r in rows),
        "y_min": min(r[1] for r in rows),
        "y_max": max(r[1] for r in rows),
        "min_dsigma_dpT_dy_pb_per_GeV": min(vals),
        "max_dsigma_dpT_dy_pb_per_GeV": max(vals),
        "trapezoid_integral_pb": integral,
    }


def run_fragmentation(
    pdf_key: str,
    quark: str,
    run_dir: Path,
    grid_file: Path,
    frag_mode: int,
    frag_param: float,
    points: list[tuple[float, float]],
    label: str,
    env: dict[str, str],
) -> tuple[dict[tuple[float, float], float], Path]:
    frag_path = run_dir / "fragmfonll.dat"
    if frag_path.exists():
        frag_path.unlink()
    run_command(
        [str(LINUX / "fragmfonll")],
        run_dir,
        frag_input_points(grid_file.name, frag_mode, frag_param, points),
        ROOT / "logs" / f"fragmfonll_{pdf_key}_{quark}_{label}.log",
        env,
    )
    if not frag_path.exists():
        raise RuntimeError(f"expected fragmentation file not found: {frag_path}")
    saved_path = run_dir / f"fragmfonll_{label}.dat"
    shutil.copy2(frag_path, saved_path)
    return parse_frag(frag_path), saved_path


def ensure_charm_calibration_grid(env: dict[str, str]) -> Path:
    run_dir = RUNDIR / "validate_cteq66_charm_y0"
    run_dir.mkdir(parents=True, exist_ok=True)
    grid_file = run_dir / "val_c_.out"
    if grid_file.exists():
        return grid_file
    run_command(
        [str(LINUX / "fonllgridlha")],
        run_dir,
        grid_input("val_c_", 10550, 1.5, y_values=[-0.5, 0.0, 0.5]),
        ROOT / "logs" / "validate_cteq66_charm_y0_fonllgrid.log",
        env,
    )
    if not grid_file.exists():
        raise RuntimeError(f"expected calibration grid file not found: {grid_file}")
    return grid_file


def combine_charm_d0_feeddown(
    direct_values: dict[tuple[float, float], float],
    vector_values: dict[tuple[float, float], float],
    vector_to_direct_weight: float,
    points: Iterable[tuple[float, float]],
) -> dict[tuple[float, float], float]:
    channels = charm_d0_feeddown_channels()
    branch_sum = sum(float(channel["branching_fraction"]) for channel in channels)
    combined: dict[tuple[float, float], float] = {}
    for pt, y in points:
        feeddown = 0.0
        for channel in channels:
            br = float(channel["branching_fraction"])
            z = float(channel["z_collinear"])
            feeddown += br * vector_values[point_key(pt / z, y)] / z
        direct = direct_values[point_key(pt, y)]
        combined[point_key(pt, y)] = (
            direct + vector_to_direct_weight * feeddown
        ) / (1.0 + vector_to_direct_weight * branch_sum)
    return combined


def fit_charm_d0_feeddown_weight(env: dict[str, str]) -> dict[str, object]:
    grid_file = ensure_charm_calibration_grid(env)
    run_dir = grid_file.parent
    reference_points = [(pt, 0.0) for pt in sorted(PUBLIC_CHARM_D0_CTEQ66_Y0)]
    parent_points = charm_feeddown_parent_points(reference_points)
    direct_values, _ = run_fragmentation(
        "cteq66",
        "charm",
        run_dir,
        grid_file,
        5,
        0.1,
        reference_points,
        "calibration_direct_pseudoscalar",
        env,
    )
    vector_values, _ = run_fragmentation(
        "cteq66",
        "charm",
        run_dir,
        grid_file,
        4,
        0.1,
        parent_points,
        "calibration_dstar_vector",
        env,
    )

    def objective(weight: float) -> float:
        combined = combine_charm_d0_feeddown(
            direct_values,
            vector_values,
            weight,
            reference_points,
        )
        return sum(
            math.log(combined[point_key(pt, 0.0)] / public_value) ** 2
            for pt, public_value in PUBLIC_CHARM_D0_CTEQ66_Y0.items()
        )

    lo, hi = 0.0, 5.0
    for _ in range(100):
        left = lo + (hi - lo) / 3.0
        right = hi - (hi - lo) / 3.0
        if objective(left) < objective(right):
            hi = right
        else:
            lo = left
    weight = 0.5 * (lo + hi)
    combined = combine_charm_d0_feeddown(direct_values, vector_values, weight, reference_points)
    point_checks = []
    max_abs_rel = 0.0
    for pt, public_value in sorted(PUBLIC_CHARM_D0_CTEQ66_Y0.items()):
        local_value = combined[point_key(pt, 0.0)]
        rel = (local_value - public_value) / public_value
        max_abs_rel = max(max_abs_rel, abs(rel))
        point_checks.append(
            {
                "pt": pt,
                "y": 0.0,
                "local": local_value,
                "public": public_value,
                "relative_difference": rel,
            }
        )
    return {
        "grid": str(grid_file),
        "vector_to_direct_weight": weight,
        "max_abs_relative_difference": max_abs_rel,
        "points": point_checks,
    }


def load_or_fit_charm_feeddown_weight(
    env: dict[str, str], cache_path: Path
) -> dict[str, object]:
    """Calibrate the D* -> D0 feeddown weight once and cache it.

    The weight is anchored to fixed public CTEQ6.6 references, so it is a
    property of the fragmentation model, not of the scale/PDF/mass variation.
    Re-fitting per variation would conflate the perturbative variation with the
    fragmentation anchor (and triple the charm compute), so every charm
    variation reuses this single calibration.
    """
    if cache_path.exists():
        return json.loads(cache_path.read_text())
    calibration = fit_charm_d0_feeddown_weight(env)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(calibration, indent=2) + "\n")
    return calibration


def _short_pdf(pdf_key: str, pdf_entry: tuple) -> str:
    if len(pdf_entry) == 3:
        return pdf_entry[2]
    return "nnpdf40_nlo_as_01180" if pdf_key == "nlo" else "nnpdf40_nnlo_as_01180"


def generate(
    pdf_key: str,
    quark: str,
    variation: Variation,
    reuse_existing_grids: bool,
    grid_workers: int,
    feeddown_calibration: dict[str, object] | None = None,
    out_root: Path | None = None,
    run_root: Path | None = None,
) -> dict[str, object]:
    pdf_entry = PDFS[pdf_key]
    pdf_name, base_pdf_id = pdf_entry[:2]
    quark_config = QUARKS[quark]
    mass = variation.mass
    lhaid = int(base_pdf_id) + int(variation.pdf_member)
    ffact, fren = variation.muF, variation.muR  # FONLL: ffact=muF, fren=muR
    short_pdf = _short_pdf(pdf_key, pdf_entry)

    out_root = OUTDIR if out_root is None else out_root
    run_root = RUNDIR if run_root is None else run_root

    final_name = (
        f"fonll_pp14tev_{short_pdf}_fonll_meson_dsdpTdy_"
        f"pt0-50_y-3to3_{variation.tag}_{quark}.dat"
    )
    out_path = out_root / final_name
    run_dir = run_root / f"{pdf_key}_{quark}_{variation.tag}"
    run_dir.mkdir(parents=True, exist_ok=True)
    env = fonll_env()
    prefix = f"{pdf_key[:2]}_{quark[:1]}_"
    grid_file = run_dir / f"{prefix}.out"
    log_label = variation.tag
    stale_files = [
        run_dir / f"{prefix}.outlog",
        run_dir / f"{prefix}fonll.log",
        run_dir / f"{prefix}hdml.tmp",
        run_dir / f"{prefix}hdmv.tmp",
        run_dir / f"{prefix}hdrs.tmp",
        run_dir / "fragmfonll.dat",
        run_dir / "fragfonll.log",
    ]
    if not reuse_existing_grids:
        stale_files.append(grid_file)
    for stale in stale_files:
        if stale.exists():
            stale.unlink()

    if reuse_existing_grids and grid_file.exists():
        print(f"[{pdf_key} {quark} {variation.tag}] reusing existing quark grid {grid_file}")
    else:
        print(f"[{pdf_key} {quark} {variation.tag}] building quark grid")
        build_quark_grid(
            run_dir,
            grid_file,
            prefix,
            lhaid,
            mass,
            pdf_key,
            quark,
            env,
            grid_workers,
            ffact=ffact,
            fren=fren,
            log_label=log_label,
        )
    if not grid_file.exists():
        raise RuntimeError(f"expected grid file not found: {grid_file}")

    points = observed_points()
    extra_headers: dict[str, object] = {
        "lhapdf_member": variation.pdf_member,
        "variation_kind": variation.kind,
        "variation_tag": variation.tag,
    }
    weight = None
    if quark_config.get("feeddown") == "public_d0":
        if feeddown_calibration is None:
            feeddown_calibration = load_or_fit_charm_feeddown_weight(
                env, out_root / "charm_feeddown_calibration.json"
            )
        weight = float(feeddown_calibration["vector_to_direct_weight"])

        print(f"[{pdf_key} {quark} {variation.tag}] applying direct pseudoscalar charm fragmentation")
        direct_values, _ = run_fragmentation(
            pdf_key,
            f"{quark}_{variation.tag}",
            run_dir,
            grid_file,
            int(quark_config["frag_mode"]),
            float(quark_config["frag_param"]),
            points,
            "direct_pseudoscalar",
            env,
        )
        print(f"[{pdf_key} {quark} {variation.tag}] applying D* vector fragmentation for feeddown")
        vector_values, _ = run_fragmentation(
            pdf_key,
            f"{quark}_{variation.tag}",
            run_dir,
            grid_file,
            4,
            float(quark_config["frag_param"]),
            charm_feeddown_parent_points(points),
            "dstar_vector_feeddown",
            env,
        )
        values = combine_charm_d0_feeddown(direct_values, vector_values, weight, points)
        extra_headers.update(
            {
                "charm_feeddown_model": "collinear two-body D* -> D0 feeddown, calibrated to public FONLL CTEQ6.6 D0",
                "charm_feeddown_vector_to_direct_weight": f"{weight:.12g}",
                "charm_feeddown_branching_sum": f"{sum(float(c['branching_fraction']) for c in charm_d0_feeddown_channels()):.12g}",
                "charm_feeddown_cteq66_reference": "public FONLL v1.3.2 CTEQ6.6, queried 2026-05-29",
                "charm_feeddown_cteq66_fit_max_abs_rel_diff": f"{float(feeddown_calibration['max_abs_relative_difference']):.6e}",
            }
        )
    else:
        print(f"[{pdf_key} {quark} {variation.tag}] applying meson fragmentation")
        values, _ = run_fragmentation(
            pdf_key,
            f"{quark}_{variation.tag}",
            run_dir,
            grid_file,
            int(quark_config["frag_mode"]),
            float(quark_config["frag_param"]),
            points,
            "central",
            env,
        )

    summary = write_final(
        out_path,
        pdf_name,
        lhaid,
        quark,
        quark_config,
        values,
        grid_file,
        extra_headers=extra_headers,
        mass=mass,
        ffact=ffact,
        fren=fren,
    )
    summary.update(
        {
            "pdf_key": pdf_key,
            "quark": quark,
            "variation_kind": variation.kind,
            "variation_tag": variation.tag,
            "muR": variation.muR,
            "muF": variation.muF,
            "lhapdf_id": lhaid,
            "lhapdf_member": variation.pdf_member,
            "heavy_quark_mass_GeV": mass,
            "raw_grid_reused": reuse_existing_grids,
        }
    )
    if weight is not None:
        summary["charm_feeddown_vector_to_direct_weight"] = weight
    if not reuse_existing_grids:
        summary["raw_grid_workers"] = grid_workers
    return summary


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fonll_version() -> str:
    readme = FONLL / "README"
    if readme.exists():
        for line in readme.read_text(errors="ignore").splitlines():
            stripped = line.strip()
            if "version" in stripped.lower() or stripped.lower().startswith("fonll"):
                return stripped
    return "FONLL (version unrecorded; see src/fonll/README)"


def write_manifest(
    summaries: list[dict[str, object]],
    pdf_key: str,
    pdf_members: list[int],
    include_scale: bool,
    include_mass: bool,
    feeddown_calibration: dict[str, object] | None,
    out_root: Path,
) -> Path:
    pdf_entry = PDFS[pdf_key]
    pdf_name, base_pdf_id = pdf_entry[:2]
    entries = []
    for summary in sorted(summaries, key=lambda s: (s["quark"], s["variation_tag"])):
        path = Path(str(summary["path"]))
        entries.append(
            {
                "quark": summary["quark"],
                "variation_kind": summary["variation_kind"],
                "variation_tag": summary["variation_tag"],
                "muR": summary["muR"],
                "muF": summary["muF"],
                "lhapdf_id": summary["lhapdf_id"],
                "lhapdf_member": summary["lhapdf_member"],
                "heavy_quark_mass_GeV": summary["heavy_quark_mass_GeV"],
                "path": str(path),
                "sha256": sha256_file(path),
                "rows": summary["rows"],
                "trapezoid_integral_pb": summary["trapezoid_integral_pb"],
            }
        )
    manifest = {
        "generated_unix": time.time(),
        "generated_iso": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "fonll_version": fonll_version(),
        "patches": [
            "misc1/fragmfonll.f: BCFY vector/pseudoscalar modes 4 and 5/8",
            "misc1/fonllgrid.f: rapidity capacity raised to 120 nodes",
        ],
        "collision": "pp",
        "ebeam1_GeV": 7000,
        "ebeam2_GeV": 7000,
        "sqrt_s_GeV": 14000,
        "pdf_set": pdf_name,
        "pdf_base_lhapdf_id": base_pdf_id,
        "pdf_members": list(pdf_members),
        "scale_points_muR_muF": SCALE_POINTS_7 if include_scale else [(1.0, 1.0)],
        "mass_variations_GeV": MASS_VARIATIONS if include_mass else {},
        "grid_bounds": {
            "pt_min_GeV": PT_VALUES[0],
            "pt_max_GeV": PT_VALUES[-1],
            "y_min": Y_VALUES[0],
            "y_max": Y_VALUES[-1],
            "n_pt": len(PT_VALUES),
            "n_y": len(Y_VALUES),
        },
        "charm_feeddown": feeddown_calibration,
        "grids": entries,
    }
    manifest_path = out_root / "variation_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest_path


def parse_pdf_members(spec: str) -> list[int]:
    if spec.strip().lower() in {"", "none"}:
        return [0]
    members: list[int] = []
    for part in spec.split(","):
        part = part.strip()
        if "-" in part:
            lo, hi = part.split("-", 1)
            members.extend(range(int(lo), int(hi) + 1))
        else:
            members.append(int(part))
    # Always include the central member so the campaign produces a central grid.
    if 0 not in members:
        members.insert(0, 0)
    return sorted(dict.fromkeys(members))


def run_campaign(
    pdf_key: str,
    quarks: list[str],
    pdf_members: list[int],
    include_scale: bool,
    include_mass: bool,
    max_parallel: int,
    grid_workers: int,
    reuse_existing_grids: bool,
) -> Path:
    env = fonll_env()
    out_root = OUTDIR
    feeddown_calibration = None
    if "charm" in quarks:
        feeddown_calibration = load_or_fit_charm_feeddown_weight(
            env, out_root / "charm_feeddown_calibration.json"
        )

    tasks: list[tuple[str, str, Variation]] = []
    for quark in quarks:
        for variation in enumerate_variations(
            quark, pdf_members, include_scale, include_mass
        ):
            tasks.append((pdf_key, quark, variation))

    print(
        f"campaign: {len(tasks)} variation grids "
        f"({', '.join(quarks)}; scale={include_scale}; "
        f"pdf_members={len(pdf_members)}; mass={include_mass}); "
        f"max_parallel={max_parallel}, grid_workers={grid_workers}"
    )

    summaries: list[dict[str, object]] = []
    failures: list[tuple[tuple[str, str, Variation], str]] = []
    if max_parallel <= 1:
        for pk, q, var in tasks:
            summaries.append(
                generate(
                    pk, q, var, reuse_existing_grids, grid_workers,
                    feeddown_calibration=feeddown_calibration if q == "charm" else None,
                )
            )
    else:
        with ThreadPoolExecutor(max_workers=max_parallel) as executor:
            futures = {
                executor.submit(
                    generate,
                    pk,
                    q,
                    var,
                    reuse_existing_grids,
                    grid_workers,
                    feeddown_calibration if q == "charm" else None,
                ): (pk, q, var)
                for (pk, q, var) in tasks
            }
            for future in as_completed(futures):
                pk, q, var = futures[future]
                try:
                    summaries.append(future.result())
                except Exception as exc:  # noqa: BLE001 - report and continue
                    failures.append(((pk, q, var), str(exc)))
                    print(f"[FAILED] {pk} {q} {var.tag}: {exc}")

    manifest_path = write_manifest(
        summaries,
        pdf_key,
        pdf_members,
        include_scale,
        include_mass,
        feeddown_calibration,
        out_root,
    )
    print(f"wrote {manifest_path} ({len(summaries)} grids)")
    if failures:
        print(f"WARNING: {len(failures)} variation grids failed:")
        for (pk, q, var), msg in failures:
            print(f"  - {pk} {q} {var.tag}: {msg}")
    return manifest_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pdf",
        action="append",
        choices=sorted(PDFS),
        help="PDF key to generate; may be repeated. Default: nlo, plus nnlo with --include-nnlo.",
    )
    parser.add_argument("--include-nnlo", action="store_true", help="also generate NNLO comparison grids")
    parser.add_argument(
        "--quark",
        action="append",
        choices=sorted(QUARKS),
        help="quark flavor to generate; may be repeated. Default: bottom and charm.",
    )
    parser.add_argument("--reuse-existing-grids", action="store_true", help="skip fonllgridlha when the raw grid already exists")
    parser.add_argument("--grid-workers", type=int, default=1, help="parallel rapidity chunks for raw FONLL grid generation")
    parser.add_argument(
        "--campaign",
        action="store_true",
        help="run the production-uncertainty variation campaign (scale + PDF members + mass)",
    )
    parser.add_argument("--scale-variations", action="store_true", help="include the 7-point scale variation grids")
    parser.add_argument(
        "--pdf-members",
        default="0",
        help="PDF members to generate, e.g. '0-100' or '0,1,5'. Default: '0' (central only).",
    )
    parser.add_argument("--mass-variations", action="store_true", help="include heavy-quark-mass variation grids")
    parser.add_argument(
        "--max-parallel",
        type=int,
        default=1,
        help="number of variation grids to compute concurrently (bounds total cores with --grid-workers)",
    )
    parser.add_argument(
        "--compress-logs",
        action="store_true",
        help="gzip each per-command log on success (recommended for large campaigns)",
    )
    args = parser.parse_args()
    if args.grid_workers < 1:
        parser.error("--grid-workers must be at least 1")
    if args.max_parallel < 1:
        parser.error("--max-parallel must be at least 1")

    global _COMPRESS_LOGS
    _COMPRESS_LOGS = args.compress_logs

    requested = args.pdf or ["nlo"]
    if args.include_nnlo:
        requested.append("nnlo")
    requested = list(dict.fromkeys(requested))

    requested_quarks = args.quark or ["bottom", "charm"]

    if args.campaign:
        narrowed = args.scale_variations or args.mass_variations or args.pdf_members != "0"
        pdf_members = parse_pdf_members(args.pdf_members)
        if narrowed:
            # Honor exactly the axes the user asked for.
            include_scale = args.scale_variations
            include_mass = args.mass_variations
        else:
            # Bare --campaign: full default set (scale + all 100 members + mass).
            include_scale = True
            include_mass = True
            pdf_members = list(range(0, 101))
        for pdf_key in requested:
            run_campaign(
                pdf_key,
                requested_quarks,
                pdf_members,
                include_scale,
                include_mass,
                args.max_parallel,
                args.grid_workers,
                args.reuse_existing_grids,
            )
        return

    # Non-campaign path: central grid per quark (backward compatible).
    env = fonll_env()
    feeddown_calibration = None
    if "charm" in requested_quarks:
        feeddown_calibration = load_or_fit_charm_feeddown_weight(
            env, OUTDIR / "charm_feeddown_calibration.json"
        )
    summaries = []
    for pdf_key in requested:
        for quark in requested_quarks:
            central = Variation("central", "central", 1.0, 1.0, 0, float(QUARKS[quark]["mass"]))
            summaries.append(
                generate(
                    pdf_key,
                    quark,
                    central,
                    args.reuse_existing_grids,
                    args.grid_workers,
                    feeddown_calibration=feeddown_calibration if quark == "charm" else None,
                )
            )

    if requested == ["nlo"] and requested_quarks == ["bottom", "charm"]:
        summary_name = "generation_summary.json"
    else:
        summary_name = f"generation_summary_{'_'.join(requested + requested_quarks)}.json"
    summary_path = OUTDIR / summary_name
    summary_path.write_text(json.dumps(summaries, indent=2) + "\n")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
