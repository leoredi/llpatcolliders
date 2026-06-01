#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
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


def fonll_env() -> dict[str, str]:
    env = os.environ.copy()
    env["PATH"] = f"{ENV / 'bin'}:{env.get('PATH', '')}"
    env["LHAPDF_DATA_PATH"] = str(ENV / "share" / "LHAPDF")
    return env


def grid_input(
    prefix: str,
    pdf_id: int,
    mass: float,
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
        " 1. 1.",
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
) -> None:
    if grid_workers == 1:
        run_command(
            [str(LINUX / "fonllgridlha")],
            run_dir,
            grid_input(prefix, pdf_id, mass),
            ROOT / "logs" / f"fonllgrid_{pdf_key}_{quark}.log",
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
            grid_input(prefix, pdf_id, mass, y_values),
            ROOT / "logs" / f"fonllgrid_{pdf_key}_{quark}_chunk{index:02d}.log",
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
) -> dict[str, object]:
    out_path.parent.mkdir(parents=True, exist_ok=True)
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
        out.write(f"# heavy_quark_mass_GeV: {quark_config['mass']:.8g}\n")
        out.write("# scale: muR=muF=sqrt(m^2+pT^2), ffact=1, fren=1\n")
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


def generate(
    pdf_key: str,
    quark: str,
    reuse_existing_grids: bool,
    grid_workers: int,
) -> dict[str, object]:
    pdf_entry = PDFS[pdf_key]
    pdf_name, pdf_id = pdf_entry[:2]
    quark_config = QUARKS[quark]
    mass = quark_config["mass"]
    short_pdf = pdf_entry[2] if len(pdf_entry) == 3 else (
        "nnpdf40_nlo_as_01180" if pdf_key == "nlo" else "nnpdf40_nnlo_as_01180"
    )
    final_name = (
        f"fonll_pp14tev_{short_pdf}_fonll_meson_dsdpTdy_"
        f"pt0-50_y-3to3_central_{quark}.dat"
    )
    out_path = OUTDIR / final_name
    run_dir = RUNDIR / f"{pdf_key}_{quark}"
    run_dir.mkdir(parents=True, exist_ok=True)
    env = fonll_env()
    prefix = f"{pdf_key[:2]}_{quark[:1]}_"
    grid_file = run_dir / f"{prefix}.out"
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
        print(f"[{pdf_key} {quark}] reusing existing quark grid {grid_file}")
    else:
        print(f"[{pdf_key} {quark}] building quark grid")
        build_quark_grid(
            run_dir,
            grid_file,
            prefix,
            pdf_id,
            mass,
            pdf_key,
            quark,
            env,
            grid_workers,
        )
    if not grid_file.exists():
        raise RuntimeError(f"expected grid file not found: {grid_file}")

    points = observed_points()
    extra_headers: dict[str, object] = {}
    if quark_config.get("feeddown") == "public_d0":
        backup_path = OUTDIR / f"{out_path.stem}.direct_bcfy_pseudoscalar.dat"
        if out_path.exists() and not backup_path.exists():
            header = out_path.read_text(errors="ignore")[:2000]
            if "charm_feeddown_model" not in header:
                shutil.copy2(out_path, backup_path)

        print(f"[{pdf_key} {quark}] calibrating public-FONLL D0 feeddown")
        calibration = fit_charm_d0_feeddown_weight(env)
        weight = float(calibration["vector_to_direct_weight"])

        print(f"[{pdf_key} {quark}] applying direct pseudoscalar charm fragmentation")
        direct_values, direct_raw = run_fragmentation(
            pdf_key,
            quark,
            run_dir,
            grid_file,
            int(quark_config["frag_mode"]),
            float(quark_config["frag_param"]),
            points,
            "direct_pseudoscalar",
            env,
        )
        print(f"[{pdf_key} {quark}] applying D* vector fragmentation for feeddown")
        vector_values, vector_raw = run_fragmentation(
            pdf_key,
            quark,
            run_dir,
            grid_file,
            4,
            float(quark_config["frag_param"]),
            charm_feeddown_parent_points(points),
            "dstar_vector_feeddown",
            env,
        )
        values = combine_charm_d0_feeddown(direct_values, vector_values, weight, points)
        shutil.copy2(direct_raw, OUTDIR / f"{out_path.stem}.direct_pseudoscalar.fragmfonll_raw.dat")
        shutil.copy2(vector_raw, OUTDIR / f"{out_path.stem}.dstar_vector_feeddown.fragmfonll_raw.dat")
        stale_generic_raw = OUTDIR / f"{out_path.stem}.fragmfonll_raw.dat"
        if stale_generic_raw.exists():
            stale_generic_raw.unlink()
        extra_headers = {
            "charm_feeddown_model": "collinear two-body D* -> D0 feeddown, calibrated to public FONLL CTEQ6.6 D0",
            "charm_feeddown_vector_to_direct_weight": f"{weight:.12g}",
            "charm_feeddown_branching_sum": f"{sum(float(c['branching_fraction']) for c in charm_d0_feeddown_channels()):.12g}",
            "charm_feeddown_cteq66_reference": "public FONLL v1.3.2 CTEQ6.6, queried 2026-05-29",
            "charm_feeddown_cteq66_fit_max_abs_rel_diff": f"{float(calibration['max_abs_relative_difference']):.6e}",
        }
    else:
        print(f"[{pdf_key} {quark}] applying meson fragmentation")
        values, frag_path = run_fragmentation(
            pdf_key,
            quark,
            run_dir,
            grid_file,
            int(quark_config["frag_mode"]),
            float(quark_config["frag_param"]),
            points,
            "central",
            env,
        )
        shutil.copy2(frag_path, OUTDIR / f"{out_path.stem}.fragmfonll_raw.dat")

    summary = write_final(
        out_path,
        pdf_name,
        pdf_id,
        quark,
        quark_config,
        values,
        grid_file,
        extra_headers=extra_headers,
    )
    if extra_headers:
        summary.update(extra_headers)
    summary["raw_grid_reused"] = reuse_existing_grids
    if not reuse_existing_grids:
        summary["raw_grid_workers"] = grid_workers
    return summary


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
    args = parser.parse_args()
    if args.grid_workers < 1:
        parser.error("--grid-workers must be at least 1")

    requested = args.pdf or ["nlo"]
    if args.include_nnlo:
        requested.append("nnlo")
    requested = list(dict.fromkeys(requested))

    requested_quarks = args.quark or ["bottom", "charm"]

    summaries = []
    for pdf_key in requested:
        for quark in requested_quarks:
            summaries.append(generate(pdf_key, quark, args.reuse_existing_grids, args.grid_workers))

    if requested == ["nlo"] and requested_quarks == ["bottom", "charm"]:
        summary_name = "generation_summary.json"
    else:
        summary_name = f"generation_summary_{'_'.join(requested + requested_quarks)}.json"
    summary_path = OUTDIR / summary_name
    summary_path.write_text(json.dumps(summaries, indent=2) + "\n")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
