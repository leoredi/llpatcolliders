#!/usr/bin/env python3
"""Combine FONLL variation grids into scale/PDF/mass/combined uncertainty bands.

Reads ``output/variation_manifest.json`` (written by ``generate_meson_grids.py
--campaign``), partitions each quark's grids by variation axis, and writes
envelope grids in the same three-column ``pT y dsigma/dpT/dy`` format the
downstream HNL sampler consumes. Band definitions:

* scale: pointwise max/min over the 7-point (muR, muF) set;
* pdf:   NNPDF Monte-Carlo prescription -- central +/- std dev (ddof=1) over
         the replica members only (member 0 is the replica mean, not a
         replica, and is excluded from the statistics); the replica mean is
         also recorded;
* mass:  pointwise max/min over the m_b/m_c up/down set;
* combined: central +/- quadrature of the three (uncorrelated) deviations.

The scale and mass bands include the central grid in their max/min, so a
quark needs the central grid plus one non-central point on an axis for that
axis's band to be emitted; the pdf band needs at least two replica members.

These envelopes are pointwise constructions, not physical cross sections:
they carry an ``envelope_band`` header and must never be fed to the HNL
sampler (`hnl/production/fonll/fonll_parser.py` rejects them). Sample only
coherent individual grids (one scale point, one replica, one mass).
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "output"
ENVDIR = OUTDIR / "envelopes"


def _load_generator_module():
    """Import generate_meson_grids.py for shared node definitions and helpers."""
    path = Path(__file__).resolve().parent / "generate_meson_grids.py"
    spec = importlib.util.spec_from_file_location("generate_meson_grids", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


GMG = _load_generator_module()
PT_VALUES = GMG.PT_VALUES
Y_VALUES = GMG.Y_VALUES
QUARKS = GMG.QUARKS
ATOL = 1e-9


def _close(a: float, b: float) -> bool:
    return abs(float(a) - float(b)) < ATOL


def load_manifest(out_root: Path) -> dict:
    manifest_path = out_root / "variation_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"{manifest_path} not found; run generate_meson_grids.py --campaign first"
        )
    return json.loads(manifest_path.read_text())


def reference_grid_coords() -> tuple[np.ndarray, np.ndarray]:
    """Canonical (pT, y) column order: pT outermost, y innermost, both ascending."""
    pt = np.repeat(np.array(PT_VALUES, dtype=float), len(Y_VALUES))
    y = np.tile(np.array(Y_VALUES, dtype=float), len(PT_VALUES))
    return pt, y


def load_grid_column(path: Path, ref_pt: np.ndarray, ref_y: np.ndarray) -> np.ndarray:
    """Load the dsigma column, asserting the grid matches the canonical ordering."""
    data = np.loadtxt(path, comments="#")
    if data.shape != (len(ref_pt), 3):
        raise ValueError(f"{path}: expected {len(ref_pt)}x3 grid, got {data.shape}")
    if not (np.allclose(data[:, 0], ref_pt, atol=1e-4)
            and np.allclose(data[:, 1], ref_y, atol=1e-4)):
        raise ValueError(f"{path}: grid coordinates do not match the canonical pT/y ordering")
    col = data[:, 2]
    if not np.all(np.isfinite(col)):
        raise ValueError(f"{path}: non-finite dsigma values")
    return col


def partition_entries(entries: list[dict], quark: str) -> dict[str, object]:
    default_mass = float(QUARKS[quark]["mass"])
    q_entries = [e for e in entries if e["quark"] == quark]
    central = None
    for e in q_entries:
        if (e["lhapdf_member"] == 0 and _close(e["muR"], 1.0) and _close(e["muF"], 1.0)
                and _close(e["heavy_quark_mass_GeV"], default_mass)):
            central = e
            break
    if central is None:
        raise ValueError(f"{quark}: no central grid (member 0, scale 1,1, mass {default_mass}) in manifest")

    scale_set = [
        e for e in q_entries
        if e["lhapdf_member"] == 0 and _close(e["heavy_quark_mass_GeV"], default_mass)
    ]
    pdf_set = [
        e for e in q_entries
        if _close(e["muR"], 1.0) and _close(e["muF"], 1.0)
        and _close(e["heavy_quark_mass_GeV"], default_mass)
    ]
    mass_set = [
        e for e in q_entries
        if e["lhapdf_member"] == 0 and _close(e["muR"], 1.0) and _close(e["muF"], 1.0)
    ]
    return {
        "central": central,
        "scale": sorted(scale_set, key=lambda e: (e["muR"], e["muF"])),
        "pdf": sorted(pdf_set, key=lambda e: e["lhapdf_member"]),
        "mass": sorted(mass_set, key=lambda e: e["heavy_quark_mass_GeV"]),
        "default_mass": default_mass,
    }


def write_envelope_grid(
    out_path: Path,
    quark: str,
    band: str,
    column: np.ndarray,
    ref_pt: np.ndarray,
    ref_y: np.ndarray,
    provenance: dict,
) -> dict:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    column = np.maximum(column, 0.0)
    dpt = PT_VALUES[1] - PT_VALUES[0]
    dy = Y_VALUES[1] - Y_VALUES[0]
    integral = float(np.trapz(
        np.trapz(column.reshape(len(PT_VALUES), len(Y_VALUES)), Y_VALUES, axis=1),
        PT_VALUES,
    ))
    with out_path.open("w") as out:
        out.write("# FONLL heavy-flavor meson grid (variation envelope)\n")
        out.write("# columns: pT y dsigma/dpT/dy\n")
        out.write("# units: GeV 1 pb/GeV\n")
        out.write("# collision: pp\n")
        out.write("# sqrt_s_GeV: 14000\n")
        out.write(f"# quark: {quark}\n")
        out.write(f"# envelope_band: {band}\n")
        for key, value in provenance.items():
            out.write(f"# {key}: {value}\n")
        out.write(f"# trapezoid_integral_pb_y-3to3_pt0to50: {integral:.12e}\n")
        out.write("# pT y dsigma/dpT/dy\n")
        for pt, y, val in zip(ref_pt, ref_y, column):
            out.write(f"{pt:.8g} {y:.8g} {val:.12e}\n")
    return {
        "band": band,
        "path": str(out_path),
        "sha256": GMG.sha256_file(out_path),
        "trapezoid_integral_pb": integral,
    }


def combine_quark(quark: str, manifest: dict, out_root: Path) -> list[dict]:
    entries = manifest["grids"]
    parts = partition_entries(entries, quark)
    ref_pt, ref_y = reference_grid_coords()

    central_col = load_grid_column(Path(parts["central"]["path"]), ref_pt, ref_y)
    pdf_name = manifest.get("pdf_set", "")
    stem = (
        f"fonll_pp14tev_{pdf_name.lower().replace('.', '')}_fonll_meson_dsdpTdy_"
        f"pt0-50_y-3to3_envelope"
    )
    written: list[dict] = []

    # Always emit the central reference for downstream convenience.
    written.append(write_envelope_grid(
        out_root / f"{stem}_central_{quark}.dat", quark, "central", central_col,
        ref_pt, ref_y,
        {"source": parts["central"]["path"], "n_grids": 1},
    ))

    # Scale band: pointwise max/min over the (muR, muF) set (includes central).
    scaleup = scaledn = None
    if len(parts["scale"]) >= 2:
        cols = np.stack([load_grid_column(Path(e["path"]), ref_pt, ref_y) for e in parts["scale"]])
        scaleup = cols.max(axis=0)
        scaledn = cols.min(axis=0)
        prov = {"n_grids": len(parts["scale"]),
                "points": "; ".join(f"({e['muR']},{e['muF']})" for e in parts["scale"])}
        written.append(write_envelope_grid(
            out_root / f"{stem}_scaleup_{quark}.dat", quark, "scale_up", scaleup, ref_pt, ref_y, prov))
        written.append(write_envelope_grid(
            out_root / f"{stem}_scaledn_{quark}.dat", quark, "scale_dn", scaledn, ref_pt, ref_y, prov))
    else:
        print(f"[{quark}] scale band skipped (only {len(parts['scale'])} grid)")

    # PDF band: NNPDF Monte-Carlo mean +/- std dev over the replica members.
    # Member 0 is the replica average, not a replica; including it in the
    # statistics would bias the std low, so it is excluded.
    pdf_sigma = None
    replicas = [e for e in parts["pdf"] if int(e["lhapdf_member"]) >= 1]
    if len(replicas) >= 2:
        cols = np.stack([load_grid_column(Path(e["path"]), ref_pt, ref_y) for e in replicas])
        pdf_mean = cols.mean(axis=0)
        pdf_sigma = cols.std(axis=0, ddof=1)
        prov = {"n_replicas": len(replicas),
                "prescription": ("NNPDF Monte-Carlo: central +/- std(ddof=1) over replica "
                                 "members >= 1 (member 0 excluded); replica mean also written")}
        written.append(write_envelope_grid(
            out_root / f"{stem}_pdfmean_{quark}.dat", quark, "pdf_mean", pdf_mean, ref_pt, ref_y, prov))
        written.append(write_envelope_grid(
            out_root / f"{stem}_pdfup_{quark}.dat", quark, "pdf_up", central_col + pdf_sigma, ref_pt, ref_y, prov))
        written.append(write_envelope_grid(
            out_root / f"{stem}_pdfdn_{quark}.dat", quark, "pdf_dn", central_col - pdf_sigma, ref_pt, ref_y, prov))
    else:
        print(f"[{quark}] pdf band skipped (only {len(replicas)} replica members)")

    # Mass band: pointwise max/min over the heavy-quark-mass set (includes central).
    massup = massdn = None
    if len(parts["mass"]) >= 2:
        cols = np.stack([load_grid_column(Path(e["path"]), ref_pt, ref_y) for e in parts["mass"]])
        massup = cols.max(axis=0)
        massdn = cols.min(axis=0)
        prov = {"n_grids": len(parts["mass"]),
                "masses_GeV": "; ".join(f"{e['heavy_quark_mass_GeV']:g}" for e in parts["mass"])}
        written.append(write_envelope_grid(
            out_root / f"{stem}_massup_{quark}.dat", quark, "mass_up", massup, ref_pt, ref_y, prov))
        written.append(write_envelope_grid(
            out_root / f"{stem}_massdn_{quark}.dat", quark, "mass_dn", massdn, ref_pt, ref_y, prov))
    else:
        print(f"[{quark}] mass band skipped (only {len(parts['mass'])} grid)")

    # Combined band: central +/- quadrature of the three uncorrelated deviations.
    components = []
    up_var = np.zeros_like(central_col)
    dn_var = np.zeros_like(central_col)
    if scaleup is not None:
        up_var += (scaleup - central_col) ** 2
        dn_var += (central_col - scaledn) ** 2
        components.append("scale")
    if pdf_sigma is not None:
        up_var += pdf_sigma ** 2
        dn_var += pdf_sigma ** 2
        components.append("pdf")
    if massup is not None:
        up_var += (massup - central_col) ** 2
        dn_var += (central_col - massdn) ** 2
        components.append("mass")
    if components:
        prov = {"components": "+".join(components),
                "treatment": "uncorrelated quadrature of per-axis deviations from central"}
        written.append(write_envelope_grid(
            out_root / f"{stem}_combup_{quark}.dat", quark, "combined_up",
            central_col + np.sqrt(up_var), ref_pt, ref_y, prov))
        written.append(write_envelope_grid(
            out_root / f"{stem}_combdn_{quark}.dat", quark, "combined_dn",
            central_col - np.sqrt(dn_var), ref_pt, ref_y, prov))
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quark", action="append", choices=sorted(QUARKS),
        help="quark to combine; may be repeated. Default: bottom and charm.",
    )
    args = parser.parse_args()
    quarks = args.quark or ["bottom", "charm"]

    manifest = load_manifest(OUTDIR)
    available = {e["quark"] for e in manifest["grids"]}
    ENVDIR.mkdir(parents=True, exist_ok=True)

    envelope_summary = {}
    for quark in quarks:
        if quark not in available:
            print(f"[{quark}] no grids in manifest; skipping")
            continue
        envelope_summary[quark] = combine_quark(quark, manifest, ENVDIR)

    summary = {
        "generated_iso": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "source_manifest": str(OUTDIR / "variation_manifest.json"),
        "envelopes": envelope_summary,
    }
    summary_path = ENVDIR / "envelope_manifest.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"wrote {summary_path}")
    for quark, grids in envelope_summary.items():
        print(f"  {quark}: {len(grids)} envelope grids -> {', '.join(g['band'] for g in grids)}")


if __name__ == "__main__":
    main()
