#!/usr/bin/env python3
"""Run and combine the independent full-statistics BC4 uncertainty campaign.

Every physics variation receives its own freshly sampled FONLL parent pool,
scalar four-vectors, detector ray casts, decay/reconstruction Monte Carlo, and
coupling scan.  No events, weights, geometry, or reconstruction outcomes are
shared between variations.

The full campaign consists of:

* one central FONLL grid;
* six non-central members of the coherent seven-point scale set;
* 100 NNPDF4.0 NLO Monte-Carlo PDF replicas;
* two bottom-mass grids (4.5 and 5.0 GeV); and
* one independently simulated LO-ChPT/spectator decay-model alternate using a
  fresh central-FONLL production sample.

Heavy products are written outside the repository.  The default can be
overridden either by ``--scratch-dir`` or ``BC4_UNCERTAINTY_DIR``::

    /Volumes/GRENDEL/extra_space/bc4_uncertainty/
      runs/<variation>/llp_4vectors/*.csv
      runs/<variation>/geometry_cache/*.npz
      runs/<variation>/results/*.json

Each four-vector file, mass result, and aggregate curve is committed atomically.
Rerunning ``run`` skips completed results and resumes incomplete variations.
Only compact aggregate curves, the combined band, and provenance belong under
``scalar/data/published/bundle``.

Examples
--------
    python -m scalar.uncertainty_band status
    python -m scalar.uncertainty_band run --workers 2
    python -m scalar.uncertainty_band collect
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import platform
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

_SCALAR_ROOT = Path(__file__).resolve().parent
_REPO_ROOT = _SCALAR_ROOT.parent
for _p in (str(_REPO_ROOT), str(_REPO_ROOT / "hnl"), str(_REPO_ROOT / "higgs")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from grendel_geometry import mesh_fiducial                            # noqa: E402
from production.fonll.fonll_parser import get_sigma_total             # noqa: E402
from production.fonll.meson_sampler import sample_meson_4vectors      # noqa: E402

from scalar import acceptance, production                             # noqa: E402


LOG_S2T_MIN, LOG_S2T_MAX, N_S2T = -12.0, -2.0, 200
CENTRAL_SCHEME = "winkler"
DECAY_SCHEME = "chpt_spectator"
DECAY_VARIATION = "decay_chpt_spectator"
GRID_ENV = "HNL_FONLL_BOTTOM_GRID"
GRID_STEM = "fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3"

DEFAULT_SCRATCH = Path(os.environ.get(
    "BC4_UNCERTAINTY_DIR", "/Volumes/GRENDEL/extra_space/bc4_uncertainty"))
DEFAULT_GRID_DIR = (_REPO_ROOT.parents[2] / "shared" / "NNPDF40"
                    / "fonll-local" / "output")
PUBLISHED_DIR = _SCALAR_ROOT / "data" / "published" / "bundle"
PUBLISHED_CENTRAL = _SCALAR_ROOT / "data" / "published" / "bc4_island.csv"
PUBLISHED_CURVES = PUBLISHED_DIR / "bc4_uncertainty_curves.csv"
PUBLISHED_BAND = PUBLISHED_DIR / "bc4_uncertainty_band.csv"
PUBLISHED_MANIFEST = PUBLISHED_DIR / "bc4_uncertainty_variations.json"

BOUNDARIES = (("u2_min", "u2_min_open"), ("u2_max", "u2_max_open"))
CODE_INPUTS = (
    "scalar/model.py",
    "scalar/production.py",
    "scalar/acceptance.py",
    "scalar/uncertainty_band.py",
    "hnl/analysis/decay_reco_acceptance.py",
    "hnl/analysis/exclusion.py",
    "hnl/analysis/_engine.py",
    "higgs/reco_common.py",
    "higgs/grendel_geometry.py",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, allow_nan=True)
        fh.write("\n")
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def _atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp.{os.getpid()}")
    frame.to_csv(tmp, index=False)
    with open(tmp, "rb") as fh:
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def _jsonable(value):
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def _grid_path(grid_dir: Path, tag: str) -> Path:
    return grid_dir / f"{GRID_STEM}_{tag}_bottom.dat"


def _code_hashes() -> dict[str, str]:
    return {rel: _sha256(_REPO_ROOT / rel) for rel in CODE_INPUTS}


def _git_head() -> str:
    proc = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=_REPO_ROOT,
        text=True, capture_output=True, check=False)
    return proc.stdout.strip() if proc.returncode == 0 else "unknown"


def _stable_seed(base_seed: int, *parts) -> int:
    message = ":".join([str(base_seed), *map(str, parts)]).encode()
    return int.from_bytes(hashlib.sha256(message).digest()[:8], "big") % (2**32)


def discover_variations(grid_dir: Path, validate_hashes=True):
    """Discover and validate the complete coherent bottom-grid campaign."""
    grid_dir = Path(grid_dir)
    manifest_path = grid_dir / "variation_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    bottom = [g for g in manifest["grids"] if g["quark"] == "bottom"]
    central_entry = next(
        (entry for entry in bottom if entry["variation_tag"] == "central"), None)
    if central_entry is None:
        raise ValueError("FONLL variation manifest has no central bottom grid")
    central_mass = float(central_entry["heavy_quark_mass_GeV"])

    records = []
    for entry in bottom:
        tag = entry["variation_tag"]
        kind = entry["variation_kind"]
        if tag == "central":
            name, axis = "central", "central"
        elif kind == "scale":
            name, axis = tag, "scale"
        elif kind == "pdf":
            if int(entry.get("lhapdf_member", 0)) == 0:
                continue
            name, axis = tag, "pdf"
        elif kind == "mass":
            mass = float(entry["heavy_quark_mass_GeV"])
            name, axis = ("mb_dn" if mass < central_mass else "mb_up"), "mass"
        else:
            continue

        path = _grid_path(grid_dir, tag)
        if not path.exists():
            raise FileNotFoundError(path)
        actual_sha = _sha256(path)
        if validate_hashes and actual_sha != entry["sha256"]:
            raise ValueError(
                f"checksum mismatch for {path.name}: {actual_sha} != {entry['sha256']}")
        records.append({
            "name": name,
            "axis": axis,
            "width_scheme": CENTRAL_SCHEME,
            "variation_tag": tag,
            "path": path,
            "sha256": actual_sha,
            "trapezoid_integral_pb": float(entry["trapezoid_integral_pb"]),
            "muR": entry.get("muR"),
            "muF": entry.get("muF"),
            "lhapdf_member": entry.get("lhapdf_member"),
            "heavy_quark_mass_GeV": entry.get("heavy_quark_mass_GeV"),
        })

    counts = {axis: sum(v["axis"] == axis for v in records)
              for axis in ("central", "scale", "pdf", "mass")}
    expected = {"central": 1, "scale": 6, "pdf": 100, "mass": 2}
    if counts != expected:
        raise ValueError(f"incomplete FONLL campaign: found {counts}, expected {expected}")

    order = {"central": 0, "scale": 2, "pdf": 3, "mass": 4}
    records.sort(key=lambda record: (order[record["axis"]], record["name"]))
    central = records[0]
    decay = {
        **central,
        "name": DECAY_VARIATION,
        "axis": "decay_model",
        "width_scheme": DECAY_SCHEME,
    }
    # Run the two independent central-grid simulations first so the decay-model
    # comparison completes before the long PDF campaign.
    records = [central, decay, *records[1:]]
    return manifest, records


def _select_variations(variations, selected):
    if not selected:
        return variations
    wanted = set(selected)
    selected_records = [record for record in variations
                        if record["name"] in wanted or record["axis"] in wanted]
    if not selected_records:
        raise ValueError(f"variation selection matched nothing: {sorted(wanted)}")
    return selected_records


@contextmanager
def _bottom_grid(path: Path):
    old = os.environ.get(GRID_ENV)
    os.environ[GRID_ENV] = str(path)
    try:
        yield
    finally:
        if old is None:
            os.environ.pop(GRID_ENV, None)
        else:
            os.environ[GRID_ENV] = old


@contextmanager
def _variation_lock(run_dir: Path):
    run_dir.mkdir(parents=True, exist_ok=True)
    lock_path = run_dir / ".lock"
    with open(lock_path, "a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(f"variation already running: {run_dir.name}") from exc
        lock.seek(0)
        lock.truncate()
        lock.write(f"pid={os.getpid()} host={platform.node()}\n")
        lock.flush()
        try:
            yield
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def _campaign_config(variation, masses, n_pool, n_samples, seed):
    payload = {
        "schema_version": 1,
        "variation": variation["name"],
        "axis": variation["axis"],
        "width_scheme": variation["width_scheme"],
        "grid_file": variation["path"].name,
        "grid_sha256": variation["sha256"],
        "n_parent_pool": int(n_pool),
        "n_scalar_events_expected": int(n_pool * len(production.B_SPECIES)),
        "n_decay_samples_per_hit": int(n_samples),
        "base_seed": int(seed),
        "parent_pool_seed": _stable_seed(seed, variation["name"], "parent_pool"),
        "mass_grid_GeV": [float(mass) for mass in masses],
        "coupling_grid": {
            "log10_min": LOG_S2T_MIN,
            "log10_max": LOG_S2T_MAX,
            "n_points": N_S2T,
        },
        "code_sha256": _code_hashes(),
    }
    stable = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["config_sha256"] = hashlib.sha256(stable).hexdigest()
    return payload


def _ensure_run_metadata(run_dir: Path, config) -> None:
    path = run_dir / "run_metadata.json"
    if path.exists():
        old = json.loads(path.read_text())
        if old.get("config_sha256") != config["config_sha256"]:
            raise RuntimeError(
                f"campaign configuration changed for {run_dir.name}; use a new "
                "--scratch-dir or restore the recorded inputs")
        return
    payload = {
        **config,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "producer_git_head_at_start": _git_head(),
        "independence_policy": (
            "fresh variation-specific FONLL parent pool, scalar four-vectors, "
            "geometry, decay/reconstruction MC, and sensitivity scan; no reuse "
            "or reweighting across variations"
        ),
    }
    _atomic_json(path, payload)


def _vector_paths(run_dir: Path, mass: float):
    label = production._mass_label(mass)
    vector = run_dir / "llp_4vectors" / f"mS_{label}.csv"
    return vector, vector.with_suffix(".meta.json")


def _valid_vector(vector: Path, metadata: Path, config, mass: float):
    if not vector.exists() or not metadata.exists():
        return None
    try:
        meta = json.loads(metadata.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if (meta.get("config_sha256") != config["config_sha256"]
            or float(meta.get("mass_GeV", -1)) != float(mass)
            or meta.get("bytes") != vector.stat().st_size):
        return None
    if meta.get("sha256") != _sha256(vector):
        return None
    return meta


def _write_vector(vector: Path, generated, config, variation, mass, seed):
    weight, energy, px, py, pz = generated
    matrix = np.column_stack([weight, energy, px, py, pz])
    vector.parent.mkdir(parents=True, exist_ok=True)
    tmp = vector.with_suffix(vector.suffix + f".tmp.{os.getpid()}")
    with open(tmp, "w") as fh:
        np.savetxt(fh, matrix, delimiter=",", fmt="%.8e")
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, vector)
    meta = {
        "config_sha256": config["config_sha256"],
        "variation": variation["name"],
        "grid_sha256": variation["sha256"],
        "mass_GeV": float(mass),
        "production_seed": int(seed),
        "rows": int(len(matrix)),
        "bytes": vector.stat().st_size,
        "sha256": _sha256(vector),
    }
    _atomic_json(vector.with_suffix(".meta.json"), meta)
    return meta


def _result_path(run_dir: Path, mass: float):
    return run_dir / "results" / f"mS_{production._mass_label(mass)}.json"


def _valid_result(path: Path, config):
    if not path.exists():
        return None
    try:
        result = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    return result if result.get("config_sha256") == config["config_sha256"] else None


def run_variation(variation, scratch_dir: Path, masses, n_pool, n_samples, seed):
    """Run one fully independent physics variation, resuming per mass."""
    scratch_dir = Path(scratch_dir)
    run_dir = scratch_dir / "runs" / variation["name"]
    config = _campaign_config(variation, production.MASS_GRID, n_pool, n_samples, seed)
    with _variation_lock(run_dir):
        _ensure_run_metadata(run_dir, config)
        complete = sum(_valid_result(_result_path(run_dir, mass), config) is not None
                       for mass in masses)
        print(f"[{variation['name']}] resume {complete}/{len(masses)} masses", flush=True)
        if complete == len(masses):
            return variation["name"]

        pool_seed = config["parent_pool_seed"]
        with _bottom_grid(variation["path"]):
            sigma = get_sigma_total("bottom")
            pool = sample_meson_4vectors(
                n_pool, "bottom", rng=np.random.default_rng(pool_seed))
        expected_sigma = variation["trapezoid_integral_pb"]
        # The manifest records a rounded decimal integral, while the parser
        # integrates the stored grid values at full float precision.
        if not np.isclose(sigma, expected_sigma, rtol=1e-8, atol=0.0):
            raise RuntimeError(
                f"{variation['name']}: grid integral {sigma} != manifest {expected_sigma}")

        s2_grid = np.logspace(LOG_S2T_MIN, LOG_S2T_MAX, N_S2T)
        for index, mass in enumerate(masses):
            mass = float(mass)
            result_path = _result_path(run_dir, mass)
            if _valid_result(result_path, config) is not None:
                continue

            production_seed = _stable_seed(
                seed, variation["name"], production._mass_label(mass), "production")
            reco_seed = _stable_seed(
                seed, variation["name"], production._mass_label(mass), "reconstruction")
            vector, vector_meta_path = _vector_paths(run_dir, mass)
            vector_meta = _valid_vector(vector, vector_meta_path, config, mass)
            if vector_meta is None:
                generated = production.generate_scalar_4vectors(
                    mass, n_pool, np.random.default_rng(production_seed),
                    sigma_bottom=sigma, pool=pool)
                vector_meta = _write_vector(
                    vector, generated, config, variation, mass, production_seed)
            if vector_meta["rows"] != config["n_scalar_events_expected"]:
                raise RuntimeError(
                    f"{variation['name']} m={mass}: expected "
                    f"{config['n_scalar_events_expected']} scalar events, got "
                    f"{vector_meta['rows']}")

            output = acceptance.process_mass_point(
                mass, mesh_fiducial, vector, s2_grid, n_samples=n_samples,
                rng=np.random.default_rng(reco_seed),
                width_scheme=variation["width_scheme"])
            if output is None:
                raise RuntimeError(f"{variation['name']} m={mass}: no production output")
            band = output[0] if isinstance(output, tuple) else output
            record = {
                **_jsonable(band),
                "variation": variation["name"],
                "axis": variation["axis"],
                "width_scheme": variation["width_scheme"],
                "mass_GeV": mass,
                "config_sha256": config["config_sha256"],
                "grid_sha256": variation["sha256"],
                "vector_sha256": vector_meta["sha256"],
                "parent_pool_seed": int(pool_seed),
                "production_seed": int(production_seed),
                "reconstruction_seed": int(reco_seed),
                "n_parent_pool": int(n_pool),
                "n_decay_samples_per_hit": int(n_samples),
            }
            _atomic_json(result_path, record)
            print(
                f"[{variation['name']}] {index + 1}/{len(masses)} m={mass:.3f} "
                f"hits={record['n_hits']} peak={record['peak_N']:.3g} "
                f"island=[{record['u2_min']}, {record['u2_max']}]", flush=True)

        rows = [json.loads(_result_path(run_dir, mass).read_text()) for mass in masses]
        curve = pd.DataFrame(rows).sort_values("mass_GeV")
        curve_path = run_dir / "curve.csv"
        _atomic_csv(curve, curve_path)
        completion = {
            "variation": variation["name"],
            "config_sha256": config["config_sha256"],
            "completed_utc": datetime.now(timezone.utc).isoformat(),
            "n_masses": len(curve),
            "curve_file": curve_path.name,
            "curve_sha256": _sha256(curve_path),
        }
        _atomic_json(run_dir / "complete.json", completion)
        return variation["name"]


def _boundary(row, boundary, open_col):
    if row is None or not bool(row.get("has_sensitivity", False)):
        return None, False
    is_open = bool(row.get(open_col, False))
    value = row.get(boundary, np.nan)
    if not np.isfinite(value) or value <= 0:
        return None, is_open
    return float(np.log10(value)), is_open


def combine_band(raw, central_curve):
    """Combine exact variation curves and rebase shifts onto the published central."""
    raw = pd.read_csv(raw) if isinstance(raw, (str, Path)) else raw.copy()
    reference = (pd.read_csv(central_curve) if isinstance(central_curve, (str, Path))
                 else central_curve.copy())
    axes = {axis: sorted(raw.loc[raw["axis"] == axis, "variation"].unique())
            for axis in ("scale", "pdf", "mass", "decay_model")}
    rows = []
    for _, ref in reference.sort_values("mass_GeV").iterrows():
        mass = float(ref["mass_GeV"])
        group = raw[np.isclose(raw["mass_GeV"], mass, rtol=0, atol=5e-10)]
        by_name = {row["variation"]: row for _, row in group.iterrows()}
        campaign_central = by_name.get("central")
        rec = {
            "mass_GeV": mass,
            "has_sensitivity": bool(ref.get("has_sensitivity", False)),
            "campaign_has_sensitivity": bool(
                campaign_central is not None
                and campaign_central.get("has_sensitivity", False)),
        }
        for boundary, open_col in BOUNDARIES:
            ref_value = ref.get(boundary, np.nan)
            ref_open = bool(ref.get(open_col, False))
            rec[f"{boundary}_central"] = ref_value
            rec[f"{boundary}_open"] = ref_open
            xc, campaign_open = _boundary(campaign_central, boundary, open_col)
            if (not rec["has_sensitivity"] or ref_open or not np.isfinite(ref_value)
                    or ref_value <= 0 or xc is None):
                rec[f"{boundary}_band_lo"] = np.nan
                rec[f"{boundary}_band_hi"] = np.nan
                rec[f"{boundary}_open"] = ref_open or campaign_open
                continue

            xref = float(np.log10(ref_value))
            any_open = ref_open or campaign_open
            missing = False

            scale_values = [xc]
            for name in axes["scale"]:
                value, is_open = _boundary(by_name.get(name), boundary, open_col)
                any_open |= is_open
                missing |= value is None and not is_open
                if value is not None:
                    scale_values.append(value)
            scale_up = max(scale_values) - xc
            scale_dn = xc - min(scale_values)

            pdf_values = []
            for name in axes["pdf"]:
                value, is_open = _boundary(by_name.get(name), boundary, open_col)
                any_open |= is_open
                missing |= value is None and not is_open
                if value is not None:
                    pdf_values.append(value)
            pdf_sigma = (float(np.std(pdf_values, ddof=1))
                         if len(pdf_values) >= 2 else 0.0)

            mass_dev = 0.0
            for name in axes["mass"]:
                value, is_open = _boundary(by_name.get(name), boundary, open_col)
                any_open |= is_open
                missing |= value is None and not is_open
                if value is not None:
                    mass_dev = max(mass_dev, abs(value - xc))

            alt, alt_open = _boundary(
                by_name.get(DECAY_VARIATION), boundary, open_col)
            any_open |= alt_open
            missing |= alt is None and not alt_open
            model_shift = 0.0 if alt is None else alt - xc
            model_up, model_dn = max(model_shift, 0.0), max(-model_shift, 0.0)

            fonll_up = float(np.sqrt(scale_up**2 + pdf_sigma**2 + mass_dev**2))
            fonll_dn = float(np.sqrt(scale_dn**2 + pdf_sigma**2 + mass_dev**2))
            total_up = float(np.hypot(fonll_up, model_up))
            total_dn = float(np.hypot(fonll_dn, model_dn))
            rec.update({
                f"{boundary}_campaign_central": 10.0**xc,
                f"{boundary}_fonll_band_lo": 10.0**(xref - fonll_dn),
                f"{boundary}_fonll_band_hi": 10.0**(xref + fonll_up),
                f"{boundary}_decay_model_alt": (
                    10.0**(xref + model_shift) if alt is not None else np.nan),
                f"{boundary}_band_lo": 10.0**(xref - total_dn),
                f"{boundary}_band_hi": 10.0**(xref + total_up),
                f"{boundary}_scale_up_dex": scale_up,
                f"{boundary}_scale_dn_dex": scale_dn,
                f"{boundary}_pdf_sigma_dex": pdf_sigma,
                f"{boundary}_mb_dev_dex": mass_dev,
                f"{boundary}_decay_model_up_dex": model_up,
                f"{boundary}_decay_model_dn_dex": model_dn,
                f"{boundary}_variation_missing": missing,
                f"{boundary}_open": any_open,
            })
        rows.append(rec)
    return pd.DataFrame(rows)


def collect_campaign(variations, scratch_dir, masses, n_pool, n_samples, seed,
                     central_curve, curves_out, band_out, manifest_out):
    """Require a complete campaign, publish compact curves, and combine the band."""
    scratch_dir = Path(scratch_dir)
    frames, run_records, missing = [], [], []
    for variation in variations:
        run_dir = scratch_dir / "runs" / variation["name"]
        config = _campaign_config(variation, production.MASS_GRID, n_pool, n_samples, seed)
        rows = []
        for mass in masses:
            result = _valid_result(_result_path(run_dir, mass), config)
            if result is None:
                missing.append(f"{variation['name']}:{mass:.3f}")
            else:
                rows.append(result)
        if len(rows) == len(masses):
            frame = pd.DataFrame(rows).sort_values("mass_GeV")
            frames.append(frame)
            run_records.append({
                "variation": variation["name"],
                "config_sha256": config["config_sha256"],
                "n_masses": len(frame),
                "parent_pool_seed": config["parent_pool_seed"],
            })
    if missing:
        preview = ", ".join(missing[:20])
        raise RuntimeError(
            f"campaign incomplete: {len(missing)} missing mass results; first: {preview}")

    raw = pd.concat(frames, ignore_index=True).sort_values(
        ["axis", "variation", "mass_GeV"])
    band = combine_band(raw, central_curve)
    _atomic_csv(raw, Path(curves_out))
    _atomic_csv(band, Path(band_out))

    source_manifest = Path(variations[0]["path"]).parent / "variation_manifest.json"
    manifest = {
        "artifact": "GRENDEL BC4 independent full-statistics uncertainty campaign",
        "schema_version": 1,
        "published_utc": datetime.now(timezone.utc).isoformat(),
        "producer_repo": "llpatcolliders (branch bc4-scalar)",
        "producer_git_head": _git_head(),
        "producer_code_sha256": _code_hashes(),
        "scratch_root": str(Path(scratch_dir).resolve()),
        "independence_policy": (
            "Every variation has a fresh FONLL parent sample, scalar four-vectors, "
            "geometry, decay/reconstruction MC, and sensitivity scan. No event, "
            "weight, geometry, or detector-outcome reuse/reweighting across variations."
        ),
        "campaign": {
            "n_fonll_variations": 109,
            "n_decay_model_variations": 1,
            "n_total_variations": len(variations),
            "n_masses": len(masses),
            "n_parent_pool_per_variation": n_pool,
            "n_scalar_events_per_mass": n_pool * len(production.B_SPECIES),
            "n_decay_samples_per_hit": n_samples,
            "base_seed": seed,
            "mass_grid_GeV": [float(mass) for mass in masses],
            "source_fonll_manifest_sha256": _sha256(source_manifest),
        },
        "combination": {
            "space": "log10(sin^2 theta)",
            "scale": "asymmetric envelope of central plus six coherent scale grids",
            "pdf": "sample standard deviation (ddof=1) over 100 NNPDF replicas",
            "bottom_mass": "maximum absolute displacement from mb=4.5/5.0 GeV",
            "decay_model": (
                "independent central-FONLL simulation with LO-ChPT widths below "
                "2 GeV and perturbative spectator widths above"
            ),
            "total": "independent source components added in quadrature per direction",
            "rebase": "component dex shifts applied to canonical published central curve",
        },
        "variations": [
            {
                "name": variation["name"],
                "axis": variation["axis"],
                "width_scheme": variation["width_scheme"],
                "grid_file": variation["path"].name,
                "grid_sha256": variation["sha256"],
                "trapezoid_integral_pb": variation["trapezoid_integral_pb"],
                **next(record for record in run_records
                       if record["variation"] == variation["name"]),
            }
            for variation in variations
        ],
        "outputs": {
            Path(curves_out).name: {
                "sha256": _sha256(Path(curves_out)), "rows": len(raw)},
            Path(band_out).name: {
                "sha256": _sha256(Path(band_out)), "rows": len(band)},
            Path(central_curve).name: {
                "sha256": _sha256(Path(central_curve)), "role": "canonical central"},
        },
        "software": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "platform": platform.platform(),
        },
        "limitations": [
            "No FONLL alpha_s companion grids are included.",
            "The decay alternate is a model envelope, not a Gaussian error.",
            "Detector-response and background systematics are outside this theory band.",
        ],
    }
    _atomic_json(Path(manifest_out), manifest)
    return raw, band, manifest


def campaign_status(variations, scratch_dir, masses, n_pool, n_samples, seed):
    total = len(variations) * len(masses)
    done = 0
    for variation in variations:
        config = _campaign_config(variation, production.MASS_GRID, n_pool, n_samples, seed)
        run_dir = Path(scratch_dir) / "runs" / variation["name"]
        count = sum(_valid_result(_result_path(run_dir, mass), config) is not None
                    for mass in masses)
        done += count
        print(f"{variation['name']:28s} {count:3d}/{len(masses)}")
    print(f"TOTAL {done}/{total} mass-variation results ({100 * done / total:.2f}%)")
    return done, total


def _common_options(parser):
    parser.add_argument("--grid-dir", type=Path, default=DEFAULT_GRID_DIR)
    parser.add_argument("--scratch-dir", type=Path, default=DEFAULT_SCRATCH)
    parser.add_argument("--n-pool", type=int, default=production.N_POOL_DEFAULT,
                        help="bottom parent events; three B species give 3*n-pool scalars")
    parser.add_argument("--n-samples", type=int, default=100)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--masses", type=float, nargs="+", default=None)
    parser.add_argument("--variations", nargs="+", default=None,
                        help="variation names or axes; default complete campaign")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    run_parser = sub.add_parser("run", help="run/resume selected independent variations")
    _common_options(run_parser)
    run_parser.add_argument("--workers", type=int, default=1,
                            help="independent variation processes")
    status_parser = sub.add_parser("status", help="report completed mass checkpoints")
    _common_options(status_parser)
    collect_parser = sub.add_parser("collect", help="require completeness and publish band")
    _common_options(collect_parser)
    collect_parser.add_argument("--central-curve", type=Path, default=PUBLISHED_CENTRAL)
    collect_parser.add_argument("--curves-out", type=Path, default=PUBLISHED_CURVES)
    collect_parser.add_argument("--band-out", type=Path, default=PUBLISHED_BAND)
    collect_parser.add_argument("--manifest-out", type=Path, default=PUBLISHED_MANIFEST)
    args = parser.parse_args(argv)

    source_manifest, variations = discover_variations(args.grid_dir.resolve())
    del source_manifest
    variations = _select_variations(variations, args.variations)
    masses = [float(mass) for mass in (args.masses or production.MASS_GRID)]

    if args.command == "status":
        campaign_status(
            variations, args.scratch_dir, masses,
            args.n_pool, args.n_samples, args.seed)
        return 0
    if args.command == "collect":
        if args.masses or args.variations:
            raise ValueError("published collection requires the complete mass/variation campaign")
        raw, band, _ = collect_campaign(
            variations, args.scratch_dir, masses, args.n_pool, args.n_samples,
            args.seed, args.central_curve, args.curves_out, args.band_out,
            args.manifest_out)
        print(f"published {len(raw)} variation rows -> {args.curves_out}")
        print(f"published {len(band)} band rows -> {args.band_out}")
        print(f"provenance -> {args.manifest_out}")
        return 0

    if args.workers < 1:
        raise ValueError("--workers must be positive")
    args.scratch_dir.mkdir(parents=True, exist_ok=True)
    if args.workers == 1:
        for variation in variations:
            run_variation(
                variation, args.scratch_dir, masses,
                args.n_pool, args.n_samples, args.seed)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(
                    run_variation, variation, args.scratch_dir, masses,
                    args.n_pool, args.n_samples, args.seed): variation["name"]
                for variation in variations
            }
            for future in as_completed(futures):
                name = futures[future]
                try:
                    future.result()
                    print(f"[orchestrator] complete: {name}", flush=True)
                except Exception as exc:
                    print(f"[orchestrator] FAILED: {name}: {exc}", file=sys.stderr,
                          flush=True)
                    raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
