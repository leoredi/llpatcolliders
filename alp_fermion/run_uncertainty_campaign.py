#!/usr/bin/env python3
"""Run exact, independently generated BC10 uncertainty variations.

Every FONLL member and C_bs point receives a fresh 600k production sample and
its own geometry/reconstruction/sensitivity tree. Gluon-surrogate variants use
the central production vectors (their production physics is identical) but own
independently simulated templates, geometry caches, and reconstruction scans.

Heavy artifacts live under ``--scratch-root``. The runner is interruption-safe:
mass CSVs and sensitivity checkpoints are atomic, and stage-completion JSON is
written only after the complete output set validates and hashes successfully.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

from alp_fermion.mass_grid import ALP_MASS_GRID
from alp_fermion.paths import TEMPLATE_DIR, add_hnl_to_path
from alp_fermion.uncertainty_campaign import (
    all_variations,
    atomic_json,
    sha256_file,
    sha256_tree,
)

add_hnl_to_path()
from config_mass_grid import format_mass_for_filename  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent


def _git_state() -> dict:
    commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
    ).strip()
    diff = subprocess.check_output(
        ["git", "diff", "--binary", "HEAD"], cwd=REPO_ROOT
    )
    import hashlib
    return {
        "commit": commit,
        "tracked_diff_sha256": hashlib.sha256(diff).hexdigest(),
        "tracked_tree_clean": not bool(diff),
    }


def _expected_vector_paths(directory: Path) -> list[Path]:
    excluded = {0.54, 0.96}
    return [
        directory / f"mA_{format_mass_for_filename(mass)}.csv"
        for mass in ALP_MASS_GRID if float(mass) not in excluded
    ]


def _expected_template_paths(directory: Path) -> list[Path]:
    excluded = {0.54, 0.96}
    return [
        directory / f"templates_{format_mass_for_filename(mass)}.npz"
        for mass in ALP_MASS_GRID if float(mass) not in excluded
    ]


def _validate_vectors(directory: Path, hash_outputs=True) -> dict:
    paths = _expected_vector_paths(directory)
    missing = [str(path) for path in paths if not path.exists() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(
            f"production tree is incomplete ({len(missing)} missing/empty); "
            f"first: {missing[0]}"
        )
    temporary = list(directory.glob("*.tmp"))
    if temporary:
        raise RuntimeError(f"production tree has partial temporary file: {temporary[0]}")
    return {
        "n_vector_files": len(paths),
        "total_vector_bytes": sum(path.stat().st_size for path in paths),
        "vector_tree_sha256": sha256_tree(paths) if hash_outputs else None,
    }


def _validate_sensitivity_csv(path: Path) -> dict:
    if not path.exists():
        raise RuntimeError(f"missing sensitivity CSV: {path}")
    frame = pd.read_csv(path)
    expected = [float(value) for value in ALP_MASS_GRID]
    observed = sorted(frame["mass_GeV"].astype(float))
    if len(frame) != len(expected) or observed != sorted(expected):
        raise RuntimeError(
            f"sensitivity checkpoint has {len(frame)} rows; expected {len(expected)}"
        )
    temporary = path.with_suffix(path.suffix + ".tmp")
    if temporary.exists():
        raise RuntimeError(f"partial sensitivity checkpoint remains: {temporary}")
    return {
        "sensitivity_csv_sha256": sha256_file(path),
        "n_sensitivity_rows": len(frame),
        "n_sensitive_rows": int(frame["has_sensitivity"].sum()),
    }


def _validate_geometry(geometry_dir: Path) -> dict:
    geometry = sorted(geometry_dir.glob("geom_*.npz"))
    if len(geometry) != 97:
        raise RuntimeError(f"geometry cache has {len(geometry)} files; expected 97")
    return {
        "n_geometry_files": len(geometry),
        "total_geometry_bytes": sum(path.stat().st_size for path in geometry),
        "geometry_tree_sha256": sha256_tree(geometry),
    }


def _validate_sensitivity(path: Path, geometry_dir: Path) -> dict:
    return {
        **_validate_sensitivity_csv(path),
        **_validate_geometry(geometry_dir),
    }


def _tree_usage(path: Path) -> dict:
    files = [item for item in path.rglob("*") if item.is_file()]
    return {
        "path": str(path),
        "n_files": len(files),
        "bytes": sum(item.stat().st_size for item in files),
    }


def _remove_generated_tree(path: Path, run_dir: Path) -> None:
    path = path.resolve()
    run_dir = run_dir.resolve()
    if path.parent != run_dir and path.parent.parent != run_dir:
        raise RuntimeError(f"refusing to compact path outside run directory: {path}")
    if path.exists():
        shutil.rmtree(path)


def _template_provenance(directory: Path) -> dict:
    paths = _expected_template_paths(directory)
    missing = [path for path in paths if not path.exists()]
    if missing:
        raise RuntimeError(
            f"template directory {directory} is missing {len(missing)} campaign "
            f"files; first: {missing[0]}"
        )
    return {
        "path": str(directory.resolve()),
        "n_files": len(paths),
        "tree_sha256": sha256_tree(paths),
        "total_bytes": sum(path.stat().st_size for path in paths),
    }


def _run_logged(command, env, log_path, label):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[{label}] {' '.join(command)}", flush=True)
    with log_path.open("a") as log:
        log.write(f"\n=== {time.strftime('%Y-%m-%dT%H:%M:%S%z')} ===\n")
        log.write("$ " + " ".join(command) + "\n")
        log.flush()
        result = subprocess.run(
            command,
            cwd=REPO_ROOT,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    if result.returncode:
        raise RuntimeError(
            f"{label} failed with exit {result.returncode}; see {log_path}"
        )


def _marker_matches(path: Path, variation, n_pool, git_state) -> bool:
    if not path.exists():
        return False
    payload = json.loads(path.read_text())
    return (
        payload.get("variation", {}).get("name") == variation["name"]
        and payload.get("variation", {}).get("grid_sha256") == variation["grid_sha256"]
        and payload.get("n_pool") == n_pool
        and payload.get("code", {}) == git_state
    )


def _validate_completed_run(run_dir: Path, marker: dict) -> None:
    output = Path(marker["sensitivity_csv"])
    current = _validate_sensitivity_csv(output)
    if current["sensitivity_csv_sha256"] != marker["sensitivity_csv_sha256"]:
        raise RuntimeError(
            f"completed {marker['variation']['name']} sensitivity checksum changed"
        )

    state = marker.get("storage_state", "full")
    if state == "full":
        geometry = _validate_geometry(run_dir / "analysis" / "geometry_cache")
        if geometry["geometry_tree_sha256"] != marker["geometry_tree_sha256"]:
            raise RuntimeError(
                f"completed {marker['variation']['name']} geometry checksum changed"
            )
        if marker["variation"]["production_mode"] == "fresh_600k":
            vectors = _validate_vectors(run_dir / "llp_4vectors")
            production = json.loads(Path(marker["production_marker"]).read_text())
            if vectors["vector_tree_sha256"] != production["vector_tree_sha256"]:
                raise RuntimeError(
                    f"completed {marker['variation']['name']} vector checksum changed"
                )
    elif state == "compacted":
        for path in (
            run_dir / "llp_4vectors",
            run_dir / "analysis" / "geometry_cache",
        ):
            if path.exists():
                raise RuntimeError(
                    f"compacted run unexpectedly retains generated tree: {path}"
                )
    elif state != "compacting":
        raise RuntimeError(f"unknown completion storage state: {state}")


def _compact_completed_run(
    completion: Path,
    marker: dict,
    run_dir: Path,
    keep_intermediates: bool,
) -> dict:
    variation = marker["variation"]
    state = marker.get("storage_state", "full")
    if variation["name"] == "central" or (
        keep_intermediates and state == "full"
    ):
        return marker
    if state == "compacted":
        return marker
    if state not in {"full", "compacting"}:
        raise RuntimeError(f"cannot compact completion in storage state {state}")

    vector_dir = run_dir / "llp_4vectors"
    geometry_dir = run_dir / "analysis" / "geometry_cache"
    targets = [geometry_dir]
    if variation["production_mode"] == "fresh_600k":
        targets.insert(0, vector_dir)

    if state == "full":
        marker = dict(marker)
        marker["storage_state"] = "compacting"
        marker["compaction"] = {
            "started_unix": time.time(),
            "targets": [_tree_usage(path) for path in targets if path.exists()],
            "retained": [
                marker["sensitivity_csv"],
                marker["production_marker"],
                str(run_dir / "production.log"),
                str(run_dir / "sensitivity.log"),
            ],
        }
        atomic_json(completion, marker)

    for target in targets:
        _remove_generated_tree(target, run_dir)

    marker = dict(marker)
    marker["storage_state"] = "compacted"
    marker["compaction"] = dict(marker["compaction"])
    marker["compaction"]["completed_unix"] = time.time()
    marker["compaction"]["reclaimed_bytes"] = sum(
        item["bytes"] for item in marker["compaction"]["targets"]
    )
    atomic_json(completion, marker)
    print(
        f"[{variation['name']}] compacted "
        f"{marker['compaction']['reclaimed_bytes'] / 2**30:.1f} GiB",
        flush=True,
    )
    return marker


def _environment(run_dir, vector_dir, analysis_dir, template_dir, variation):
    env = os.environ.copy()
    env.update({
        "ALP_TMP_DIR": str(run_dir),
        "ALP_LLP_VECTORS_DIR": str(vector_dir),
        "ALP_ANALYSIS_DIR": str(analysis_dir),
        "ALP_GEOM_CACHE_DIR": str(analysis_dir / "geometry_cache"),
        "ALP_TEMPLATE_DIR": str(template_dir),
        "HNL_FONLL_BOTTOM_GRID": variation["grid_path"],
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
    })
    return env


def run_variation(variation, args, git_state, template_cache):
    scratch = Path(args.scratch_root).expanduser().resolve()
    run_dir = scratch / "runs" / variation["name"]
    run_dir.mkdir(parents=True, exist_ok=True)
    analysis_dir = run_dir / "analysis"
    central_dir = scratch / "runs" / "central"
    if variation["template_variant"] == "central":
        template_dir = Path(args.central_template_dir).expanduser().resolve()
    else:
        template_dir = (
            Path(args.gluon_template_root).expanduser().resolve()
            / variation["template_variant"]
        )
    template_key = str(template_dir)
    if template_key not in template_cache:
        template_cache[template_key] = _template_provenance(template_dir)

    completion = run_dir / "variation.complete.json"
    if _marker_matches(completion, variation, args.n_pool, git_state):
        marker = json.loads(completion.read_text())
        _validate_completed_run(run_dir, marker)
        _compact_completed_run(
            completion, marker, run_dir, args.keep_intermediates
        )
        print(f"[{variation['name']}] variation already complete", flush=True)
        return
    if completion.exists():
        raise RuntimeError(
            f"stale completion marker for {variation['name']}; use a new scratch root"
        )

    if variation["production_mode"] == "fresh_600k":
        vector_dir = run_dir / "llp_4vectors"
        production_marker = run_dir / "production.complete.json"
        if _marker_matches(
            production_marker, variation, args.n_pool, git_state
        ):
            info = _validate_vectors(vector_dir, hash_outputs=True)
            marker = json.loads(production_marker.read_text())
            if info["vector_tree_sha256"] != marker["vector_tree_sha256"]:
                raise RuntimeError(
                    f"completed production checksum changed for {variation['name']}"
                )
            print(f"[{variation['name']}] production already complete", flush=True)
        else:
            if production_marker.exists():
                raise RuntimeError(
                    f"stale production marker for {variation['name']}; "
                    "use a new scratch root or remove that generated run explicitly"
                )
            env = _environment(
                run_dir, vector_dir, analysis_dir, template_dir, variation
            )
            command = [
                sys.executable, "-u", "-m", "alp_fermion.alp_production",
                "--n-pool", str(args.n_pool),
                "--seed", str(variation["production_seed"]),
                "--cbs-amplitude-scale", str(variation["cbs_amplitude_scale"]),
                "--resume",
            ]
            started = time.time()
            production_log = run_dir / "production.log"
            _run_logged(
                command, env, production_log,
                f"{variation['name']}:production",
            )
            vector_info = _validate_vectors(vector_dir)
            atomic_json(production_marker, {
                "variation": variation,
                "n_pool": args.n_pool,
                "code": git_state,
                "command": command,
                "log": str(production_log),
                "log_sha256": sha256_file(production_log),
                "started_unix": started,
                "completed_unix": time.time(),
                **vector_info,
            })
    else:
        central_marker = central_dir / "production.complete.json"
        if not central_marker.exists():
            raise RuntimeError(
                f"{variation['name']} requires completed central production: "
                f"{central_marker}"
            )
        vector_dir = central_dir / "llp_4vectors"
        _validate_vectors(vector_dir, hash_outputs=False)
        atomic_json(run_dir / "production.reference.json", {
            "variation": variation,
            "central_production_marker": str(central_marker),
            "central_production_marker_sha256": sha256_file(central_marker),
            "vector_dir": str(vector_dir),
        })
    env = _environment(run_dir, vector_dir, analysis_dir, template_dir, variation)
    output = analysis_dir / "bc10_sensitivity.csv"
    command = [
        sys.executable, "-u", "-m", "alp_fermion.sensitivity",
        "--output", str(output), "--resume",
    ]
    started = time.time()
    sensitivity_log = run_dir / "sensitivity.log"
    _run_logged(
        command, env, sensitivity_log,
        f"{variation['name']}:sensitivity",
    )
    analysis_info = _validate_sensitivity(output, analysis_dir / "geometry_cache")
    production_marker = (
        run_dir / "production.complete.json"
        if variation["production_mode"] == "fresh_600k"
        else run_dir / "production.reference.json"
    )
    atomic_json(completion, {
        "variation": variation,
        "n_pool": args.n_pool,
        "code": git_state,
        "python": {
            "executable": sys.executable,
            "version": platform.python_version(),
        },
        "template": template_cache[template_key],
        "production_marker": str(production_marker),
        "production_marker_sha256": sha256_file(production_marker),
        "sensitivity_csv": str(output),
        "command": command,
        "sensitivity_log": str(sensitivity_log),
        "sensitivity_log_sha256": sha256_file(sensitivity_log),
        "storage_state": "full",
        "regeneration": {
            "grid_path": variation["grid_path"],
            "grid_sha256": variation["grid_sha256"],
            "production_seed": variation["production_seed"],
            "n_pool": args.n_pool,
            "cbs_amplitude_scale": variation["cbs_amplitude_scale"],
            "production_mode": variation["production_mode"],
            "template_path": template_cache[template_key]["path"],
            "template_tree_sha256": template_cache[template_key]["tree_sha256"],
            "production_command": (
                json.loads(production_marker.read_text()).get("command")
                if variation["production_mode"] == "fresh_600k"
                else "exact reuse of central production vectors"
            ),
            "sensitivity_command": command,
        },
        "started_unix": started,
        "completed_unix": time.time(),
        **analysis_info,
    })
    marker = json.loads(completion.read_text())
    _compact_completed_run(
        completion, marker, run_dir, args.keep_intermediates
    )
    print(f"[{variation['name']}] COMPLETE", flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scratch-root",
        default=os.environ.get("ALP_UNCERTAINTY_SCRATCH"),
        required=os.environ.get("ALP_UNCERTAINTY_SCRATCH") is None,
    )
    parser.add_argument("--grid-dir", type=Path, required=True)
    parser.add_argument(
        "--central-template-dir", type=Path,
        default=os.environ.get("ALP_CENTRAL_TEMPLATE_DIR", TEMPLATE_DIR),
    )
    parser.add_argument(
        "--gluon-template-root", type=Path,
        default=os.environ.get("ALP_GLUON_TEMPLATE_ROOT"),
        required=os.environ.get("ALP_GLUON_TEMPLATE_ROOT") is None,
    )
    parser.add_argument(
        "--axes", nargs="+", choices=("fonll", "decay_gg", "cbs"),
        default=["fonll"],
    )
    parser.add_argument("--only", nargs="+", default=None)
    parser.add_argument("--n-pool", type=int, default=600_000)
    parser.add_argument("--worker-index", type=int, default=0)
    parser.add_argument("--worker-count", type=int, default=1)
    parser.add_argument("--max-variations", type=int, default=None)
    parser.add_argument(
        "--keep-intermediates", action="store_true",
        help=(
            "retain non-central production vectors and geometry caches after "
            "their validated completion marker is written"
        ),
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="accepted for explicit restart intent; stages always resume atomically",
    )
    args = parser.parse_args(argv)
    if not 0 <= args.worker_index < args.worker_count:
        parser.error("worker-index must satisfy 0 <= index < worker-count")
    scratch = Path(args.scratch_root).expanduser().resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    variations = [
        variation for variation in all_variations(args.grid_dir)
        if variation["campaign_axis"] in args.axes
    ]
    if args.only:
        selected = set(args.only)
        variations = [v for v in variations if v["name"] in selected]
        missing = selected - {v["name"] for v in variations}
        if missing:
            parser.error(f"unknown/ineligible variations: {sorted(missing)}")
    variations = [
        variation for index, variation in enumerate(variations)
        if index % args.worker_count == args.worker_index
    ]
    if args.max_variations is not None:
        variations = variations[:args.max_variations]
    git_state = _git_state()
    print(
        f"BC10 exact uncertainty worker {args.worker_index}/{args.worker_count}: "
        f"{len(variations)} variations under {scratch}",
        flush=True,
    )
    template_cache = {}
    for variation in variations:
        run_variation(variation, args, git_state, template_cache)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
