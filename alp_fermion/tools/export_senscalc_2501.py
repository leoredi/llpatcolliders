#!/usr/bin/env python3
"""Verify and export the pinned SensCalc arXiv:2501.04525 decay inputs.

The upstream files are Wolfram MX dumps, so they must be decoded by the
official Wolfram Engine.  This wrapper verifies the exact SensCalc release and
input hashes before invoking ``export_senscalc_2501.wls``.

Example:

    python alp_fermion/tools/export_senscalc_2501.py /path/to/SensCalc
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


SENSCALC_REPOSITORY = "https://github.com/maksymovchynnikov/SensCalc"
SENSCALC_TAG = "v.1.3.3"
SENSCALC_COMMIT = "0bca050633aae16e148d47f21840fa07ff4b8724"

DECAY_DATA_SUBDIR = Path("phenomenology/ALP-fermion/decay widths")
SOURCE_FILES = {
    "widths": (
        DECAY_DATA_SUBDIR
        / "Widths-model-ALP-fermion-scale-1000.-GeV-2501.04525.m",
        "d12fb78d28edff0aa081c9fb66d829b42b4ec71202684019d7b9047ecb40b869",
    ),
    "branching_ratios": (
        DECAY_DATA_SUBDIR
        / "Br-ratios-SensCalc-model-ALP-fermion-scale-1000.-GeV-2501.04525.m",
        "34a09eed87d081bfffe79b464094741454022f79478c5b28bc236bc361049286",
    ),
    "matrix_elements": (
        DECAY_DATA_SUBDIR
        / "Matrix-elements-squared-model-ALP-fermion-scale-1000.-GeV-2501.04525.m",
        "f959c2257fa349e5af6966795db8cbf0da2e3ff8d097b5fccc3316d270ff17ed",
    ),
    "acceptance_notebook": (
        Path("codes/Acceptances/ALP-fermion.nb"),
        "8d205b65da456fb1a8fac8d1803eaa73ae839128a60d0624eeea17cacee5d4b4",
    ),
}

DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parents[1] / "data" / "senscalc_2501"
EXPORTER = Path(__file__).with_suffix(".wls")
EXPECTED_OUTPUT_FILES = (
    "widths_raw_senscalc.csv",
    "widths_bnt.csv",
    "widths_metadata.json",
    "branching_ratios.csv",
    "decay_channels.json",
    "matrix_elements.json",
    "EXPORT_MANIFEST.json",
)
MACOS_ENGINE_KERNEL = Path(
    "/Applications/Wolfram Engine.app/Contents/Resources/"
    "Wolfram Player.app/Contents/MacOS/WolframKernel"
)
MX_SYSTEM_ID_RE = re.compile(
    rb"(?:Windows|MacOSX|Linux)-[A-Za-z0-9-]+"
)


class InputError(RuntimeError):
    """The local SensCalc checkout does not match the publication pin."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def mx_system_id(path: Path) -> str | None:
    """Return the platform recorded in a Wolfram MX header, if present."""
    with path.open("rb") as handle:
        match = MX_SYSTEM_ID_RE.search(handle.read(256))
    return match.group().decode("ascii") if match else None


def source_mx_system_ids(paths: dict[str, Path]) -> dict[str, str]:
    """Map binary source keys to the system IDs embedded in their headers."""
    return {
        key: system_id
        for key, path in paths.items()
        if (system_id := mx_system_id(path)) is not None
    }


def git_head(repo: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.STDOUT,
        ).strip()
    except (OSError, subprocess.CalledProcessError) as exc:
        raise InputError(f"not a readable Git checkout: {repo}") from exc


def verify_sources(repo: Path) -> dict[str, Path]:
    repo = repo.resolve()
    head = git_head(repo)
    if head != SENSCALC_COMMIT:
        raise InputError(
            f"SensCalc HEAD is {head}, expected {SENSCALC_TAG} "
            f"({SENSCALC_COMMIT})"
        )

    verified = {}
    for key, (relative_path, expected_hash) in SOURCE_FILES.items():
        path = repo / relative_path
        if not path.is_file():
            raise InputError(f"missing SensCalc input: {path}")
        actual_hash = sha256(path)
        if actual_hash != expected_hash:
            raise InputError(
                f"SHA-256 mismatch for {relative_path}: {actual_hash}; "
                f"expected {expected_hash}"
            )
        verified[key] = path
    return verified


def wolfram_environment() -> dict[str, str]:
    env = os.environ.copy()
    if "WOLFRAMSCRIPT_KERNELPATH" not in env and MACOS_ENGINE_KERNEL.is_file():
        env["WOLFRAMSCRIPT_KERNELPATH"] = str(MACOS_ENGINE_KERNEL)
    return env


def validate_export(export_dir: Path) -> None:
    """Reject incomplete or internally inconsistent Wolfram output."""
    missing = [name for name in EXPECTED_OUTPUT_FILES
               if not (export_dir / name).is_file()]
    if missing:
        raise InputError(f"Wolfram export is missing: {', '.join(missing)}")

    try:
        manifest = json.loads((export_dir / "EXPORT_MANIFEST.json").read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise InputError("EXPORT_MANIFEST.json is not valid JSON") from exc

    if manifest.get("senscalc_tag") != SENSCALC_TAG:
        raise InputError("export manifest has the wrong SensCalc tag")
    if manifest.get("senscalc_commit") != SENSCALC_COMMIT:
        raise InputError("export manifest has the wrong SensCalc commit")
    if manifest.get("phenomenology") != "arXiv:2501.04525":
        raise InputError("export manifest has the wrong phenomenology source")
    if not isinstance(manifest.get("wolfram_system_id"), str):
        raise InputError("export manifest has no Wolfram system ID")
    expected_source_hashes = {
        key: expected_hash
        for key, (_, expected_hash) in SOURCE_FILES.items()
    }
    if manifest.get("source_sha256") != expected_source_hashes:
        raise InputError("export manifest has the wrong source hashes")
    output_hashes = manifest.get("output_sha256")
    if not isinstance(output_hashes, dict):
        raise InputError("export manifest has no output_sha256 mapping")
    for name in EXPECTED_OUTPUT_FILES:
        if name == "EXPORT_MANIFEST.json":
            continue
        expected = output_hashes.get(name)
        actual = sha256(export_dir / name)
        if expected != actual:
            raise InputError(
                f"export hash mismatch for {name}: {actual}; expected {expected}"
            )


def install_export(staging_dir: Path, output_dir: Path) -> None:
    """Atomically replace each validated data product in the destination."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in EXPECTED_OUTPUT_FILES:
        os.replace(staging_dir / name, output_dir / name)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("senscalc_root", type=Path, help="SensCalc Git checkout")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"export destination (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="verify the commit and source hashes without running Wolfram",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        verified = verify_sources(args.senscalc_root)
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print(f"SensCalc {SENSCALC_TAG} verified at {SENSCALC_COMMIT}")
    mx_system_ids = source_mx_system_ids(verified)
    for key, path in verified.items():
        suffix = f"; MX {mx_system_ids[key]}" if key in mx_system_ids else ""
        print(f"  {key}: {path.name} ({sha256(path)}{suffix})")

    if mx_system_ids:
        systems = ", ".join(sorted(set(mx_system_ids.values())))
        print(
            "note: Wolfram MX is system-dependent; these inputs record "
            f"{systems}. The export manifest records the decoder system ID."
        )

    if args.check_only:
        return 0

    wolframscript = shutil.which("wolframscript")
    if wolframscript is None:
        print(
            "error: wolframscript is required to decode the upstream MX files",
            file=sys.stderr,
        )
        return 3

    output_dir = args.output_dir.resolve()
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    try:
        with tempfile.TemporaryDirectory(
            prefix=f".{output_dir.name}-", dir=output_dir.parent
        ) as staging_name:
            staging_dir = Path(staging_name)
            command = [
                wolframscript,
                "-file",
                str(EXPORTER),
                str(args.senscalc_root.resolve()),
                str(staging_dir),
            ]
            subprocess.run(command, check=True, env=wolfram_environment())
            validate_export(staging_dir)
            install_export(staging_dir, output_dir)
    except subprocess.CalledProcessError as exc:
        print(
            "error: Wolfram export failed. If the engine reports a license "
            "problem, obtain the free developer license and run "
            "`wolframscript -activate` locally before retrying.",
            file=sys.stderr,
        )
        return exc.returncode or 4
    except InputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 5

    print(f"Exported SensCalc 2501 inputs to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
