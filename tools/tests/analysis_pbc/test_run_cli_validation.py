#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve()
REPO_ROOT = HERE.parents[3]
RUN_SCRIPT = REPO_ROOT / "analysis_pbc" / "limits" / "run.py"


def _run_cli(args: list[str]) -> subprocess.CompletedProcess:
    cmd = [sys.executable, str(RUN_SCRIPT)] + args
    return subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True)


def _combined_output(result: subprocess.CompletedProcess) -> str:
    return f"{result.stdout}\n{result.stderr}"


def test_rejects_max_separation_not_greater_than_min():
    result = _run_cli(["--separation-mm", "1.0", "--max-separation-mm", "0.5"])
    assert result.returncode != 0
    assert "strictly greater than --separation-mm" in _combined_output(result)


def test_rejects_brvis_kappa_with_max_separation():
    result = _run_cli(
        [
            "--decay-mode",
            "brvis-kappa",
            "--max-separation-mm",
            "100.0",
        ]
    )
    assert result.returncode != 0
    assert "--max-separation-mm is not supported with --decay-mode brvis-kappa" in _combined_output(result)


def test_rejects_brvis_kappa_with_any_pair_policy():
    result = _run_cli(
        [
            "--decay-mode",
            "brvis-kappa",
            "--separation-policy",
            "any-pair-window",
        ]
    )
    assert result.returncode != 0
    assert "--separation-policy must be 'all-pairs-min'" in _combined_output(result)


def test_invalid_geometry_model_fails_parser():
    result = _run_cli(["--geometry-model", "invalid"])
    assert result.returncode != 0
    assert "invalid choice" in _combined_output(result).lower()
