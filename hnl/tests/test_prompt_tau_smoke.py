"""End-to-end smoke test for the prompt-tau driver.

Skipped automatically if MG5 or LHAPDF aren't available in this environment
(matches the existing test_smoke.py pattern for the meson side). Runs a tiny
1k-event Stage 1 + single mass Stage 2, asserts the output schemas.

Note: this test mutates the on-disk Stage 1 pool at vendored/tau_pool.csv
(re-generates if not present). The pool is small (~1k events ~ 150 KB) and
deterministic from MG5's internal seeding, so subsequent normal runs are
unaffected.
"""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from production.madgraph._mg5_common import MG5_EXE, LHAPDF_CONFIG
from production.madgraph.run_tau_production import POOL_CSV

needs_mg5 = pytest.mark.skipif(
    not MG5_EXE.exists() or not LHAPDF_CONFIG.exists(),
    reason="MG5 or LHAPDF not available in this environment",
)

HNL_ROOT = Path(__file__).resolve().parent.parent


@needs_mg5
def test_prompt_tau_smoke_e2e(tmp_path):
    """Drive run_tau_production --test; assert pool + per-flavor CSV."""
    out_csv = HNL_ROOT / "output" / "llp_4vectors" / "Umu" / "tau" / "mN_1p000.csv"
    # Snapshot what's there so we can avoid polluting the actual production
    # output dir with a 1000-event smoke if something else cares.
    pre_existed = out_csv.exists()
    pre_mtime = out_csv.stat().st_mtime if pre_existed else None

    env = os.environ.copy()
    result = subprocess.run(
        [sys.executable, "-m", "production.madgraph.run_tau_production",
         "--test", "--nb-core", "2"],
        cwd=HNL_ROOT, env=env, capture_output=True, text=True, timeout=900,
    )
    assert result.returncode == 0, (
        f"run_tau_production --test failed:\nSTDOUT: {result.stdout}\n"
        f"STDERR: {result.stderr}"
    )

    # Stage 1 product: 6-column CSV with origin column.
    assert POOL_CSV.exists() and POOL_CSV.stat().st_size > 0
    pool = np.loadtxt(POOL_CSV, delimiter=",")
    if pool.ndim == 1:
        pool = pool.reshape(1, -1)
    assert pool.shape[1] == 6, f"pool schema is {pool.shape[1]} cols; expected 6"
    origins = np.unique(pool[:, 5].astype(int))
    # Every origin should be a real PDG (W±=24, Z=23, or initial-state parton
    # / gluon). Allow gluon (21), quarks (±1..±6), W±, Z, photon (22) — but
    # reject zeros (placeholder for malformed lines) and unknown garbage.
    allowed = {21, 22, 23, 24, -24} | {q for q in range(-6, 7) if q != 0}
    assert all(int(o) in allowed for o in origins), f"unexpected origin PDGs: {origins}"
    # Must contain at least one W± (from W -> tau nu) — guards against MG5
    # silently dropping that branch.
    assert any(abs(int(o)) == 24 for o in origins)

    # Stage 2 product for Umu: 5-column CSV, on-shell N.
    assert out_csv.exists() and out_csv.stat().st_size > 0
    csv = np.loadtxt(out_csv, delimiter=",")
    if csv.ndim == 1:
        csv = csv.reshape(1, -1)
    assert csv.shape[1] == 5
    w, E, px, py, pz = csv.T
    m2 = E ** 2 - px ** 2 - py ** 2 - pz ** 2
    # m_N = 1 GeV → m^2 = 1; the spread is dominated by the 3-body component.
    assert np.allclose(m2.mean(), 1.0, atol=5e-3)
    # All weights positive (no negative-weight LHE for tree-level W/Z + DY).
    assert (w > 0).all()
