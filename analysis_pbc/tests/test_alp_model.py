"""
Minimal ALPModel validation checks.

Run:
  cd analysis_pbc
  python tests/test_alp_model.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ANALYSIS_ROOT = Path(__file__).resolve().parents[1]
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from models.alp_model import ALPModel


def _rel_close(a: float, b: float, rel: float = 1e-3) -> bool:
    if a == b:
        return True
    if a == 0 or b == 0:
        return False
    return abs(a - b) / max(abs(a), abs(b)) <= rel


def test_lifetime_scales_as_fa_squared() -> None:
    m_a = 1.0
    ct1 = ALPModel(m_a, 1e3, "BC10").ctau0_m
    ct2 = ALPModel(m_a, 1e5, "BC10").ctau0_m
    ratio = ct2 / ct1
    assert _rel_close(ratio, 1e4, rel=1e-6), f"Expected ~1e4, got {ratio}"


def test_gammagamma_width_scales_as_m_cubed() -> None:
    fa = 1e6
    g1 = ALPModel(0.1, fa, "BC9").width_to_gammagamma()
    g2 = ALPModel(1.0, fa, "BC9").width_to_gammagamma()
    g3 = ALPModel(10.0, fa, "BC9").width_to_gammagamma()
    assert _rel_close(g2 / g1, 1e3, rel=1e-10), f"Expected ~1e3, got {g2/g1}"
    assert _rel_close(g3 / g2, 1e3, rel=1e-10), f"Expected ~1e3, got {g3/g2}"


def test_branching_ratios_sum_to_one() -> None:
    for benchmark in ["BC9", "BC10", "BC11"]:
        alp = ALPModel(5.0, 1e6, benchmark)
        total = (
            alp.branching_ratio("gamma_gamma")
            + alp.branching_ratio("ee")
            + alp.branching_ratio("mumu")
            + alp.branching_ratio("tautau")
            + alp.branching_ratio("hadrons")
        )
        assert _rel_close(total, 1.0, rel=1e-12), f"{benchmark}: sum BR={total}"


if __name__ == "__main__":
    test_lifetime_scales_as_fa_squared()
    test_gammagamma_width_scales_as_m_cubed()
    test_branching_ratios_sum_to_one()
    print("OK: ALPModel sanity checks passed")
