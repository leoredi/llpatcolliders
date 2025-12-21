"""
production/alp_production.py

Analytic production branching ratios for ALP sensitivity projections.

These are lightweight, order-of-magnitude expressions suitable for fast scans.
For precision studies, replace these with a dedicated ALP EFT implementation.
"""

from __future__ import annotations

import numpy as np

from config.production_xsecs import (
    SIGMA_H_TOTAL_PB,
    SIGMA_Z_TOTAL_PB,
)

# ============================================================
# PHYSICAL CONSTANTS - Use these exact values
# ============================================================

ALPHA_EM = 7.2973525693e-3
HBAR_GEV_S = 6.582119e-25

M_PION_CHARGED = 0.13957
M_PION_NEUTRAL = 0.13498
M_KAON_CHARGED = 0.49368
M_KAON_NEUTRAL = 0.49761
M_B_CHARGED = 5.27934
M_B_NEUTRAL = 5.27965
M_BS = 5.36688
M_HIGGS = 125.25
M_Z = 91.1876
M_W = 80.377

V_TB = 0.999105
V_TS = 0.04110
V_US = 0.2243
V_UD = 0.97435

SIN2_THETA_W = 0.23122
COS2_THETA_W = 0.76878
G_FERMI = 1.1663788e-5
VEV = 246.22

TAU_B_CHARGED = 1.638e-12
TAU_B_NEUTRAL = 1.519e-12
TAU_BS = 1.515e-12
TAU_KAON_CHARGED = 1.2380e-8
TAU_KAON_LONG = 5.116e-8


def _kallen(a: float, b: float, c: float) -> float:
    return a**2 + b**2 + c**2 - 2 * a * b - 2 * b * c - 2 * c * a


def br_higgs_to_aa(m_a_GeV: float, coupling_lambda_aha: float) -> float:
    if m_a_GeV > M_HIGGS / 2.0:
        return 0.0

    beta = np.sqrt(1.0 - 4.0 * m_a_GeV**2 / M_HIGGS**2)
    gamma_h_aa = (coupling_lambda_aha**2 * VEV**2 * beta) / (32.0 * np.pi * M_HIGGS)
    gamma_h_sm = 4.07e-3
    return float(gamma_h_aa / (gamma_h_sm + gamma_h_aa))


def br_higgs_to_Za(m_a_GeV: float, C_aZh: float, *, Lambda_GeV: float = 1000.0) -> float:
    if m_a_GeV > M_HIGGS - M_Z:
        return 0.0

    x_Z = M_Z**2 / M_HIGGS**2
    x_a = m_a_GeV**2 / M_HIGGS**2
    lam = _kallen(1.0, x_Z, x_a)
    if lam < 0.0:
        return 0.0

    gamma_h_Za = (C_aZh**2 * M_HIGGS**3 * lam**1.5) / (16.0 * np.pi * Lambda_GeV**2)
    gamma_h_sm = 4.07e-3
    return float(gamma_h_Za / (gamma_h_sm + gamma_h_Za))


def br_Z_to_gamma_a(m_a_GeV: float, C_gamma_Z: float, fa_GeV: float) -> float:
    if m_a_GeV > M_Z:
        return 0.0

    x = m_a_GeV**2 / M_Z**2
    phase_space = (1.0 - x) ** 3
    gamma_Z_ga = (
        (ALPHA_EM**3 * M_Z**3 * phase_space * C_gamma_Z**2)
        / (48.0 * np.pi**2 * SIN2_THETA_W * COS2_THETA_W * fa_GeV**2)
    )
    gamma_Z_sm = 2.4952
    return float(gamma_Z_ga / gamma_Z_sm)


def br_B_to_Ka(m_a_GeV: float, fa_GeV: float, C_WW: float, *, meson: str = "B+") -> float:
    m_B = {"B+": M_B_CHARGED, "B0": M_B_NEUTRAL, "Bs": M_BS}[meson]
    m_K = M_KAON_CHARGED if meson == "B+" else M_KAON_NEUTRAL

    if m_a_GeV > m_B - m_K:
        return 0.0
    if C_WW == 0.0:
        return 0.0

    f_plus_0 = 0.33
    M_TOP = 172.76
    C_bs = (
        (G_FERMI * ALPHA_EM / (4.0 * np.sqrt(2) * np.pi))
        * V_TB
        * V_TS
        * C_WW
        * M_TOP**2
        / fa_GeV
    )

    lam = _kallen(m_B**2, m_K**2, m_a_GeV**2)
    if lam < 0.0:
        return 0.0

    tau_B = {"B+": TAU_B_CHARGED, "B0": TAU_B_NEUTRAL, "Bs": TAU_BS}[meson]
    gamma_B_Ka = (C_bs**2 * f_plus_0**2 * lam**1.5) / (16.0 * np.pi * m_B**3)
    return float(gamma_B_Ka * tau_B / HBAR_GEV_S)


def br_K_to_pi_a(m_a_GeV: float, fa_GeV: float, C_WW: float, *, meson: str = "K+") -> float:
    m_K = M_KAON_CHARGED if meson == "K+" else M_KAON_NEUTRAL
    m_pi = M_PION_CHARGED if meson == "K+" else M_PION_NEUTRAL

    if m_a_GeV > m_K - m_pi:
        return 0.0
    if C_WW == 0.0:
        return 0.0

    # TODO: Implement full formula from arXiv:2110.10698
    return 0.0


def n_parents_hllhc(parent_type: str, *, lumi_fb: float = 3000.0) -> float:
    """
    Number of parent particles produced at HL-LHC (integrated luminosity in fb^-1).
    """
    lumi_pb = float(lumi_fb) * 1000.0
    if parent_type == "higgs":
        return lumi_pb * SIGMA_H_TOTAL_PB
    if parent_type == "Z":
        return lumi_pb * SIGMA_Z_TOTAL_PB
    raise ValueError(f"Unknown parent_type: {parent_type}")

