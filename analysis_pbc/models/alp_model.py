"""
models/alp_model.py

ALP Model for transverse LLP detector sensitivity calculations.

This module provides the ALPModel class that computes:
- Decay widths for all relevant channels
- Proper decay length (cτ0)
- Effective couplings for different PBC benchmarks
"""

from __future__ import annotations

import numpy as np

# ============================================================
# PHYSICAL CONSTANTS - Use these exact values
# ============================================================

# Fundamental constants
ALPHA_EM = 7.2973525693e-3  # Fine structure constant (1/137.036)
HBAR_C_GEV_M = 1.97327e-16  # ℏc in GeV·m (for lifetime conversion)
HBAR_GEV_S = 6.582119e-25  # ℏ in GeV·s

# Particle masses in GeV
M_ELECTRON = 0.000510999
M_MUON = 0.105658
M_TAU = 1.77686
M_PION_CHARGED = 0.13957

# Decay constants in GeV
F_PI = 0.0922


class ALPModel:
    """
    Axion-Like Particle model for sensitivity calculations.

    Parameters
    ----------
    mass_GeV : float
        ALP mass in GeV
    fa_GeV : float
        Decay constant f_a in GeV (coupling strength = 1/f_a)
    benchmark : str
        PBC benchmark: 'BC9' (photon), 'BC10' (fermion), 'BC11' (gluon)
    C_gg : float, optional
        Photon Wilson coefficient (default from benchmark)
    C_GG : float, optional
        Gluon Wilson coefficient (default from benchmark)
    c_ff : float, optional
        Universal fermion coefficient (default from benchmark)
    """

    def __init__(
        self,
        mass_GeV: float,
        fa_GeV: float,
        benchmark: str = "BC10",
        *,
        C_gg: float | None = None,
        C_GG: float | None = None,
        c_ff: float | None = None,
    ) -> None:
        self.mass = float(mass_GeV)
        self.fa = float(fa_GeV)
        self.benchmark = str(benchmark)

        benchmarks = {
            "BC9": {"C_gg": 1.0, "C_GG": 0.0, "c_ff": 0.0},
            "BC10": {"C_gg": 0.0, "C_GG": 0.0, "c_ff": 1.0},
            "BC11": {"C_gg": 0.0, "C_GG": 1.0, "c_ff": 0.0},
        }
        bm = benchmarks.get(self.benchmark, benchmarks["BC10"])

        self.C_gg = float(C_gg) if C_gg is not None else float(bm["C_gg"])
        self.C_GG = float(C_GG) if C_GG is not None else float(bm["C_GG"])
        self.c_ff = float(c_ff) if c_ff is not None else float(bm["c_ff"])

    # ========== EFFECTIVE COUPLINGS ==========

    @property
    def g_agg(self) -> float:
        """Effective photon coupling g_aγγ in GeV⁻¹."""
        return self.C_gg * ALPHA_EM / (np.pi * self.fa)

    @property
    def g_aGG(self) -> float:
        """Effective gluon coupling g_aGG in GeV⁻¹."""
        alpha_s = 0.118  # Approximate αs at O(1–10) GeV
        return self.C_GG * alpha_s / (np.pi * self.fa)

    # ========== DECAY WIDTHS ==========

    def width_to_gammagamma(self) -> float:
        """
        Γ(a → γγ) in GeV.

        Γ = g_aγγ² × m_a³ / (64π)
        """
        if self.C_gg == 0.0:
            return 0.0
        return (self.g_agg**2 * self.mass**3) / (64.0 * np.pi)

    def width_to_ll(self, m_lepton: float) -> float:
        """
        Γ(a → ℓ⁺ℓ⁻) in GeV.

        Γ = |c_ℓ|² × m_a × m_ℓ² / (8π f_a²) × √(1 - 4m_ℓ²/m_a²)
        """
        if self.c_ff == 0.0:
            return 0.0
        if self.mass < 2.0 * m_lepton:
            return 0.0

        beta = np.sqrt(1.0 - 4.0 * m_lepton**2 / self.mass**2)
        return (self.c_ff**2 * self.mass * m_lepton**2 * beta) / (8.0 * np.pi * self.fa**2)

    def width_to_ee(self) -> float:
        return self.width_to_ll(M_ELECTRON)

    def width_to_mumu(self) -> float:
        return self.width_to_ll(M_MUON)

    def width_to_tautau(self) -> float:
        return self.width_to_ll(M_TAU)

    def width_to_gg_perturbative(self) -> float:
        """
        Γ(a → gg) in GeV (perturbative QCD; valid for m_a ≳ 2–3 GeV).

        Γ = α_s² × |C_GG|² × m_a³ / (8π³ f_a²)
        """
        if self.C_GG == 0.0:
            return 0.0
        if self.mass < 2.0:
            return 0.0
        alpha_s = 0.118
        return (alpha_s**2 * self.C_GG**2 * self.mass**3) / (8.0 * np.pi**3 * self.fa**2)

    def width_to_3pi_chpt(self) -> float:
        """
        Γ(a → 3π) in GeV (ChPT-inspired; for m_a < ~2 GeV).

        This is an approximate scaling: Γ ∝ C_GG² m_a^5 / (f_a² f_π⁴) × phase_space.
        """
        if self.C_GG == 0.0:
            return 0.0

        m_pi = M_PION_CHARGED
        if self.mass < 3.0 * m_pi:
            return 0.0

        phase_space = (1.0 - 9.0 * m_pi**2 / self.mass**2) ** 2.5
        prefactor = 1.0 / (256.0 * np.pi**3)
        return prefactor * (self.C_GG**2 * self.mass**5 * phase_space) / (self.fa**2 * F_PI**4)

    def width_to_hadrons(self) -> float:
        """Inclusive hadronic width proxy (gg at high mass; 3π below)."""
        if self.C_GG == 0.0:
            return 0.0
        if self.mass > 3.0:
            return self.width_to_gg_perturbative()
        if self.mass > 3.0 * M_PION_CHARGED:
            return self.width_to_3pi_chpt()
        return 0.0

    # ========== TOTAL WIDTH AND LIFETIME ==========

    @property
    def total_width(self) -> float:
        width = 0.0
        width += self.width_to_gammagamma()
        width += self.width_to_ee()
        width += self.width_to_mumu()
        width += self.width_to_tautau()
        width += self.width_to_hadrons()
        return float(width)

    @property
    def lifetime_seconds(self) -> float:
        if self.total_width == 0.0:
            return float("inf")
        return float(HBAR_GEV_S / self.total_width)

    @property
    def ctau0_m(self) -> float:
        if self.total_width == 0.0:
            return float("inf")
        return float(HBAR_C_GEV_M / self.total_width)

    # ========== BRANCHING RATIOS ==========

    def branching_ratio(self, channel: str) -> float:
        if self.total_width == 0.0:
            return 0.0
        widths = {
            "gamma_gamma": self.width_to_gammagamma(),
            "ee": self.width_to_ee(),
            "mumu": self.width_to_mumu(),
            "tautau": self.width_to_tautau(),
            "hadrons": self.width_to_hadrons(),
        }
        return float(widths.get(channel, 0.0) / self.total_width)

    # ========== VISIBLE BRANCHING RATIO ==========

    @property
    def visible_br(self) -> float:
        """
        Approximate fraction of decays producing ≥2 charged tracks.

        For a tracking-only detector, γγ is mostly invisible (except conversions).
        """
        br_visible = 0.0
        br_visible += self.branching_ratio("ee")
        br_visible += self.branching_ratio("mumu")
        br_visible += self.branching_ratio("tautau")
        br_visible += 0.9 * self.branching_ratio("hadrons")
        br_visible += 0.01 * self.branching_ratio("gamma_gamma")
        return float(br_visible)

    def __repr__(self) -> str:
        return (
            f"ALPModel(m_a={self.mass:.4g} GeV, f_a={self.fa:.4g} GeV, "
            f"benchmark='{self.benchmark}', cτ0={self.ctau0_m:.4g} m)"
        )

