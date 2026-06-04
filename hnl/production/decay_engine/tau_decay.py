"""
production/decay_engine/tau_decay.py

Shared τ → N + X sampler used by both τ-source drivers:

  - generate_induced_tau.py        (τ from Ds, B+ meson decays)
  - madgraph/run_tau_production.py (τ from prompt W → τν, Z → ττ)

Holds the HNLCalc-driven branching-fraction computation and the
polarisation-aware 2-body / 3-body decay sampling logic. Callers supply the
τ four-vectors and a polarisation `asymmetry` value; this module is
independent of where the τ came from.
"""

import numpy as np
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "vendored" / "HNLCalc"))

from production.constants import M_TAU
from production.decay_engine.kinematics import (
    decay_2body_polarized, decay_3body_flat,
)


# --- Tau polarization model -------------------------------------------------
# A pseudoscalar P± → τ± ν (Ds, B+, K+, ...) or W± → τ± ν fixes the τ helicity
# exactly: with a left-handed ν and the (V-A) charged current, τ+ comes out
# helicity -1 and τ- helicity +1 in the parent rest frame. By CP the N energy
# spectrum is the same for both charges, so a single asymmetry value covers
# both (see the comment in compute_tau_production_br_components for the
# explicit (1 + α·P·cosθ) form).
#
# The spin-analyzing power α is mass- and channel-dependent and would need
# the full decay matrix element to pin down exactly; we use the chiral
# (maximal) limit |α| = 1, exact as m_N → 0 and a documented upper bound on
# the polarisation effect at finite m_N. Set TAU_2BODY_ANALYZING_POWER = 0
# to recover the previous isotropic treatment.
#
# For Drell-Yan-like sources (Z → τ+τ-) the τ is *not* fully polarised; pass
# asymmetry=0 to sample_hnl_from_tau for those.
TAU_POLARIZATION = -1.0
TAU_2BODY_ANALYZING_POWER = 1.0
TAU_2BODY_ASYMMETRY = TAU_POLARIZATION * TAU_2BODY_ANALYZING_POWER


# Tau 2-body hadronic channels into N (parent_pdg = 15 (tau-))
# daughter is meson (anti-meson for tau- decay)
TAU_2BODY_MESON_PDGS = [211, 321, 213, 323]   # pi+, K+, rho+, K*+ (HNLCalc convention)


def init_hnlcalc(flavor):
    """Initialize HNLCalc with unit coupling for the given flavor."""
    from HNLCalc import HNLCalc
    if flavor == "Ue":
        return HNLCalc(ve=1, vmu=0, vtau=0)
    elif flavor == "Umu":
        return HNLCalc(ve=0, vmu=1, vtau=0)
    elif flavor == "Utau":
        return HNLCalc(ve=0, vmu=0, vtau=1)
    else:
        raise ValueError(f"Unknown flavor: {flavor}")


def _eval_tau_2body_br(hnl, meson_pdg, m_N):
    """BR(tau -> meson N) at U^2 = 1 via HNLCalc."""
    # All TAU_2BODY_MESON_PDGS daughters are charged and supported by HNLCalc;
    # no catch-all guard, so any unexpected error surfaces instead of returning 0.
    br_expr = hnl.get_2body_br_tau(15, meson_pdg)
    mass = m_N        # noqa: F841 - referenced by eval'd HNLCalc expression
    coupling = 1.0    # noqa: F841
    br_val = eval(br_expr)
    if np.isnan(br_val) or br_val < 0:
        return 0.0
    return float(br_val)


def _eval_tau_3body_br(hnl, lep_pid, nu_pid, m_N):
    """BR(tau- -> lep- nu N) at U^2 = 1 via HNLCalc, integrated over phase space."""
    # Leptonic tau channels are fully parameterized in HNLCalc; no catch-all guard,
    # so any unexpected error surfaces instead of being silently turned into 0.
    dbr = hnl.get_3body_dbr_tau(15, -lep_pid, nu_pid)
    m_lep = hnl.masses(lep_pid)
    br_val = hnl.integrate_3body_br(
        dbr, m_N, M_TAU, m_lep, 0.0,
        coupling=1.0, nsample=500, integration="dE",
    )
    if br_val is None or np.isnan(br_val) or br_val < 0:
        return 0.0
    return float(br_val)


def compute_tau_production_br_components(hnl, m_N):
    """
    Total BR(tau -> N + X) at U^2 = 1, plus per-channel breakdown for sampling.

    Returns
    -------
    br_2body_channels : list of (meson_mass, br_value)
        2-body hadronic channels with non-zero BR.
    br_3body_channels : list of (lep_mass, br_value)
        3-body leptonic channels (tau- -> lep- nu N) with non-zero BR.
    br_total : float
        Sum of all BRs.
    """
    br_2body_channels = []
    for meson_pdg in TAU_2BODY_MESON_PDGS:
        if m_N >= M_TAU - hnl.masses(meson_pdg):
            continue
        br = _eval_tau_2body_br(hnl, meson_pdg, m_N)
        if br > 0:
            br_2body_channels.append((hnl.masses(meson_pdg), br))

    br_3body_channels = []
    for lep_pid in [11, 13]:
        m_lep = hnl.masses(lep_pid)
        if m_N >= M_TAU - m_lep:
            continue
        # tau- -> lep- nu_tau N   (neutrino is nu_tau, pdg=16)
        br_nt = _eval_tau_3body_br(hnl, lep_pid, 16, m_N)
        if br_nt > 0:
            br_3body_channels.append((m_lep, br_nt))
        # tau- -> lep- nu_lep_bar N   (neutrino is anti-nu_lep, pdg = lep_pid+1)
        br_nl = _eval_tau_3body_br(hnl, lep_pid, lep_pid + 1, m_N)
        if br_nl > 0:
            br_3body_channels.append((m_lep, br_nl))

    br_total = sum(br for _, br in br_2body_channels) + sum(br for _, br in br_3body_channels)
    return br_2body_channels, br_3body_channels, br_total


def sample_hnl_from_tau(tau_E, tau_px, tau_py, tau_pz, m_N,
                        br_2body_channels, br_3body_channels, br_total, rng,
                        asymmetry=TAU_2BODY_ASYMMETRY):
    """Decay each tau to (N + X) with channel selection weighted by BR.

    ``asymmetry`` (= analyzing_power × P_tau) sets the longitudinal-polarization
    angular weight ``1 + asymmetry·cosθ`` for the N in the tau rest frame of the
    2-body hadronic modes; 0.0 reproduces isotropic decay. Scalar only — for
    a mix of W (polarised) and Z (unpolarised) taus, call once per category
    with the appropriate scalar.
    """
    n_events = len(tau_E)
    hnl_4v = np.empty((n_events, 4))

    br_2body_total = sum(br for _, br in br_2body_channels)
    p_2body = br_2body_total / br_total if br_total > 0 else 0.0
    use_2body = rng.random(n_events) < p_2body
    use_3body = ~use_2body

    # 2-body: tau -> meson + N. Pick channel weighted by per-channel BR.
    if use_2body.any() and br_2body_channels:
        br_arr = np.array([br for _, br in br_2body_channels], dtype=float)
        prob = br_arr / br_arr.sum()
        ch_idx = rng.choice(len(br_2body_channels), size=use_2body.sum(), p=prob)
        evt_idx = np.where(use_2body)[0]
        for k in np.unique(ch_idx):
            m_meson, _ = br_2body_channels[int(k)]
            sel = evt_idx[ch_idx == k]
            # N (second daughter) is the spin-analyzed particle; its polar angle
            # follows the tau longitudinal polarization (scalar asymmetry).
            _, hnl_2b = decay_2body_polarized(
                tau_E[sel], tau_px[sel], tau_py[sel], tau_pz[sel],
                M_TAU, m_meson, m_N, asymmetry=asymmetry, rng=rng,
            )
            hnl_4v[sel] = hnl_2b

    # 3-body: tau -> lep + nu + N. Flat phase space.
    if use_3body.any() and br_3body_channels:
        br_arr = np.array([br for _, br in br_3body_channels], dtype=float)
        prob = br_arr / br_arr.sum()
        ch_idx = rng.choice(len(br_3body_channels), size=use_3body.sum(), p=prob)
        evt_idx = np.where(use_3body)[0]
        for k in np.unique(ch_idx):
            m_lep, _ = br_3body_channels[int(k)]
            sel = evt_idx[ch_idx == k]
            _, _, hnl_3b = decay_3body_flat(
                tau_E[sel], tau_px[sel], tau_py[sel], tau_pz[sel],
                M_TAU, m_lep, 0.0, m_N, rng=rng,
            )
            hnl_4v[sel] = hnl_3b

    return hnl_4v, use_2body, use_3body
