#!/usr/bin/env python3
"""
production/decay_engine/generate_kaon_csvs.py

Driver for kaon -> HNL production:  K+ -> l+ N  (2-body, dominant) and
K+ -> pi0 l+ N (3-body), the leading HNL source below ~0.5 GeV.

FONLL only covers heavy quarks, so unlike the B/D/Bc path there is no tabulated
kaon spectrum. We sample charged kaons from a parametrized soft-QCD flux
(Tsallis pT, Gaussian rapidity; see production/constants.py) and then reuse the
exact 2-body/3-body decay + BR machinery from generate_meson_csvs.

Weight per HNL 4-vector i, in pb at U^2 = 1:

    w_i = SIGMA_KAON_PB * BR(K -> N + X | U^2=1) / N_pool

SIGMA_KAON_PB already bundles both charges (K+ + K-), matching the
particle+antiparticle convention of the meson channels, so no extra factor 2.
The absolute kaon flux is approximate (see constants.py) and should be treated
as a systematic.

Output: output/llp_4vectors/{Ue,Umu,Utau}/Kmeson/mN_{mass}.csv
Format: headerless, 5 columns: weight,E,px,py,pz
"""

import random
import sys
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "vendored" / "HNLCalc"))

from config_mass_grid import MASS_GRID, format_mass_for_filename
from production.constants import (
    M_KAON, LEPTON_MASSES, FLAVOR_TO_LEPTON_PDG, SIGMA_KAON_PB,
    KAON_TSALLIS_T, KAON_TSALLIS_N, KAON_PT_MAX, KAON_RAPIDITY_SIGMA,
    KAON_E_MAX,
)
from production.decay_engine.generate_meson_csvs import (
    _init_hnlcalc, compute_production_br_components, _sample_hnl_from_mesons,
    _write_csv, _write_empty_csv,
)

OUTPUT_BASE = PROJECT_ROOT / "output" / "llp_4vectors"

N_POOL = 100_000
KAON_PDG = 321


def sample_kaon_4vectors(n_pool, rng):
    """Sample charged-kaon 4-vectors from a parametrized soft-QCD flux.

    Transverse momentum follows a Tsallis/Hagedorn shape
    ``pdf(pT) ∝ pT * (1 + (mT - m_K)/(n T))**(-n)`` (inverse-CDF sampled), and
    rapidity is Gaussian about mid-rapidity. This is an approximate stand-in for
    a measured K± spectrum, not a precision input.

    Returns a dict with the same keys as
    ``meson_sampler.sample_meson_4vectors`` so the decay code is reused as-is.
    """
    m_k = M_KAON

    # Inverse-CDF sampling of the Tsallis pT pdf on a fine grid.
    n_grid = 4000
    pt_edges = np.linspace(0.0, KAON_PT_MAX, n_grid + 1)
    pt_centers = 0.5 * (pt_edges[:-1] + pt_edges[1:])
    mt = np.sqrt(pt_centers**2 + m_k**2)
    pdf = pt_centers * (1.0 + (mt - m_k) / (KAON_TSALLIS_N * KAON_TSALLIS_T)) ** (-KAON_TSALLIS_N)
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]

    u = rng.random(n_pool)
    idx = np.searchsorted(cdf, u)
    idx = np.clip(idx, 0, n_grid - 1)
    # Uniform sub-bin smearing inside the drawn pT bin.
    pt = rng.uniform(pt_edges[idx], pt_edges[idx + 1])

    # Sample rapidity, then accept-reject on the physical energy ceiling
    # E_K = mT * cosh(y) < KAON_E_MAX. A plain Gaussian tail extends to
    # |y| ~ 5 sigma where E_K reaches 30-70 TeV and the 2-body boost numerics
    # break down (events with m_HNL^2 < 0). The per-event accept-reject is
    # tighter than a static |y| cap because the limit depends on pT.
    y = rng.normal(0.0, KAON_RAPIDITY_SIGMA, n_pool)
    mt = np.sqrt(pt * pt + m_k * m_k)
    bad = mt * np.cosh(y) >= KAON_E_MAX
    while bad.any():
        y[bad] = rng.normal(0.0, KAON_RAPIDITY_SIGMA, int(bad.sum()))
        bad = mt * np.cosh(y) >= KAON_E_MAX
    phi = rng.uniform(0.0, 2.0 * np.pi, n_pool)

    mt_samp = np.sqrt(pt**2 + m_k**2)
    pz = mt_samp * np.sinh(y)
    E = mt_samp * np.cosh(y)
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)

    return {
        'E': E, 'px': px, 'py': py, 'pz': pz,
        'species_pdg': np.full(n_pool, KAON_PDG, dtype=int),
        'pt': pt, 'y': y, 'phi': phi,
    }


def process_flavor(flavor, pool, masses, rng):
    """Decay the kaon pool into N+X for one flavor across all mass points."""
    hnl = _init_hnlcalc(flavor)
    lepton_pdg = FLAVOR_TO_LEPTON_PDG[flavor]
    m_lepton = LEPTON_MASSES[flavor]
    out_dir = OUTPUT_BASE / flavor / "Kmeson"
    out_dir.mkdir(parents=True, exist_ok=True)

    n_pool = len(pool['E'])
    m_parent = M_KAON

    for m_N in masses:
        mass_label = format_mass_for_filename(m_N)
        csv_path = out_dir / f"mN_{mass_label}.csv"

        # Whole channel closes once even the 2-body mode is kinematically shut.
        if m_N >= m_parent - m_lepton:
            _write_empty_csv(csv_path)
            continue

        _, br_3body_channels, br = compute_production_br_components(
            hnl, KAON_PDG, lepton_pdg, m_N
        )
        if br <= 0:
            _write_empty_csv(csv_path)
            continue

        hnl_4v = _sample_hnl_from_mesons(
            pool['E'], pool['px'], pool['py'], pool['pz'],
            m_parent, m_lepton, m_N,
            br_3body_channels, br, hnl, rng,
        )
        w = SIGMA_KAON_PB * br / n_pool
        weights = np.full(n_pool, w)
        _write_csv(csv_path, weights, hnl_4v[:, 0], hnl_4v[:, 1], hnl_4v[:, 2], hnl_4v[:, 3])
        print(f"    {csv_path.name}: {n_pool} events, BR={br:.3e}, w_sum={weights.sum():.3e}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Generate kaon -> HNL CSVs (K+ -> l N)")
    parser.add_argument("--flavor", choices=["Ue", "Umu", "Utau"], nargs="+",
                        default=["Ue", "Umu", "Utau"])
    parser.add_argument("--n-pool", type=int, default=N_POOL)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--masses", type=float, nargs="+", default=None,
                        help="Custom mass list (default: full grid)")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    random.seed(args.seed)  # HNLCalc's 3-body BR integrator uses stdlib random
    masses = args.masses if args.masses else MASS_GRID

    print(f"Sampling charged-kaon pool ({args.n_pool} events)...")
    pool = sample_kaon_4vectors(args.n_pool, rng)

    for flavor in args.flavor:
        print(f"\n{'='*60}\nFlavor: {flavor}\n{'='*60}")
        process_flavor(flavor, pool, masses, rng)

    print("\nDone.")


if __name__ == "__main__":
    main()
