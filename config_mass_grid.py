#!/usr/bin/env python3
"""
HNL Mass Grid Configuration
============================

SINGLE SOURCE OF TRUTH for all mass points.

To add/modify mass points: Edit the lists below.
Meson (Pythia) and EW (MadGraph) grids are kept side-by-side and can be
combined when both are needed.

Adaptive spacing optimized for physics:
- Dense near kinematic thresholds (K, D, B, tau masses)
"""

# ===========================================================================
# KINEMATIC LIMITS FOR MESON DECAYS
# ===========================================================================

# Kinematic limits for meson decays M → ℓ N (2-body)
# m_N < m_M - m_ℓ
KINEMATIC_LIMITS = {
    'K': {'electron': 0.493 - 0.000511, 'muon': 0.493 - 0.106},  # ~0.493, ~0.387 GeV
    'D': {'electron': 1.870 - 0.000511, 'muon': 1.870 - 0.106},  # ~1.87, ~1.76 GeV
    'B': {'electron': 5.279 - 0.000511, 'muon': 5.279 - 0.106, 'tau': 5.279 - 1.777},  # ~5.28, ~5.17, ~3.50 GeV
}

# ===========================================================================
# MESON (PYTHIA) MASS GRIDS  (validated, low-mass <~8 GeV)
# ===========================================================================

# Dense, inclusive mass grid for all regimes (meson + EW):
# Union of legacy + new production points to avoid sensitivity gaps
#
# Mass grid usage by meson type:
#   K mesons: m_N < 0.39 GeV (muon), < 0.49 GeV (electron)
#   D mesons: m_N < 1.76 GeV (muon), < 1.87 GeV (electron)
#   B mesons: m_N < 5.17 GeV (muon), < 5.28 GeV (electron), < 3.50 GeV (tau)
# The C++ Pythia code handles kinematic filtering; mass points outside valid ranges produce zero events.
_COMMON_GRID = [
    # K-regime (0.2-0.5 GeV): 0.05 GeV steps near kaon threshold
    0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,

    # D-regime (0.5-2.0 GeV): 0.1-0.2 GeV steps through charm transition
    0.60, 0.70, 0.80, 0.90, 1.00,
    1.10, 1.20, 1.25, 1.30, 1.40,
    1.50, 1.60, 1.70, 1.80, 1.90, 2.00,

    # B-regime (2.0-5.0 GeV): 0.2 GeV steps through beauty transition
    2.20, 2.30, 2.40, 2.60, 2.80, 3.00,
    3.20, 3.40, 3.60, 3.80, 4.00,
    4.20, 4.40, 4.60, 4.80, 5.00,

    # EW transition (5.0-8.0 GeV): 0.2 GeV steps through critical beauty→EW transition
    5.20, 5.40, 5.50, 5.60, 5.80, 6.00,
    6.20, 6.40, 6.60, 6.80, 7.00,
    7.20, 7.40, 7.60, 7.70, 7.80, 8.00,

    # High-mass EW (8-17 GeV): 0.5 GeV steps through tail
    8.50, 8.70, 9.00, 9.50, 10.00,
    10.50, 11.00, 11.50, 12.00, 12.50,
    13.00, 13.50, 14.00, 14.50, 15.00,
    15.50, 16.00, 16.50, 17.00,
]

# Meson grids now use the common grid
ELECTRON_MASSES_MESON = _COMMON_GRID
MUON_MASSES_MESON = _COMMON_GRID

# Tau: only B mesons can produce tau-coupled HNL
# Kinematic constraint: B → τ N requires m_N + m_τ < m_B, so m_N < m_B - m_τ ≈ 3.5 GeV
# D mesons are kinematically forbidden (m_D - m_τ < 0.2 GeV)
TAU_MASSES_MESON = [m for m in _COMMON_GRID if m < 3.5]

# ===========================================================================
# ELECTROWEAK (MADGRAPH) MASS GRIDS (W/Z-mediated, high-mass)
# ===========================================================================

# Electroweak production (W/Z bosons) via MadGraph uses the same grid
_EW_LOW_EDGE = _COMMON_GRID  # naming kept for backward compatibility; splits no longer used
_EW_CORE = []  # unused; full grid is _COMMON_GRID

ELECTRON_MASSES_EW = _EW_LOW_EDGE + _EW_CORE
MUON_MASSES_EW = _EW_LOW_EDGE + _EW_CORE

# Tau: EW production W → τ N
# Kinematic constraint: m_N + m_τ < m_W, so m_N < m_W - m_τ ≈ 78.5 GeV
# All masses in the common grid (up to 17 GeV) are kinematically allowed
TAU_MASSES_EW = _COMMON_GRID

# ===========================================================================
# COMBINED MASS GRIDS
# ===========================================================================

# Full mass grids (meson + EW)
ELECTRON_MASSES_COMBINED = sorted(set(ELECTRON_MASSES_MESON + ELECTRON_MASSES_EW))
MUON_MASSES_COMBINED = sorted(set(MUON_MASSES_MESON + MUON_MASSES_EW))
TAU_MASSES_COMBINED = sorted(set(TAU_MASSES_MESON + TAU_MASSES_EW))

# ===========================================================================
# HELPER FUNCTIONS
# ===========================================================================

def get_mass_grid(flavour, mode='production'):
    """
    Get mass grid for a specific flavour.

    Parameters:
    -----------
    flavour : str
        'electron', 'muon', or 'tau'
    mode : str
        'meson'      - Only meson-driven points (Pythia, low mass)
        'ew'         - Only electroweak points (MadGraph, high mass)
        'combined'   - Union of meson + ew
        Backward-compatible aliases:
          'production'/'base' -> meson

    Returns:
    --------
    list of float : Mass points in GeV
    """
    flavour = flavour.lower()
    mode = mode.lower()

    # Backward-compatible aliases
    if mode in ['production', 'base']:
        mode = 'meson'

    if flavour == 'electron':
        if mode == 'meson':
            return ELECTRON_MASSES_MESON
        elif mode == 'ew':
            return ELECTRON_MASSES_EW
        elif mode == 'combined':
            return ELECTRON_MASSES_COMBINED
        else:
            raise ValueError(f"Unknown mode: {mode}")

    elif flavour == 'muon':
        if mode == 'meson':
            return MUON_MASSES_MESON
        elif mode == 'ew':
            return MUON_MASSES_EW
        elif mode == 'combined':
            return MUON_MASSES_COMBINED
        else:
            raise ValueError(f"Unknown mode: {mode}")

    elif flavour == 'tau':
        if mode == 'meson':
            return TAU_MASSES_MESON
        elif mode == 'ew':
            return TAU_MASSES_EW
        elif mode == 'combined':
            return TAU_MASSES_COMBINED
        else:
            raise ValueError(f"Unknown mode: {mode}")

    else:
        raise ValueError(f"Unknown flavour: {flavour}")


def format_mass_for_filename(mass):
    """
    Convert mass (float) to filename format (e.g., 2.6 -> '2p60', 0.25 -> '0p25')

    Parameters:
    -----------
    mass : float
        Mass in GeV

    Returns:
    --------
    str : Mass formatted for filename with 2 decimal places (e.g., '2p60', '0p25', '15p00')
    """
    return f"{mass:.2f}".replace('.', 'p')


def export_to_bash(flavour, mode='production'):
    """
    Export mass grid as bash array string for use in shell scripts.

    Parameters:
    -----------
    flavour : str
        'electron', 'muon', or 'tau'
    mode : str
        'meson', 'ew', or 'combined'

    Returns:
    --------
    str : Bash array definition (e.g., "MASSES=(0.2 0.22 0.25 ...)")
    """
    masses = get_mass_grid(flavour, mode)
    mass_str = ' '.join(f"{m:.2f}" for m in masses)
    return f"MASSES=({mass_str})"


# ===========================================================================
# VALIDATION
# ===========================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("HNL MASS GRID CONFIGURATION")
    print("=" * 70)
    print()

    for flavour in ['electron', 'muon', 'tau']:
        print(f"{flavour.upper()} COUPLING:")
        print(f"  Meson (Pythia):     {len(get_mass_grid(flavour, 'meson'))} points")
        print(f"  EW (MadGraph):      {len(get_mass_grid(flavour, 'ew'))} points")
        print(f"  Combined:           {len(get_mass_grid(flavour, 'combined'))} points")
        print()

    print("=" * 70)
    print("EXAMPLE USAGE:")
    print("=" * 70)
    print()
    print("# Python:")
    print("from config_mass_grid import get_mass_grid")
    print("masses = get_mass_grid('electron', 'meson')")
    print()
    print("# Bash:")
    print("python -c 'from config_mass_grid import export_to_bash; print(export_to_bash(\"electron\", \"meson\"))'")
    print()
