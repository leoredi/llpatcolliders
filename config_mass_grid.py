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
# MESON (PYTHIA) MASS GRIDS  (validated, low-mass <~8 GeV)
# ===========================================================================

# Unified mass grid for all regimes (meson + EW):
# 0.2 → 1.0 GeV in 0.2 GeV steps; 1.0 → 11.0 GeV in 0.25 GeV steps; 12.0, 13.0 GeV.
_COMMON_GRID = (
    [round(x, 2) for x in [0.2 + 0.2 * i for i in range(0, 5)]] +
    [round(x, 2) for x in [1.0 + 0.25 * i for i in range(0, 41) if 1.0 + 0.25 * i <= 11.0]] +
    [12.0, 13.0]
)

# Meson grids now use the common grid
ELECTRON_MASSES_MESON = _COMMON_GRID
MUON_MASSES_MESON = _COMMON_GRID
TAU_MASSES_MESON = _COMMON_GRID

# ===========================================================================
# ELECTROWEAK (MADGRAPH) MASS GRIDS (W/Z-mediated, high-mass)
# ===========================================================================

# Electroweak production (W/Z bosons) via MadGraph uses the same grid
_EW_LOW_EDGE = _COMMON_GRID  # naming kept for backward compatibility; splits no longer used
_EW_CORE = []  # unused; full grid is _COMMON_GRID

ELECTRON_MASSES_EW = _EW_LOW_EDGE + _EW_CORE
MUON_MASSES_EW = _EW_LOW_EDGE + _EW_CORE
TAU_MASSES_EW = [4.0, 4.2, 4.4, 4.6, 4.8] + _EW_CORE

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
    Convert mass (float) to filename format (e.g., 2.6 -> '2p6', 0.25 -> '0p25')

    Parameters:
    -----------
    mass : float
        Mass in GeV

    Returns:
    --------
    str : Mass formatted for filename (e.g., '2p6', '0p25')
    """
    return f"{mass:.2f}".replace('.', 'p').rstrip('0').rstrip('p')


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
