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

# Electron coupling (Benchmark 100: Ue²=x, Uμ²=0, Uτ²=0)
# Extended to 8 GeV for smooth meson→EW transition
ELECTRON_MASSES_MESON = [
    # Kaon regime (0.2-0.5 GeV): Dense sampling near threshold
    0.20, 0.22, 0.25, 0.28, 0.30, 0.32, 0.35, 0.38, 0.40, 0.42, 0.45, 0.48,
    # D-meson regime (0.5-2.0 GeV): Peak sensitivity
    0.50, 0.52, 0.55, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40,
    # D-meson threshold region (1.5-2.0 GeV): Dense near mD
    1.50, 1.60, 1.70, 1.75, 1.80, 1.82, 1.85, 1.90, 2.00,
    # B-meson regime (2.0-8.0 GeV): Extended for better transition coverage
    2.30, 2.60, 3.00, 3.40, 3.80, 4.20, 4.60, 4.80, 5.00, 5.20, 5.50, 6.0, 6.5, 7.0, 7.5, 8.0
]

# Muon coupling (Benchmark 010: Ue²=0, Uμ²=x, Uτ²=0)
# Extended to 8 GeV for smooth meson→EW transition
MUON_MASSES_MESON = [
    # Kaon regime (0.2-0.5 GeV): Dense sampling near threshold
    0.20, 0.22, 0.25, 0.28, 0.30, 0.32, 0.35, 0.37, 0.38, 0.39, 0.40, 0.42, 0.45, 0.48,
    # D-meson regime (0.5-2.0 GeV): Peak sensitivity
    0.50, 0.55, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20, 1.40,
    # D-meson threshold region (1.6-2.0 GeV): Dense near mD
    1.60, 1.65, 1.70, 1.75, 1.76, 1.78, 1.80, 1.85, 1.90, 2.00,
    # B-meson regime (2.0-8.0 GeV): Extended for better transition coverage
    2.30, 2.60, 3.00, 3.40, 3.80, 4.20, 4.60, 4.80, 5.00, 5.20, 5.50, 6.0, 6.5, 7.0, 7.5, 8.0
]

# Tau coupling (Benchmark 001: Ue²=0, Uμ²=0, Uτ²=x)
# Extended to 8 GeV for smooth meson→EW transition
# Note: Below 1.64 GeV, tau simulations generate BOTH "_direct" and "_fromTau" files
TAU_MASSES_MESON = [
    # D-meson regime (0.5-2.0 GeV): Tau threshold at 1.777 GeV
    0.50, 0.55, 0.60, 0.65, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30,
    # Near tau threshold (1.4-2.0 GeV): Dense sampling
    1.40, 1.45, 1.50, 1.55, 1.60, 1.62, 1.64, 1.66, 1.70, 1.74, 1.78, 1.80, 1.85, 1.90, 2.00,
    # B-meson regime (2.0-8.0 GeV): Extended for better transition coverage
    2.40, 2.80, 3.00, 3.20, 3.60, 4.00, 4.50, 5.00, 5.50, 6.0, 6.5, 7.0, 7.5, 8.0
]

# ===========================================================================
# ELECTROWEAK (MADGRAPH) MASS GRIDS (W/Z-mediated, high-mass)
# ===========================================================================

# Electroweak production (W/Z bosons) via MadGraph
# Extended DOWN to 1 GeV (0.5 GeV steps) for smooth overlap with meson production, UP to 80 GeV

_EW_LOW_EDGE = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]  # overlap with meson
_EW_CORE = [5.0, 5.2, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]  # stop above 20 GeV

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
