#!/usr/bin/env python3
"""
HNL Mass Grid Configuration
============================

SINGLE SOURCE OF TRUTH for all mass points.

To add/modify mass points: Edit the lists below.
Production script automatically uses PRODUCTION grids (base + closure).

Adaptive spacing optimized for physics:
- Dense near kinematic thresholds (K, D, B, tau masses)
- Closure points (4-5.5 GeV) ensure smooth exclusion island boundaries
"""

# ===========================================================================
# BASE MASS GRIDS (Low-mass regime < 5 GeV - VALIDATED)
# ===========================================================================

# Electron coupling (Benchmark 100: Ue²=x, Uμ²=0, Uτ²=0)
# Extended to 8 GeV for smooth meson→EW transition
ELECTRON_MASSES_BASE = [
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
MUON_MASSES_BASE = [
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
TAU_MASSES_BASE = [
    # D-meson regime (0.5-2.0 GeV): Tau threshold at 1.777 GeV
    0.50, 0.55, 0.60, 0.65, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30,
    # Near tau threshold (1.4-2.0 GeV): Dense sampling
    1.40, 1.45, 1.50, 1.55, 1.60, 1.62, 1.64, 1.66, 1.70, 1.74, 1.78, 1.80, 1.85, 1.90, 2.00,
    # B-meson regime (2.0-8.0 GeV): Extended for better transition coverage
    2.40, 2.80, 3.00, 3.20, 3.60, 4.00, 4.50, 5.00, 5.50, 6.0, 6.5, 7.0, 7.5, 8.0
]

# ===========================================================================
# ISLAND CLOSURE POINTS (4-5.5 GeV - NOW COVERED BY EXTENDED GRIDS)
# ===========================================================================

# No longer needed - transition region now covered by:
# - Pythia meson production extended to 8 GeV
# - MadGraph EW production extended down to 4 GeV

ELECTRON_MASSES_CLOSURE = []
MUON_MASSES_CLOSURE = []
TAU_MASSES_CLOSURE = []

# ===========================================================================
# HIGH-MASS REGIME (≥ 5 GeV - CURRENTLY FAILING IN PYTHIA)
# ===========================================================================

# Electroweak production (W/Z bosons) via MadGraph
# Extended DOWN to 4 GeV for smooth overlap with meson production

ELECTRON_MASSES_EW = [4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 25.0, 28.0, 30.0, 32.0, 35.0, 38.0, 40.0, 45.0, 50.0]
MUON_MASSES_EW = [4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 25.0, 28.0, 30.0, 32.0, 35.0, 38.0, 40.0, 45.0, 50.0]
TAU_MASSES_EW = [4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 25.0, 28.0, 30.0, 32.0, 35.0, 38.0, 40.0, 45.0, 50.0]

# ===========================================================================
# COMBINED MASS GRIDS
# ===========================================================================

# Full mass grids (base + closure + EW)
ELECTRON_MASSES_ALL = sorted(set(ELECTRON_MASSES_BASE + ELECTRON_MASSES_CLOSURE + ELECTRON_MASSES_EW))
MUON_MASSES_ALL = sorted(set(MUON_MASSES_BASE + MUON_MASSES_CLOSURE + MUON_MASSES_EW))
TAU_MASSES_ALL = sorted(set(TAU_MASSES_BASE + TAU_MASSES_CLOSURE + TAU_MASSES_EW))

# Recommended for production (base + closure only, skip failing EW regime)
ELECTRON_MASSES_PRODUCTION = sorted(set(ELECTRON_MASSES_BASE + ELECTRON_MASSES_CLOSURE))
MUON_MASSES_PRODUCTION = sorted(set(MUON_MASSES_BASE + MUON_MASSES_CLOSURE))
TAU_MASSES_PRODUCTION = sorted(set(TAU_MASSES_BASE + TAU_MASSES_CLOSURE))

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
        'base' - Only validated low-mass points (< 5 GeV)
        'closure' - Only island closure points (4-5.5 GeV)
        'production' - Base + closure (recommended for new simulations)
        'all' - Base + closure + EW (includes failing high-mass points)

    Returns:
    --------
    list of float : Mass points in GeV
    """
    flavour = flavour.lower()
    mode = mode.lower()

    if flavour == 'electron':
        if mode == 'base':
            return ELECTRON_MASSES_BASE
        elif mode == 'closure':
            return ELECTRON_MASSES_CLOSURE
        elif mode == 'production':
            return ELECTRON_MASSES_PRODUCTION
        elif mode == 'all':
            return ELECTRON_MASSES_ALL
        else:
            raise ValueError(f"Unknown mode: {mode}")

    elif flavour == 'muon':
        if mode == 'base':
            return MUON_MASSES_BASE
        elif mode == 'closure':
            return MUON_MASSES_CLOSURE
        elif mode == 'production':
            return MUON_MASSES_PRODUCTION
        elif mode == 'all':
            return MUON_MASSES_ALL
        else:
            raise ValueError(f"Unknown mode: {mode}")

    elif flavour == 'tau':
        if mode == 'base':
            return TAU_MASSES_BASE
        elif mode == 'closure':
            return TAU_MASSES_CLOSURE
        elif mode == 'production':
            return TAU_MASSES_PRODUCTION
        elif mode == 'all':
            return TAU_MASSES_ALL
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
        'base', 'closure', 'production', or 'all'

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
        print(f"  Base (< 5 GeV):     {len(get_mass_grid(flavour, 'base'))} points")
        print(f"  Closure (4-5.5):    {len(get_mass_grid(flavour, 'closure'))} points")
        print(f"  Production (rec):   {len(get_mass_grid(flavour, 'production'))} points")
        print(f"  All (inc. EW):      {len(get_mass_grid(flavour, 'all'))} points")
        print()

    print("=" * 70)
    print("EXAMPLE USAGE:")
    print("=" * 70)
    print()
    print("# Python:")
    print("from config_mass_grid import get_mass_grid")
    print("masses = get_mass_grid('electron', 'production')")
    print()
    print("# Bash:")
    print("python -c 'from config_mass_grid import export_to_bash; print(export_to_bash(\"electron\", \"production\"))'")
    print()
