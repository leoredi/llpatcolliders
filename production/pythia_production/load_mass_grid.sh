#!/bin/bash
# Mass Grid Loader for Production Scripts
# Sources mass grids from central config_mass_grid.py

# Find project root (two levels up from production/pythia_production/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Python executable (use conda env if available)
if [ -n "$CONDA_DEFAULT_ENV" ] && [ "$CONDA_DEFAULT_ENV" = "llpatcolliders" ]; then
    PYTHON="python"
else
    PYTHON="/opt/homebrew/Caskroom/miniconda/base/envs/llpatcolliders/bin/python"
fi

# Helper function to load mass grid from Python config
load_mass_grid() {
    local flavour=$1
    local mode=${2:-meson}  # Default to meson mode

    # Capture stdout only, redirect stderr to /dev/null to avoid polluting output
    local masses_str=$(cd "$PROJECT_ROOT" && $PYTHON -c "
from config_mass_grid import get_mass_grid
masses = get_mass_grid('$flavour', '$mode')
print(' '.join(f'{m:.2f}' for m in masses))
" 2>/dev/null)

    # Validate output contains only numbers, dots, and spaces
    if [[ -z "$masses_str" ]]; then
        echo "ERROR: Empty mass grid returned for $flavour $mode" >&2
        echo "()"
        return 1
    fi

    if [[ ! "$masses_str" =~ ^[0-9.\ ]+$ ]]; then
        echo "ERROR: Invalid mass grid format for $flavour $mode: $masses_str" >&2
        echo "()"
        return 1
    fi

    echo "($masses_str)"
}

# Export mass grids as arrays
# Usage in production scripts:
#   source load_mass_grid.sh
#   eval "ELECTRON_MASSES=$(load_mass_grid electron meson)"

# For convenience, also provide ready-to-use variables:
eval "ELECTRON_MASSES_MESON=$(load_mass_grid electron meson)"
eval "MUON_MASSES_MESON=$(load_mass_grid muon meson)"
eval "TAU_MASSES_MESON=$(load_mass_grid tau meson)"

# Validate arrays were loaded successfully
if [ ${#ELECTRON_MASSES_MESON[@]} -eq 0 ]; then
    echo "FATAL: Failed to load ELECTRON_MASSES_MESON from config_mass_grid.py" >&2
    exit 1
fi
if [ ${#MUON_MASSES_MESON[@]} -eq 0 ]; then
    echo "FATAL: Failed to load MUON_MASSES_MESON from config_mass_grid.py" >&2
    exit 1
fi
if [ ${#TAU_MASSES_MESON[@]} -eq 0 ]; then
    echo "FATAL: Failed to load TAU_MASSES_MESON from config_mass_grid.py" >&2
    exit 1
fi

# Print loaded grids (for debugging)
if [ "${VERBOSE_MASS_GRID:-0}" -eq 1 ]; then
    echo "Loaded mass grids from: $PROJECT_ROOT/config_mass_grid.py"
    echo "  ELECTRON_MASSES_MESON: ${#ELECTRON_MASSES_MESON[@]} points"
    echo "  MUON_MASSES_MESON: ${#MUON_MASSES_MESON[@]} points"
    echo "  TAU_MASSES_MESON: ${#TAU_MASSES_MESON[@]} points"
    echo ""
fi
