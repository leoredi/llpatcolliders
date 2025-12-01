#!/bin/bash
# Mass Grid Loader for Production Scripts
# Sources mass grids from central config_mass_grid.py

# Find project root (one level up from production/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Python executable (use conda env if available)
if [ -n "$CONDA_DEFAULT_ENV" ] && [ "$CONDA_DEFAULT_ENV" = "llpatcolliders" ]; then
    PYTHON="python"
else
    PYTHON="/opt/homebrew/Caskroom/miniconda/base/envs/llpatcolliders/bin/python"
fi

# Helper function to load mass grid from Python config
load_mass_grid() {
    local flavour=$1
    local mode=${2:-production}  # Default to 'production' mode

    # Get masses from Python config and convert to bash array
    local masses_str=$(cd "$PROJECT_ROOT" && $PYTHON -c "
from config_mass_grid import get_mass_grid
masses = get_mass_grid('$flavour', '$mode')
print(' '.join(f'{m:.2f}' for m in masses))
")

    # Convert space-separated string to array
    echo "($masses_str)"
}

# Export mass grids as arrays
# Usage in production scripts:
#   source load_mass_grid.sh
#   eval "ELECTRON_MASSES=$(load_mass_grid electron production)"

# For convenience, also provide ready-to-use variables:
eval "ELECTRON_MASSES_BASE=$(load_mass_grid electron base)"
eval "ELECTRON_MASSES_CLOSURE=$(load_mass_grid electron closure)"
eval "ELECTRON_MASSES_PRODUCTION=$(load_mass_grid electron production)"

eval "MUON_MASSES_BASE=$(load_mass_grid muon base)"
eval "MUON_MASSES_CLOSURE=$(load_mass_grid muon closure)"
eval "MUON_MASSES_PRODUCTION=$(load_mass_grid muon production)"

eval "TAU_MASSES_BASE=$(load_mass_grid tau base)"
eval "TAU_MASSES_CLOSURE=$(load_mass_grid tau closure)"
eval "TAU_MASSES_PRODUCTION=$(load_mass_grid tau production)"

# Print loaded grids (for debugging)
if [ "${VERBOSE_MASS_GRID:-0}" -eq 1 ]; then
    echo "Loaded mass grids from: $PROJECT_ROOT/config_mass_grid.py"
    echo "  ELECTRON_MASSES_BASE: ${#ELECTRON_MASSES_BASE[@]} points"
    echo "  ELECTRON_MASSES_CLOSURE: ${#ELECTRON_MASSES_CLOSURE[@]} points"
    echo "  ELECTRON_MASSES_PRODUCTION: ${#ELECTRON_MASSES_PRODUCTION[@]} points"
    echo "  MUON_MASSES_BASE: ${#MUON_MASSES_BASE[@]} points"
    echo "  MUON_MASSES_CLOSURE: ${#MUON_MASSES_CLOSURE[@]} points"
    echo "  MUON_MASSES_PRODUCTION: ${#MUON_MASSES_PRODUCTION[@]} points"
    echo "  TAU_MASSES_BASE: ${#TAU_MASSES_BASE[@]} points"
    echo "  TAU_MASSES_CLOSURE: ${#TAU_MASSES_CLOSURE[@]} points"
    echo "  TAU_MASSES_PRODUCTION: ${#TAU_MASSES_PRODUCTION[@]} points"
    echo ""
fi
