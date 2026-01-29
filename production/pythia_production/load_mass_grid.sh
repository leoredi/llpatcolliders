#!/bin/bash
# Load MASS_GRID from config_mass_grid.py

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

if [ -n "$CONDA_DEFAULT_ENV" ] && [ "$CONDA_DEFAULT_ENV" = "llpatcolliders" ]; then
    PYTHON="python"
elif command -v conda >/dev/null 2>&1; then
    PYTHON="conda run -n llpatcolliders python"
else
    PYTHON="python3"
fi

MASS_GRID_STR=$(cd "$PROJECT_ROOT" && $PYTHON -c "
from config_mass_grid import MASS_GRID
print(' '.join(f'{m:.2f}' for m in MASS_GRID))
" 2>/dev/null)

eval "MASS_GRID=($MASS_GRID_STR)"

if [ ${#MASS_GRID[@]} -eq 0 ]; then
    echo "FATAL: Failed to load MASS_GRID from config_mass_grid.py" >&2
    exit 1
fi

[ "${VERBOSE_MASS_GRID:-0}" -eq 1 ] && echo "MASS_GRID: ${#MASS_GRID[@]} points" || true
