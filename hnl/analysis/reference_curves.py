"""
hnl/analysis/reference_curves.py

Load reference (m_N, U^2) exclusion contours from other PBC experiments
(MATHUSLA, ANUBIS, CODEX-b, SHiP, ...) for the money plot.

File format: space-separated, 3 columns (mass_GeV, u2_min, u2_max).
Lines starting with '#' are comments.

Files live in hnl/data/reference_curves/{experiment}_{flavor}.dat.
"""

import numpy as np
from pathlib import Path

# hnl/data/reference_curves/
DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "reference_curves"


def load_reference_curve(experiment, flavor):
    """Load a (m_N, U^2_min, U^2_max) reference contour for one (exp, flavor)."""
    path = DATA_DIR / f"{experiment}_{flavor}.dat"
    if not path.exists():
        return None
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] < 3:
        return None
    return {
        "mass": data[:, 0],
        "u2_min": data[:, 1],
        "u2_max": data[:, 2],
    }


def load_all_references(flavors=None):
    """Discover and load every reference curve in DATA_DIR."""
    if flavors is None:
        flavors = ["Ue", "Umu", "Utau"]

    if not DATA_DIR.exists():
        return {}

    experiments = set()
    for f in DATA_DIR.glob("*.dat"):
        name = f.stem  # e.g. "MATHUSLA_Ue"
        for flav in flavors:
            if name.endswith(f"_{flav}"):
                exp = name[: -(len(flav) + 1)]
                experiments.add(exp)

    result = {}
    for exp in sorted(experiments):
        curves = {}
        for flav in flavors:
            c = load_reference_curve(exp, flav)
            if c is not None:
                curves[flav] = c
        if curves:
            result[exp] = curves
    return result
