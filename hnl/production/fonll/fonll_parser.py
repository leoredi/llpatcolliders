"""
production/fonll/fonll_parser.py

Parse vendored FONLL meson-level differential cross-section tables.

Files are rectangular grids in (pT, y), 3 columns:
pT  y  dσ/dpT/dy [pb/GeV].  The current NNPDF4.0 central tables use
100 pT nodes over 0..50 GeV and 100 y nodes over -3..3.
pT varies slowly (outer loop), y varies fast (inner loop).
"""

import numpy as np
import os
from pathlib import Path

# Vendored FONLL table paths (relative to project root)
_VENDORED_DIR = Path(__file__).parent.parent.parent / "vendored"

FONLL_FILE_SETS = {
    "nnpdf40_nlo": {
        "bottom": _VENDORED_DIR / "fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_bottom.dat",
        "charm": _VENDORED_DIR / "fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_charm.dat",
    },
    "cteq66_legacy": {
        "bottom": _VENDORED_DIR / "fonll_pp14tev_cteq66_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_bottom.dat",
        "charm": _VENDORED_DIR / "fonll_pp14tev_cteq66_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_charm.dat",
    },
}

# HNL_FONLL_SET is read once at module import. To switch backends, export
# it in the shell before starting Python; mutating os.environ at runtime
# has no effect on FONLL_FILES.
FONLL_DEFAULT_SET = os.environ.get("HNL_FONLL_SET", "nnpdf40_nlo")
if FONLL_DEFAULT_SET not in FONLL_FILE_SETS:
    valid = ", ".join(sorted(FONLL_FILE_SETS))
    raise ValueError(f"unknown HNL_FONLL_SET={FONLL_DEFAULT_SET!r}; valid values: {valid}")

FONLL_FILES = FONLL_FILE_SETS[FONLL_DEFAULT_SET]


def parse_fonll_file(path):
    """
    Parse a FONLL meson-level dσ/dpT/dy table.

    Tables are required to be ordered with pT outermost (slow-varying) and
    y innermost (fast-varying), with y in ascending order, and to be written
    with a single consistent decimal format so that repeated pT and y values
    round-trip as bit-identical floats through ``np.loadtxt``. The parser
    relies on exact equality of those coordinates against ``np.unique`` for
    the rectangularity check; tolerance-aware clustering is intentionally
    not attempted.

    Parameters
    ----------
    path : str or Path
        Path to the .dat file.

    Returns
    -------
    pt_arr : ndarray, shape (N_pt,)
        Unique pT values in GeV.
    y_arr : ndarray, shape (N_y,)
        Unique rapidity values.
    dsigma_2d : ndarray, shape (N_pt, N_y)
        dσ/dpT/dy in pb/GeV.
    """
    data = np.loadtxt(path, comments="#")
    pt_all = data[:, 0]
    y_all = data[:, 1]
    dsigma_all = data[:, 2]

    pt_unique = np.unique(pt_all)
    y_unique = np.unique(y_all)
    n_pt = len(pt_unique)
    n_y = len(y_unique)

    # Reshape: pT varies slowly (outer), y varies fast (inner)
    if len(data) != n_pt * n_y:
        raise ValueError(f"{path}: expected a rectangular pT-y grid")
    expected_pt = np.repeat(pt_unique, n_y)
    expected_y = np.tile(y_unique, n_pt)
    if not np.array_equal(pt_all, expected_pt) or not np.array_equal(y_all, expected_y):
        # Distinguish a genuinely descending-y file from a shuffled one: the
        # first n_y rows match y_unique reversed iff y is strictly descending
        # in the inner loop. The shuffled case (e.g. y-outer / pT-inner) does
        # not match either, so it falls through to the generic message.
        if np.array_equal(y_all[:n_y], y_unique[::-1]):
            raise ValueError(f"{path}: y column must be in ascending order")
        raise ValueError(f"{path}: rows must be ordered with pT outermost and y innermost")
    dsigma_2d = dsigma_all.reshape(n_pt, n_y)

    return pt_unique, y_unique, dsigma_2d


def get_sigma_total(quark):
    """
    Integrate FONLL dσ/dpT/dy over the full (pT, y) grid → total σ in pb.

    Uses the trapezoidal rule on the 2D grid.

    Parameters
    ----------
    quark : str
        "bottom" or "charm"

    Returns
    -------
    float
        Total meson-level cross-section in pb.
    """
    path = FONLL_FILES[quark]
    pt_arr, y_arr, dsigma_2d = parse_fonll_file(path)

    # Trapezoidal integration: first over y, then over pT.
    # np.trapz was renamed np.trapezoid in NumPy 2.0; support both.
    trapezoid = getattr(np, "trapezoid", None) or np.trapz
    integral_over_y = trapezoid(dsigma_2d, y_arr, axis=1)  # shape (N_pt,)
    sigma_total = trapezoid(integral_over_y, pt_arr)        # scalar (pb)

    return sigma_total
