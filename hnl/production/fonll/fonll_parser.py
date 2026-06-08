"""Parse vendored FONLL meson-level differential cross-section tables."""

import numpy as np
from pathlib import Path

_VENDORED_DIR = Path(__file__).parent.parent.parent / "vendored"

FONLL_FILES = {
    "bottom": _VENDORED_DIR / "fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_bottom.dat",
    "charm": _VENDORED_DIR / "fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_dsdpTdy_pt0-50_y-3to3_central_charm.dat",
}


def parse_fonll_file(path):
    data = np.loadtxt(path, comments="#")
    pt_all = data[:, 0]
    y_all = data[:, 1]
    dsigma_all = data[:, 2]

    pt_unique = np.unique(pt_all)
    y_unique = np.unique(y_all)
    n_pt = len(pt_unique)
    n_y = len(y_unique)

    if len(data) != n_pt * n_y:
        raise ValueError(f"{path}: expected a rectangular pT-y grid")
    expected_pt = np.repeat(pt_unique, n_y)
    expected_y = np.tile(y_unique, n_pt)
    if not np.array_equal(pt_all, expected_pt) or not np.array_equal(y_all, expected_y):
        if np.array_equal(y_all[:n_y], y_unique[::-1]):
            raise ValueError(f"{path}: y column must be in ascending order")
        raise ValueError(f"{path}: rows must be ordered with pT outermost and y innermost")
    dsigma_2d = dsigma_all.reshape(n_pt, n_y)

    return pt_unique, y_unique, dsigma_2d


def get_sigma_total(quark):
    path = FONLL_FILES[quark]
    pt_arr, y_arr, dsigma_2d = parse_fonll_file(path)

    trapezoid = getattr(np, "trapezoid", None) or np.trapz
    integral_over_y = trapezoid(dsigma_2d, y_arr, axis=1)  # shape (N_pt,)
    sigma_total = trapezoid(integral_over_y, pt_arr)        # scalar (pb)

    return sigma_total
