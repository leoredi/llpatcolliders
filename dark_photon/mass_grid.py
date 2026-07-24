"""Canonical BC1 sensitivity and Drell--Yan validation grids [GeV].

The central sensitivity scan uses 1 MeV spacing from 20 to 200 MeV.  The
completed extended campaign shows that all perturbative Drell--Yan points are
many orders of magnitude below the three-event threshold, so those points are
retained as a production-validation grid rather than diluting the canonical
physics scan and plot.
"""
from __future__ import annotations

import numpy as np


MESON_MASS_GRID = [milligev / 1000 for milligev in range(20, 201)]

DY_MASS_GRID = sorted({round(x, 3) for x in np.concatenate([
    np.arange(1.650, 3.001, 0.050),    # low-mass DY turn-on
    np.arange(3.100, 10.001, 0.100),
    np.arange(3.250, 10.001, 0.250),   # retain original quarter-GeV anchors
    np.arange(10.250, 50.001, 0.250),
])})

MASS_GRID = MESON_MASS_GRID

assert len(MESON_MASS_GRID) == 181
assert len(DY_MASS_GRID) == 272
assert len(MASS_GRID) == 181
