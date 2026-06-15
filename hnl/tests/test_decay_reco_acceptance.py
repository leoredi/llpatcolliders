"""Unit tests for the FairShip->GRENDEL signal acceptance core (Stage 2)."""
import sys
from pathlib import Path

import numpy as np
import pytest

_HNL = Path(__file__).resolve().parent.parent
for _p in (_HNL / "analysis", _HNL / "geometry"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import decay_reco_acceptance as dra  # noqa: E402
from grendel_geometry import mesh_fiducial  # noqa: E402


def test_boost_rest_to_lab_at_rest_is_identity():
    px = np.array([0.3, -0.3]); py = np.array([0.1, -0.1])
    pz = np.array([0.2, -0.2]); E = np.sqrt(px**2 + py**2 + pz**2 + 0.01)
    lx, ly, lz, lE = dra.boost_rest_to_lab([1.0, 0, 0, 0], px, py, pz, E)
    assert np.allclose(lx, px) and np.allclose(ly, py)
    assert np.allclose(lz, pz) and np.allclose(lE, E)


def test_boost_rest_to_lab_increases_energy_and_conserves_momentum():
    # consistent 2-body decay of a mass-1 parent: back-to-back, each E = M/2,
    # so the rest-frame daughter energies sum to the parent mass.
    px = np.array([0.5, -0.5]); py = np.array([0.0, 0.0]); pz = np.array([0.0, 0.0])
    E = np.array([0.5, 0.5])
    p4 = np.array([10.0, 0.0, 0.0, np.sqrt(100 - 1.0)])  # parent (m=1) moving in +z
    lx, ly, lz, lE = dra.boost_rest_to_lab(p4, px, py, pz, E)
    assert lE.sum() == pytest.approx(p4[0], rel=1e-9)          # energy conserved
    assert lz.sum() == pytest.approx(p4[3], rel=1e-9)          # pz conserved
    assert lE.sum() > E.sum()                                  # boosted up


def test_best_two_directions_picks_leading_charged():
    # 3 charged + 1 neutral; expect the two highest-|p| charged, unit dirs
    px = np.array([5.0, 0.0, 3.0, 0.0]); py = np.array([0.0, 4.0, 0.0, 9.0])
    pz = np.array([0.0, 0.0, 0.0, 0.0])
    charge = np.array([1.0, -1.0, 1.0, 0.0]); stable = np.array([1, 1, 1, 1])
    out = dra.best_two_directions(px, py, pz, charge, stable)
    assert out is not None
    d1, d2, p_soft = out
    assert np.allclose(np.linalg.norm(d1), 1.0)
    assert np.allclose(d1, [1, 0, 0])          # |p|=5 leading
    assert np.allclose(d2, [0, 1, 0])          # |p|=4 second (neutral 9 excluded)
    assert p_soft == pytest.approx(4.0)


def test_best_two_directions_needs_two_above_cut():
    px = np.array([5.0, 0.1]); py = np.zeros(2); pz = np.zeros(2)
    charge = np.array([1.0, -1.0]); stable = np.array([1, 1])
    assert dra.best_two_directions(px, py, pz, charge, stable) is None  # 2nd below P_CUT


def _hitting_vertex(rng):
    origin = np.zeros(3)
    for _ in range(5000):
        d = rng.normal(size=3); d[1] = abs(d[1]) + 0.5; d /= np.linalg.norm(d)
        loc, _, _ = mesh_fiducial.ray.intersects_location([origin], [d])
        if len(loc) >= 2:
            ds = np.sort(np.linalg.norm(loc - origin, axis=1))
            return origin + 0.5 * (ds[0] + ds[1]) * d, d
    pytest.skip("no hitting direction found")


def test_reconstruct_and_select_smoke():
    rng = np.random.default_rng(1)
    vtx, d = _hitting_vertex(rng)
    # two forward daughters with a small opening about the HNL direction
    e1 = np.cross(d, [0, 0, 1.0]); e1 /= np.linalg.norm(e1)
    dir1 = d + 0.05 * e1; dir1 /= np.linalg.norm(dir1)
    dir2 = d - 0.05 * e1; dir2 /= np.linalg.norm(dir2)
    N = 16
    mc = dra.reconstruct_decays(
        np.tile(vtx, (N, 1)), np.tile(dir1, (N, 1)), np.tile(dir2, (N, 1)),
        np.full(N, 5.0), dra.HIT_RESOLUTION, dra.SIGMA_T_DEFAULT, rng)
    for k in ('sep', 'sep_outer', 'open_angle', 'dca', 'collin', 'pointing',
              'vtx_in', 'on_tracker', 'timing_chi2', 'p_soft'):
        assert k in mc and len(mc[k]) == N
    sel = dra.selection_mask(mc)
    assert sel.dtype == bool and len(sel) == N
