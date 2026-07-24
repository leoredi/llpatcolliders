from pathlib import Path

import numpy as np

from dark_photon import production
from dark_photon.mass_grid import DY_MASS_GRID, MASS_GRID, MESON_MASS_GRID


SPECTRA = Path(__file__).resolve().parents[1] / "data" / "spectra" / "dy"


def test_pythia_dy_grid_is_complete_and_self_describing():
    paths = sorted(SPECTRA.glob("dy_*.npz"))
    assert len(paths) == len(DY_MASS_GRID) == 272
    rows = []
    seeds = []
    for path in paths:
        data = np.load(path)
        assert str(data["generator"]) == "Pythia 8.315 NewGaugeBoson pure-Zprime"
        assert str(data["pdf"]) == "NNPDF40_nlo_as_01180 member 0"
        assert float(data["histogram_coverage"]) == 1.0
        assert int(data["hist"].sum()) == int(data["n_sample"])
        assert float(data["sigma_dy_pb"]) > 0.0
        rows.append((
            float(data["mass_GeV"]),
            float(data["sigma_generator_pb"]),
            float(data["sigma_generator_err_pb"]),
        ))
        seeds.append(int(data["seed"]))
    rows.sort()
    assert [row[0] for row in rows] == DY_MASS_GRID
    assert rows[0][0] == 1.65 and rows[-1][0] == 50.0
    assert len(set(seeds)) == len(seeds)
    # The physical cross section decreases with mass.  Permit only fluctuations
    # compatible with five times the quoted integration uncertainty.
    for left, right in zip(rows, rows[1:]):
        tolerance = 5.0 * np.hypot(left[2], right[2])
        assert right[1] <= left[1] + tolerance


def test_canonical_mass_grid_is_one_mev_to_200_mev():
    assert len(MESON_MASS_GRID) == 181
    assert len(DY_MASS_GRID) == 272
    assert len(MASS_GRID) == 181
    assert MASS_GRID[0] == 0.02
    assert MASS_GRID[-1] == 0.2
    assert max(np.diff(MESON_MASS_GRID)) <= 0.0011
    assert max(np.diff([m for m in DY_MASS_GRID if m <= 3.0])) <= 0.0501


def test_dy_sampler_preserves_mass_and_cross_section():
    mass = 5.0
    p4, weights = production.generate_dy_channel(
        mass, np.random.default_rng(12345), spectrum_dir=SPECTRA
    )
    reconstructed = np.sqrt(
        np.maximum(p4[:, 0] ** 2 - np.sum(p4[:, 1:] ** 2, axis=1), 0.0)
    )
    data = np.load(SPECTRA / "dy_5p000.npz")
    assert np.allclose(reconstructed, mass, rtol=0.0, atol=2e-9)
    assert np.isclose(weights.sum(), float(data["sigma_dy_pb"]), rtol=1e-12)
