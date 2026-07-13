"""Tests for the SensCalc LHC parent-spectrum boundary."""

import math

import numpy as np
import pytest

from alp_fermion.production_spectra import (
    TabulatedSpectrum,
    transverse_theta_interval,
)


def _write_grid(path, theta, energy, density):
    rows = []
    for i, theta_value in enumerate(theta):
        for j, energy_value in enumerate(energy):
            rows.append((theta_value, energy_value, density[i, j]))
    np.savetxt(path, rows)


def test_load_integrate_and_select_polar_interval(tmp_path):
    path = tmp_path / "grid.txt"
    theta = np.array([0.2, 1.0, 2.0, 2.8])
    energy = np.array([1.0, 2.0, 4.0])
    density = np.ones((len(theta), len(energy)))
    _write_grid(path, theta, energy, density)

    spectrum = TabulatedSpectrum.load(path)

    assert spectrum.integral() == pytest.approx((2.8 - 0.2) * (4.0 - 1.0))
    assert spectrum.polar_fraction(1.0, 2.0) == pytest.approx(1.0 / 2.6)


def test_sample_returns_on_shell_four_vectors(tmp_path):
    path = tmp_path / "grid.txt"
    theta = np.array([0.4, 1.2, 2.0, 2.7])
    energy = np.array([0.4, 2.0, 4.0])
    density = np.ones((len(theta), len(energy)))
    _write_grid(path, theta, energy, density)
    spectrum = TabulatedSpectrum.load(path)

    mass = 0.5
    sample = spectrum.sample(20_000, mass, rng=np.random.default_rng(1234))
    momentum_squared = sample["px"] ** 2 + sample["py"] ** 2 + sample["pz"] ** 2

    assert len(sample["E"]) == 20_000
    assert np.all(sample["E"] >= mass)
    assert np.allclose(sample["E"] ** 2 - momentum_squared, mass**2)
    assert np.all((sample["theta"] >= theta.min()) & (sample["theta"] <= theta.max()))


def test_load_rejects_incomplete_cartesian_grid(tmp_path):
    path = tmp_path / "bad.txt"
    np.savetxt(path, [(0.1, 1.0, 1.0), (0.2, 2.0, 1.0)])
    with pytest.raises(ValueError, match="not a Cartesian grid"):
        TabulatedSpectrum.load(path)


def test_transverse_theta_interval_matches_pseudorapidity_definition():
    theta_min, theta_max = transverse_theta_interval(0.5)
    assert theta_min == pytest.approx(2.0 * math.atan(math.exp(-0.5)))
    assert theta_max == pytest.approx(math.pi - theta_min)
