import os

from dark_photon import production
from dark_photon.run_sensitivity import _production_is_stale


def test_empty_vector_cache_is_stale(tmp_path):
    csv = tmp_path / "mA_1p700.csv"
    csv.touch()
    assert _production_is_stale(1.7, csv)


def test_dy_vector_older_than_spectrum_is_stale(tmp_path, monkeypatch):
    spectrum_dir = tmp_path / "spectra"
    spectrum_dir.mkdir()
    spectrum = spectrum_dir / "dy_2p000.npz"
    spectrum.write_bytes(b"spectrum")
    csv = tmp_path / "mA_2p000.csv"
    csv.write_bytes(b"vectors")
    os.utime(csv, (1, 1))
    os.utime(spectrum, (2, 2))
    monkeypatch.setattr(production, "DY_SPECTRUM_DIR", spectrum_dir)
    assert _production_is_stale(2.0, csv)

    os.utime(csv, (3, 3))
    assert not _production_is_stale(2.0, csv)
