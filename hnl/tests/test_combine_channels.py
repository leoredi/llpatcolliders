"""Regression tests for production.combine_channels.

After the tau split (Dec 2026 refactor), the channel list must include
both `tau` (prompt) and `induced_tau` (from D/B meson decays). The combine
step must remain silent-skip on missing CSVs so kinematically-closed
points (Utau Dmeson at any mass, m_N > m_tau for both tau channels) do
not break it.
"""

import numpy as np
import pytest

from production import combine_channels


def test_channels_includes_both_tau_variants():
    """tau and induced_tau both in the default CHANNELS list."""
    assert "tau" in combine_channels.CHANNELS
    assert "induced_tau" in combine_channels.CHANNELS
    # Order doesn't matter to combine semantics, but the canonical inventory
    # should still be present.
    for ch in ("Bmeson", "Dmeson", "Bc", "Kmeson", "WZ"):
        assert ch in combine_channels.CHANNELS


def test_combine_handles_missing_channels(tmp_path, monkeypatch):
    """A missing channel directory contributes zero rows, no exception."""
    monkeypatch.setattr(combine_channels, "OUTPUT_BASE", tmp_path)
    flavor = "Ue"
    mass = 1.0
    mass_label = "1p000"

    # Create only ONE channel directory with a small CSV; the rest are absent.
    ch_dir = tmp_path / flavor / "Bmeson"
    ch_dir.mkdir(parents=True)
    sample = np.array([[1.0, 50.0, 1.0, 2.0, 3.0],
                       [2.0, 60.0, 4.0, 5.0, 6.0]])
    np.savetxt(ch_dir / f"mN_{mass_label}.csv", sample, delimiter=",", fmt="%.8e")

    n = combine_channels.combine_for_point(flavor, mass)
    assert n == 2

    combined_path = tmp_path / flavor / "combined" / f"mN_{mass_label}.csv"
    assert combined_path.exists()
    out = np.loadtxt(combined_path, delimiter=",")
    if out.ndim == 1:
        out = out.reshape(1, -1)
    assert out.shape == (2, 5)


def test_combine_empty_channel_csv_is_skipped(tmp_path, monkeypatch):
    """Zero-byte CSVs (kinematically-closed channels) are silently skipped."""
    monkeypatch.setattr(combine_channels, "OUTPUT_BASE", tmp_path)
    flavor = "Utau"
    mass = 0.3
    mass_label = "0p300"

    # Dmeson is empty (D -> tau N closed for Utau on the full grid).
    (tmp_path / flavor / "Dmeson").mkdir(parents=True)
    (tmp_path / flavor / "Dmeson" / f"mN_{mass_label}.csv").write_text("")

    # Bmeson has real rows.
    bmeson = tmp_path / flavor / "Bmeson"
    bmeson.mkdir(parents=True)
    np.savetxt(bmeson / f"mN_{mass_label}.csv",
               np.array([[7.0, 80.0, 0.0, 0.0, 70.0]]),
               delimiter=",", fmt="%.8e")

    n = combine_channels.combine_for_point(flavor, mass)
    assert n == 1
