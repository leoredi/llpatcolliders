from pathlib import Path

import numpy as np
import pandas as pd

from dark_photon.plot_exclusion import EXISTING_CONSTRAINTS, _segments


def test_senscalc_existing_exclusion_is_self_describing():
    assert Path(EXISTING_CONSTRAINTS).exists()
    data = np.load(EXISTING_CONSTRAINTS, allow_pickle=False)
    assert int(data["n_polygons"]) == 5
    assert str(data["mass_units"]) == "GeV"
    assert str(data["ordinate"]) == "epsilon"
    assert str(data["source_revision"]) == (
        "0bca050633aae16e148d47f21840fa07ff4b8724"
    )
    for index in range(5):
        points = data[f"polygon_{index}"]
        assert points.ndim == 2 and points.shape[1] == 2
        assert np.isfinite(points).all()
        assert (points > 0).all()


def test_plot_omits_unresolved_single_mass_islands():
    frame = pd.DataFrame({
        "mass_GeV": [0.150, 0.151, 0.152, 0.153, 0.154],
        "has_sensitivity": [True, True, False, True, False],
    })
    segments = list(_segments(frame))
    assert len(segments) == 1
    assert segments[0]["mass_GeV"].tolist() == [0.150, 0.151]
