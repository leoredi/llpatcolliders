import numpy as np

from scalar.production import B_SPECIES, generate_scalar_4vectors


def test_scalar_production_propagates_fonll_importance_weights():
    sampling_weight = np.array([0.5, 1.0, 1.5])
    pool = {
        "pt": np.array([2.0, 4.0, 8.0]),
        "y": np.array([-0.2, 0.0, 0.3]),
        "phi": np.array([0.1, 1.2, 2.4]),
        "sampling_weight": sampling_weight,
    }

    weights, *_ = generate_scalar_4vectors(
        0.5,
        len(sampling_weight),
        np.random.default_rng(7),
        sigma_bottom=100.0,
        pool=pool,
    )

    by_species = weights.reshape(len(B_SPECIES), len(sampling_weight))
    expected_shape = sampling_weight / sampling_weight.mean()
    for row in by_species:
        assert np.allclose(row / row.mean(), expected_shape)
