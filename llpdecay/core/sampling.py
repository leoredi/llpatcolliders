"""
Sampling utilities for Monte Carlo event generation.

Provides random number generation with proper seeding and common
sampling distributions.
"""

import numpy as np
from typing import Optional, List


def create_rng(seed: Optional[int] = None) -> np.random.Generator:
    """
    Create a random number generator.

    Parameters
    ----------
    seed : int, optional
        Random seed for reproducibility

    Returns
    -------
    np.random.Generator
        NumPy random generator
    """
    return np.random.default_rng(seed)


def sample_discrete(
    probabilities: List[float],
    n_samples: int = 1,
    rng: Optional[np.random.Generator] = None
) -> np.ndarray:
    """
    Sample from discrete probability distribution.

    Parameters
    ----------
    probabilities : list of float
        Probability for each outcome (will be normalized)
    n_samples : int
        Number of samples to generate
    rng : np.random.Generator, optional
        Random number generator

    Returns
    -------
    np.ndarray, shape (n_samples,)
        Sampled indices

    Examples
    --------
    >>> probs = [0.5, 0.3, 0.2]  # Three outcomes
    >>> indices = sample_discrete(probs, n_samples=100)
    >>> assert all(0 <= i < 3 for i in indices)
    """
    if rng is None:
        rng = np.random.default_rng()

    # Normalize probabilities
    probs = np.array(probabilities)
    probs = probs / np.sum(probs)

    return rng.choice(len(probs), size=n_samples, p=probs)


def rejection_sampling(
    distribution_func,
    x_range: tuple,
    n_samples: int,
    max_value: Optional[float] = None,
    rng: Optional[np.random.Generator] = None
) -> np.ndarray:
    """
    Sample from arbitrary 1D distribution using rejection sampling.

    Parameters
    ----------
    distribution_func : callable
        Function f(x) proportional to the target distribution
    x_range : tuple of float
        (x_min, x_max) sampling range
    n_samples : int
        Number of samples to generate
    max_value : float, optional
        Maximum value of distribution_func in x_range.
        If None, will be estimated.
    rng : np.random.Generator, optional
        Random number generator

    Returns
    -------
    np.ndarray, shape (n_samples,)
        Sampled values

    Examples
    --------
    >>> # Sample from f(x) = x² on [0, 1]
    >>> samples = rejection_sampling(lambda x: x**2, (0, 1), n_samples=1000)
    >>> assert all(0 <= x <= 1 for x in samples)
    """
    if rng is None:
        rng = np.random.default_rng()

    x_min, x_max = x_range

    # Estimate max if not provided
    if max_value is None:
        x_test = np.linspace(x_min, x_max, 1000)
        max_value = np.max(distribution_func(x_test)) * 1.1  # 10% safety margin

    samples = []
    while len(samples) < n_samples:
        # Propose random x
        x_prop = rng.uniform(x_min, x_max)
        # Accept with probability f(x) / f_max
        if rng.uniform(0, max_value) < distribution_func(x_prop):
            samples.append(x_prop)

    return np.array(samples)


def weighted_sample(
    values: np.ndarray,
    weights: np.ndarray,
    n_samples: int,
    replace: bool = True,
    rng: Optional[np.random.Generator] = None
) -> np.ndarray:
    """
    Sample from array with weights.

    Parameters
    ----------
    values : np.ndarray
        Array of values to sample from
    weights : np.ndarray
        Weight for each value
    n_samples : int
        Number of samples
    replace : bool
        Sample with replacement
    rng : np.random.Generator, optional
        Random number generator

    Returns
    -------
    np.ndarray
        Sampled values
    """
    if rng is None:
        rng = np.random.default_rng()

    # Normalize weights
    weights = np.array(weights)
    weights = weights / np.sum(weights)

    indices = rng.choice(len(values), size=n_samples, replace=replace, p=weights)
    return values[indices]


def uniform_sphere(n_samples: int, rng: Optional[np.random.Generator] = None) -> np.ndarray:
    """
    Sample points uniformly on unit sphere.

    Uses Marsaglia's method for uniform sampling.

    Parameters
    ----------
    n_samples : int
        Number of points
    rng : np.random.Generator, optional
        Random number generator

    Returns
    -------
    np.ndarray, shape (n_samples, 3)
        Unit vectors (x, y, z)

    Examples
    --------
    >>> points = uniform_sphere(100)
    >>> assert np.allclose(np.linalg.norm(points, axis=1), 1.0)
    """
    if rng is None:
        rng = np.random.default_rng()

    # Marsaglia method
    phi = rng.uniform(0, 2*np.pi, n_samples)
    cos_theta = rng.uniform(-1, 1, n_samples)
    sin_theta = np.sqrt(1 - cos_theta**2)

    x = sin_theta * np.cos(phi)
    y = sin_theta * np.sin(phi)
    z = cos_theta

    return np.stack([x, y, z], axis=-1)
