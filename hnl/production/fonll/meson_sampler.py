"""
production/fonll/meson_sampler.py

Inverse-CDF sampler from 2D FONLL dσ/dpT/dy grids.

Generates meson 4-vectors (E, px, py, pz) drawn from the FONLL
differential cross-section. Species are assigned over the simulated meson
subset only; absolute fragmentation fractions are applied later in the event
weights.
"""

import numpy as np
from .fonll_parser import parse_fonll_file, FONLL_FILES
from ..constants import QUARK_MESON_MAP, MESON_MASSES


def _node_bin_edges(values):
    """Build bin edges for the FONLL sampler.

    The convention is intentionally non-uniform at the table edges:
    interior grid nodes are bin centers (so the bin edges sit at the
    midpoints of consecutive nodes), but the first and last nodes are
    bin boundaries themselves. This keeps every sampled (pT, y) value
    strictly inside the tabulated grid range, with no overflow and no
    need for a post-hoc clip.
    """
    if len(values) < 2:
        raise ValueError("FONLL sampling requires at least two grid nodes")
    edges = np.empty(len(values) + 1)
    edges[0] = values[0]
    edges[-1] = values[-1]
    edges[1:-1] = 0.5 * (values[:-1] + values[1:])
    return edges


def _sample_node_intervals(edges, indices, rng):
    """Uniformly smear selected nodes inside their bounded bin intervals."""
    return rng.uniform(edges[indices], edges[indices + 1])


def _build_cdf(pt_arr, y_arr, dsigma_2d):
    """
    Build a 1D CDF from the 2D dσ/dpT/dy grid for inverse-CDF sampling.

    Each bin probability is dσ/dpT/dy × ΔpT × Δy (bin area).

    Returns
    -------
    cdf : ndarray, shape (N_pt * N_y,)
        Cumulative distribution function (normalized).
    pt_edges : ndarray
        Bin edges in pT.
    y_edges : ndarray
        Bin edges in rapidity.
    """
    n_pt = len(pt_arr)
    n_y = len(y_arr)

    # Bin scheme: interior grid nodes are bin centers, while the first and
    # last nodes are bin boundaries. Each bin's probability weight is the
    # bin width (np.diff of the edges) times dsigma at the node — half-width
    # for edge bins, full width for interior bins — so the sampler never
    # extends past the tabulated range.
    pt_edges = _node_bin_edges(pt_arr)
    y_edges = _node_bin_edges(y_arr)
    pt_widths = np.diff(pt_edges)
    y_widths = np.diff(y_edges)

    # 2D bin probabilities (unnormalized)
    prob_2d = dsigma_2d * pt_widths[:, None] * y_widths[None, :]
    prob_2d = np.maximum(prob_2d, 0.0)  # clip negatives

    # Flatten and build CDF
    prob_flat = prob_2d.ravel()
    cdf = np.cumsum(prob_flat)
    norm = cdf[-1]
    if not np.isfinite(norm) or norm <= 0:
        raise ValueError("FONLL sampling grid has no finite positive weight")
    cdf /= norm

    return cdf, pt_edges, y_edges


def sample_meson_4vectors(n_events, quark, rng=None, force_species=None):
    """
    Sample meson 4-vectors from the FONLL dσ/dpT/dy distribution.

    For each sample:
      1. Draw (pT, y) bin from inverse CDF
      2. Uniform sub-bin smearing within the drawn bin
      3. Uniform φ in [0, 2π)
      4. Assign meson species (from renormalized fragmentation fractions, or
         forced to a fixed PDG id via ``force_species``)
      5. Compute 4-vector (E, px, py, pz) from (pT, y, φ, m)

    Parameters
    ----------
    n_events : int
        Number of meson 4-vectors to generate.
    quark : str
        "bottom" or "charm"
    rng : numpy.random.Generator, optional
        Random number generator (for reproducibility).
    force_species : int, optional
        If set to a PDG id (e.g. 541 for Bc), every sampled event gets that
        PDG id and the mass used in the 4-vector reconstruction is
        ``MESON_MASSES[force_species]``. The (pT, y) sampling from the FONLL
        grid for the requested ``quark`` is unchanged. If ``None`` (default),
        species are drawn from the renormalized fragmentation cumulative.

    Returns
    -------
    dict with keys:
        'E', 'px', 'py', 'pz' : ndarray, shape (n_events,)
            Meson 4-vectors in the lab frame.
        'species_pdg' : ndarray of int, shape (n_events,)
            PDG ID of the assigned meson species.
        'pt', 'y', 'phi' : ndarray, shape (n_events,)
            Kinematic variables used to construct the 4-vector.
    """
    if rng is None:
        rng = np.random.default_rng()

    # Parse FONLL grid and build CDF
    path = FONLL_FILES[quark]
    pt_arr, y_arr, dsigma_2d = parse_fonll_file(path)
    cdf, pt_edges, y_edges = _build_cdf(pt_arr, y_arr, dsigma_2d)

    n_pt = len(pt_arr)
    n_y = len(y_arr)

    # Draw random samples from CDF
    u = rng.random(n_events)
    flat_idx = np.searchsorted(cdf, u)
    flat_idx = np.clip(flat_idx, 0, n_pt * n_y - 1)

    i_pt = flat_idx // n_y
    i_y = flat_idx % n_y

    # Sub-bin smearing: uniform within the bin
    pt_sampled = _sample_node_intervals(pt_edges, i_pt, rng)
    y_sampled = _sample_node_intervals(y_edges, i_y, rng)

    # Uniform azimuthal angle
    phi = rng.uniform(0, 2 * np.pi, n_events)

    # Assign meson species: either forced to a fixed PDG id or drawn from
    # the renormalized fragmentation cumulative.
    if force_species is not None:
        species_pdg = np.full(n_events, int(force_species), dtype=int)
        m = np.full(n_events, MESON_MASSES[int(force_species)], dtype=float)
    else:
        species_list = QUARK_MESON_MAP[quark]
        pdg_ids = np.array([s[0] for s in species_list])
        fractions = np.array([s[1] for s in species_list])
        frac_cumsum = np.cumsum(fractions)
        frac_cumsum /= frac_cumsum[-1]  # normalize

        u_species = rng.random(n_events)
        species_idx = np.searchsorted(frac_cumsum, u_species)
        species_idx = np.clip(species_idx, 0, len(pdg_ids) - 1)
        species_pdg = pdg_ids[species_idx]

        # Meson mass for each event
        m = np.array([MESON_MASSES[pdg] for pdg in species_pdg])

    # 4-vector from (pT, y, φ, m):
    #   mT = sqrt(pT² + m²)
    #   pz = mT * sinh(y)
    #   E  = mT * cosh(y)
    #   px = pT * cos(φ)
    #   py = pT * sin(φ)
    mt = np.sqrt(pt_sampled**2 + m**2)
    pz = mt * np.sinh(y_sampled)
    E = mt * np.cosh(y_sampled)
    px = pt_sampled * np.cos(phi)
    py = pt_sampled * np.sin(phi)

    return {
        'E': E, 'px': px, 'py': py, 'pz': pz,
        'species_pdg': species_pdg,
        'pt': pt_sampled, 'y': y_sampled, 'phi': phi,
    }
