from __future__ import annotations

import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import numpy as np
import trimesh

from decay.rhn_decay_library import load_decay_events, pick_decay_events, select_decay_file
from geometry.per_parent_efficiency import eta_phi_to_directions

try:
    from particle import Particle
except Exception:
    Particle = None


@dataclass
class DecaySelection:
    separation_m: float
    seed: int = 12345
    p_min_GeV: float = 0.5
    max_separation_m: float | None = None
    separation_policy: str = "all-pairs-min"


@dataclass
class DecayCache:
    charged_directions: List[List[np.ndarray]]
    decay_u: np.ndarray
    hit_indices: np.ndarray


def _normalize_separation_policy(policy: str) -> str:
    policy_norm = str(policy).strip().lower().replace("_", "-")
    if policy_norm not in {"all-pairs-min", "any-pair-window"}:
        raise ValueError(
            f"Unsupported separation policy '{policy}'. "
            "Use 'all-pairs-min' or 'any-pair-window'."
        )
    return policy_norm


def pairwise_separation_pass(
    points: np.ndarray,
    min_separation_m: float,
    max_separation_m: float | None = None,
    separation_policy: str = "all-pairs-min",
) -> bool:
    if len(points) < 2:
        return False

    policy_norm = _normalize_separation_policy(separation_policy)
    min_sep = float(min_separation_m)
    max_sep = None if max_separation_m is None else float(max_separation_m)

    diffs = points[:, None, :] - points[None, :, :]
    dists = np.linalg.norm(diffs, axis=2)
    tri = np.triu_indices(len(points), k=1)
    pair_dists = dists[tri]
    if pair_dists.size == 0:
        return False

    if policy_norm == "all-pairs-min":
        # Explicitly enforce: min(pairwise) >= min_sep and, if requested, max(pairwise) <= max_sep.
        min_pair = float(np.min(pair_dists))
        max_pair = float(np.max(pair_dists))
        if min_pair < min_sep:
            return False
        if max_sep is not None and max_pair > max_sep:
            return False
        return True

    upper = np.inf if max_sep is None else max_sep
    return bool(np.any((pair_dists >= min_sep) & (pair_dists <= upper)))


def _unit_vector(vec: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(vec)
    if norm <= 0.0:
        return np.array([0.0, 0.0, 1.0], dtype=float)
    return vec / norm


def _rotation_matrix_from_z(target: np.ndarray) -> np.ndarray:
    target = _unit_vector(target)
    z_axis = np.array([0.0, 0.0, 1.0], dtype=float)
    dot = np.clip(np.dot(z_axis, target), -1.0, 1.0)
    if np.isclose(dot, 1.0):
        return np.eye(3)
    if np.isclose(dot, -1.0):
        return np.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]])
    axis = np.cross(z_axis, target)
    axis = _unit_vector(axis)
    angle = np.arccos(dot)
    K = np.array(
        [
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ],
        dtype=float,
    )
    return np.eye(3) + np.sin(angle) * K + (1.0 - np.cos(angle)) * (K @ K)


def _boost_along_direction(
    E: float, p_vec: np.ndarray, beta: float, direction: np.ndarray
) -> Tuple[float, np.ndarray]:
    direction = _unit_vector(direction)
    if beta <= 0.0:
        return E, p_vec
    gamma = 1.0 / np.sqrt(1.0 - beta * beta)
    p_parallel = np.dot(p_vec, direction)
    p_perp_vec = p_vec - p_parallel * direction
    E_prime = gamma * (E + beta * p_parallel)
    p_parallel_prime = gamma * (p_parallel + beta * E)
    p_vec_prime = p_perp_vec + p_parallel_prime * direction
    return E_prime, p_vec_prime


def _charge_from_pdg(pid: int) -> float:
    if Particle is None:
        return 0.0
    try:
        particle = Particle.from_pdgid(int(pid))
    except Exception:
        return 0.0
    if particle is None or particle.charge is None:
        return 0.0
    return float(particle.charge)


def _is_charged(pid: int) -> bool:
    charge = _charge_from_pdg(pid)
    if charge != 0.0:
        return True
    return abs(pid) in {11, 13, 15, 211, 321, 2212, 24}


def _first_intersection_point(
    mesh: trimesh.Trimesh, origin: np.ndarray, direction: np.ndarray
) -> np.ndarray | None:
    direction = _unit_vector(direction)
    try:
        locations, index_ray, _ = mesh.ray.intersects_location(
            ray_origins=np.array([origin], dtype=float),
            ray_directions=np.array([direction], dtype=float),
        )
    except Exception:
        return None
    if len(locations) == 0:
        return None
    deltas = locations - origin[None, :]
    t = np.einsum("ij,ij->i", deltas, direction[None, :])
    mask = np.isfinite(t) & (t > 0.0)
    if not np.any(mask):
        return None
    t = t[mask]
    locations = locations[mask]
    return locations[int(np.argmin(t))]


def _batch_first_intersections(
    mesh: trimesh.Trimesh, origins: np.ndarray, directions: np.ndarray
) -> List[np.ndarray | None]:
    n = len(directions)
    if n == 0:
        return []

    norms = np.linalg.norm(directions, axis=1, keepdims=True)
    norms = np.where(norms > 0, norms, 1.0)
    directions = directions / norms

    ray_origins = np.tile(origins.reshape(1, 3), (n, 1))

    try:
        locations, index_ray, _ = mesh.ray.intersects_location(
            ray_origins=ray_origins,
            ray_directions=directions,
        )
    except Exception:
        return [None] * n

    if len(locations) == 0:
        return [None] * n

    results: List[np.ndarray | None] = [None] * n
    for ray_idx in range(n):
        mask = index_ray == ray_idx
        if not np.any(mask):
            continue
        pts = locations[mask]
        deltas = pts - origins
        t = np.einsum("ij,j->i", deltas, directions[ray_idx])
        valid = np.isfinite(t) & (t > 0.0)
        if not np.any(valid):
            continue
        t_valid = t[valid]
        pts_valid = pts[valid]
        results[ray_idx] = pts_valid[np.argmin(t_valid)]

    return results


def build_decay_cache(
    geom_df,
    mass_GeV: float,
    flavour: str,
    selection: DecaySelection,
    verbose: bool = True,
) -> DecayCache:
    rng = np.random.default_rng(selection.seed)
    decay_file = select_decay_file(flavour, mass_GeV)
    decay_events = load_decay_events(decay_file.path)

    hits_tube = geom_df["hits_tube"].to_numpy(dtype=bool)
    hit_indices = np.flatnonzero(hits_tube)
    n_hits = len(hit_indices)
    n_total = len(geom_df)

    if verbose:
        print(f"build_decay_cache: {n_hits} hits / {n_total} total", flush=True)

    beta_gamma = geom_df["beta_gamma"].to_numpy(dtype=float)
    eta = geom_df["eta"].to_numpy(dtype=float)
    phi = geom_df["phi"].to_numpy(dtype=float)
    directions = eta_phi_to_directions(eta, phi)

    selected_events = pick_decay_events(rng, decay_events, n_hits)
    decay_u = rng.uniform(0.0, 1.0, size=n_total)

    charged_directions: List[List[np.ndarray]] = []
    for i, idx in enumerate(hit_indices):
        if verbose and i > 0 and i % 500 == 0:
            print(f"  build_decay_cache: {i}/{n_hits}", flush=True)

        direction = directions[idx]
        if not np.all(np.isfinite(direction)):
            charged_directions.append([])
            continue

        bg = beta_gamma[idx]
        beta = bg / np.sqrt(1.0 + bg * bg)
        rot = _rotation_matrix_from_z(direction)

        dirs: List[np.ndarray] = []
        for E_rf, px_rf, py_rf, pz_rf, _, pid in selected_events[i]:
            if not _is_charged(pid):
                continue
            p_rf = np.array([px_rf, py_rf, pz_rf], dtype=float)
            p_rot = rot @ p_rf
            _, p_lab = _boost_along_direction(E_rf, p_rot, beta, direction)
            if not np.all(np.isfinite(p_lab)):
                continue
            if np.linalg.norm(p_lab) < selection.p_min_GeV:
                continue
            dirs.append(_unit_vector(p_lab))
        charged_directions.append(dirs)

    assert len(charged_directions) == n_hits, f"Mismatch: {len(charged_directions)} != {n_hits}"
    if verbose:
        non_empty = sum(1 for d in charged_directions if len(d) > 0)
        print(f"  build_decay_cache done: {non_empty}/{n_hits} have charged tracks", flush=True)

    return DecayCache(
        charged_directions=charged_directions,
        decay_u=decay_u,
        hit_indices=hit_indices,
    )


def save_decay_cache(cache: DecayCache, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump({
            "charged_directions": cache.charged_directions,
            "decay_u": cache.decay_u,
            "hit_indices": cache.hit_indices,
        }, f)


def load_decay_cache(path: Path) -> DecayCache:
    with open(path, "rb") as f:
        data = pickle.load(f)
    return DecayCache(
        charged_directions=data["charged_directions"],
        decay_u=data["decay_u"],
        hit_indices=data["hit_indices"],
    )


def compute_decay_acceptance(
    geom_df,
    mass_GeV: float,
    flavour: str,
    ctau0_m: float,
    mesh: trimesh.Trimesh,
    selection: DecaySelection,
    decay_cache: DecayCache | None = None,
    verbose: bool = False,
) -> np.ndarray:
    n_total = len(geom_df)
    separation_pass = np.zeros(n_total, dtype=bool)

    if decay_cache is None:
        decay_cache = build_decay_cache(geom_df, mass_GeV, flavour, selection, verbose=verbose)

    hit_indices = decay_cache.hit_indices
    charged_dirs_list = decay_cache.charged_directions
    decay_u = decay_cache.decay_u

    entry = geom_df["entry_distance"].to_numpy(dtype=float)
    path_length = geom_df["path_length"].to_numpy(dtype=float)
    beta_gamma = geom_df["beta_gamma"].to_numpy(dtype=float)
    eta = geom_df["eta"].to_numpy(dtype=float)
    phi = geom_df["phi"].to_numpy(dtype=float)
    directions = eta_phi_to_directions(eta, phi)

    for i, idx in enumerate(hit_indices):
        if not np.isfinite(entry[idx]) or not np.isfinite(path_length[idx]):
            continue
        if path_length[idx] <= 0.0:
            continue

        lam = beta_gamma[idx] * ctau0_m
        if lam <= 0.0:
            continue

        direction = directions[idx]
        if not np.all(np.isfinite(direction)):
            continue

        exp_entry = np.exp(-entry[idx] / lam)
        exp_exit = np.exp(-(entry[idx] + path_length[idx]) / lam)
        denom = exp_entry - exp_exit
        if denom <= 0.0:
            continue

        u = decay_u[idx]
        value = exp_entry - u * denom
        if value <= 0.0:
            continue
        decay_distance = -lam * np.log(value)
        decay_pos = decay_distance * direction

        dirs = charged_dirs_list[i]
        if len(dirs) < 2:
            continue

        dirs_arr = np.array(dirs)
        hits = _batch_first_intersections(mesh, decay_pos, dirs_arr)
        charged_hits = [h for h in hits if h is not None]

        if len(charged_hits) < 2:
            continue

        points = np.stack(charged_hits, axis=0)
        if pairwise_separation_pass(
            points=points,
            min_separation_m=selection.separation_m,
            max_separation_m=selection.max_separation_m,
            separation_policy=selection.separation_policy,
        ):
            separation_pass[idx] = True

    return separation_pass


def compute_separation_pass_static(
    geom_df,
    decay_cache: DecayCache,
    mesh: trimesh.Trimesh,
    separation_m: float,
    max_separation_m: float | None = None,
    separation_policy: str = "all-pairs-min",
) -> np.ndarray:
    n_total = len(geom_df)
    separation_pass = np.zeros(n_total, dtype=bool)

    hit_indices = decay_cache.hit_indices
    charged_dirs_list = decay_cache.charged_directions

    entry = geom_df["entry_distance"].to_numpy(dtype=float)
    path_length = geom_df["path_length"].to_numpy(dtype=float)
    eta = geom_df["eta"].to_numpy(dtype=float)
    phi = geom_df["phi"].to_numpy(dtype=float)
    directions = eta_phi_to_directions(eta, phi)

    for i, idx in enumerate(hit_indices):
        if not np.isfinite(entry[idx]) or not np.isfinite(path_length[idx]):
            continue
        if path_length[idx] <= 0.0:
            continue

        direction = directions[idx]
        if not np.all(np.isfinite(direction)):
            continue

        decay_distance = entry[idx] + 0.5 * path_length[idx]
        decay_pos = decay_distance * direction

        dirs = charged_dirs_list[i]
        if len(dirs) < 2:
            continue

        dirs_arr = np.array(dirs)
        hits = _batch_first_intersections(mesh, decay_pos, dirs_arr)
        charged_hits = [h for h in hits if h is not None]

        if len(charged_hits) < 2:
            continue

        points = np.stack(charged_hits, axis=0)
        if pairwise_separation_pass(
            points=points,
            min_separation_m=separation_m,
            max_separation_m=max_separation_m,
            separation_policy=separation_policy,
        ):
            separation_pass[idx] = True

    return separation_pass
