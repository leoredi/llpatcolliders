import sys

import numpy as np
import pandas as pd
import trimesh
from tqdm import tqdm

# Try to import the specific RTreeError used by trimesh's ray code.
# If it's not available for some reason, fall back to a generic Exception subclass.
try:
    from rtree.exceptions import RTreeError
except Exception:  # pragma: no cover
    class RTreeError(Exception):
        pass


def eta_phi_to_direction(eta: float, phi: float) -> np.ndarray:
    """
    Convert pseudorapidity (eta) and azimuthal angle (phi) to a 3D unit direction vector.

    Args:
        eta: Pseudorapidity
        phi: Azimuthal angle in radians

    Returns:
        Normalized 3D direction vector [dx, dy, dz]
    """
    # Convert eta to polar angle theta
    theta = 2.0 * np.arctan(np.exp(-eta))

    # Convert to Cartesian direction
    dx = np.sin(theta) * np.cos(phi)
    dy = np.sin(theta) * np.sin(phi)
    dz = np.cos(theta)

    return np.array([dx, dy, dz], dtype=float)


def eta_phi_to_directions(eta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """
    Vectorized version of eta_phi_to_direction.

    Args:
        eta: array of pseudorapidity values
        phi: array of azimuthal angles (radians)

    Returns:
        Array of shape (N, 3) with unit direction vectors.
    """
    eta = np.asarray(eta, dtype=float)
    phi = np.asarray(phi, dtype=float)

    theta = 2.0 * np.arctan(np.exp(-eta))
    dx = np.sin(theta) * np.cos(phi)
    dy = np.sin(theta) * np.sin(phi)
    dz = np.cos(theta)

    directions = np.stack([dx, dy, dz], axis=1)
    norms = np.linalg.norm(directions, axis=1)
    norms = np.where(norms > 0.0, norms, np.nan)
    return directions / norms[:, None]


def create_tube_mesh(path_points: np.ndarray,
                     radius: float = 1.0,
                     n_segments: int = 16) -> tuple[np.ndarray, np.ndarray]:
    """
    Create a tube mesh along a 3D polyline path with a circular cross-section.

    Args:
        path_points: array-like of shape (N, 3), the 3D path of the tube axis
        radius: tube radius (same units as path_points)
        n_segments: number of points around the circular cross-section

    Returns:
        vertices: (Nv, 3) array of vertex positions
        faces:    (Nf, 3) array of triangular faces (indices into vertices)
    """
    path_points = np.asarray(path_points, dtype=float)
    vertices: list[np.ndarray] = []
    faces: list[list[int]] = []

    n_points = len(path_points)
    if n_points < 2:
        raise ValueError("Need at least two path points to build a tube.")

    for i in range(n_points):
        if i == 0:
            tangent = path_points[1] - path_points[0]
        elif i == n_points - 1:
            tangent = path_points[i] - path_points[i - 1]
        else:
            tangent = path_points[i + 1] - path_points[i - 1]

        tangent = tangent / np.linalg.norm(tangent)

        # Choose an "up" vector not parallel to tangent
        if abs(tangent[2]) < 0.9:
            up = np.array([0.0, 0.0, 1.0])
        else:
            up = np.array([1.0, 0.0, 0.0])

        right = np.cross(tangent, up)
        right = right / np.linalg.norm(right)
        up = np.cross(right, tangent)
        up = up / np.linalg.norm(up)

        # Build the ring of vertices at this path point
        for j in range(n_segments):
            angle = 2.0 * np.pi * j / n_segments
            offset = radius * (np.cos(angle) * right + np.sin(angle) * up)
            vertex = path_points[i] + offset
            vertices.append(vertex)

        # Connect this ring to the previous ring
        if i > 0:
            for j in range(n_segments):
                v1 = (i - 1) * n_segments + j
                v2 = (i - 1) * n_segments + (j + 1) % n_segments
                v3 = i * n_segments + (j + 1) % n_segments
                v4 = i * n_segments + j

                # Two triangles per quad
                faces.append([v1, v4, v3])
                faces.append([v1, v3, v2])

    # Cap the ends
    center_start = len(vertices)
    vertices.append(path_points[0])
    for j in range(n_segments):
        v1 = j
        v2 = (j + 1) % n_segments
        faces.append([center_start, v1, v2])

    center_end = len(vertices)
    vertices.append(path_points[-1])
    last_ring_start = (n_points - 1) * n_segments
    for j in range(n_segments):
        v1 = last_ring_start + j
        v2 = last_ring_start + (j + 1) % n_segments
        faces.append([center_end, v2, v1])

    return np.array(vertices, dtype=float), np.array(faces, dtype=int)


def build_drainage_gallery_mesh() -> trimesh.Trimesh:
    """
    Build the CMS drainage-gallery tube mesh at z = 22 m.

    Uses the 'correctedVert' polyline and shifts/scales exactly as in the original
    Higgs→LLP script:
      - subtract 11908.8279764855 in x, add 13591.106147774964 in y
      - divide by 1000 to go from mm (?) to m
      - set z = 22 m
      - radius = 1.4 * 1.1
      - n_segments = 32

    Returns:
        mesh: trimesh.Trimesh object representing the tube volume.
    """
    correctedVert = [
        (-86.57954338701529, 0.1882163986665546),
        (-1731.590867740335, 3.764327973349282),
        (-3549.761278867689, 7.716872345365118),
        (-5887.408950317142, 12.798715109387558),
        (-8053.403266181902, -504.23173203003535),
        (-10046.991360867298, -1282.5065405198511),
        (-11783.350377373874, -2930.9057600491833),
        (-12913.652590171332, -4580.622494369192),
        (-13095.344153684957, -7536.749251839814),
        (-13099.610392054752, -9015.000846973791),
        (-13278.792403586143, -11101.567842600896),
        (-13372.39869252341, -13536.146959364076),
        (-13292.093029091975, -15710.234580371536),
        (-12779.140603923677, -17972.21925955668),
        (-11659.12755425337, -19887.69754879509),
        (-10105.714877251532, -21630.204967658145),
        (-7512.845769209047, -23201.0590309365),
        (-5262.530506741277, -23466.820585854904),
        (-2751.72374851779, -23472.278861416264),
        (-241.41890069074725, -23651.64908934632),
        (1749.6596420124115, -23742.93404270002),
        (3827.568683300815, -23747.45123626804),
        (6078.6368113632525, -23752.344862633392),
        (8502.613071001502, -23844.570897980426),
        (11446.568501358292, -23764.01427935077),
        (13438.399909656131, -23594.431304151418),
        (15777.051401898476, -23251.689242178036),
        (18289.614846509525, -22648.455684448927),
        (20889.761655300477, -21697.58643838109),
        (23143.841245741598, -20659.00835053422),
        (25486.006110759066, -19098.88262197991),
        (27742.09334278597, -17364.656724658227),
        (28871.391734790544, -16062.763895075637),
        (30781.662703665817, -14153.873179790575),
        (32518.021720172394, -12505.473960261239),
        (34513.49197884447, -11075.029330388788),
        (36636.57295581305, -10427.47081077351),
        (38759.40297758341, -9866.868267342572),
        (41357.416667189485, -9655.12481884172),
        (43694.93886103982, -9703.684649697909),
        (46379.03018363646, -9666.041369964427),
        (49409.43967978114, -9629.150955825604),
        (51660.88424064092, -9503.610617914434),
        (54258.0195870532, -9596.213086058811),
        (57028.564975437745, -9602.236010816167),
        (59539.87364405768, -9433.782334008818),
        (62050.42944708294, -9526.196585754526),
    ]

    correctedVertWithShift: list[tuple[float, float]] = []
    for x, y in correctedVert:
        correctedVertWithShift.append(
            (
                (x - 11908.8279764855) / 1000.0,
                (y + 13591.106147774964) / 1000.0,
            )
        )

    Z_POSITION = 22.0  # m
    path_3d = np.array(
        [[x, y, Z_POSITION] for x, y in correctedVertWithShift],
        dtype=float,
    )

    # Physical tube radius is ~1.4 m. Use a small envelope margin to avoid
    # edge-effects from polygonization / discretization in this fast-acceptance
    # model and to be robust against small geometry/alignment uncertainties.
    tube_radius_m = 1.4
    tube_envelope_margin = 1.1
    tube_radius = tube_radius_m * tube_envelope_margin
    vertices, faces = create_tube_mesh(
        path_3d, radius=tube_radius, n_segments=32
    )
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces)

    if mesh.volume < 0:
        mesh = mesh.copy()  # avoid modifying original in place
        mesh.invert()

    print("Tube mesh created at z=22 m")
    print(
        "Mesh bounds: "
        f"X:[{mesh.bounds[0][0]:.1f}, {mesh.bounds[1][0]:.1f}], "
        f"Y:[{mesh.bounds[0][1]:.1f}, {mesh.bounds[1][1]:.1f}], "
        f"Z:[{mesh.bounds[0][2]:.1f}, {mesh.bounds[1][2]:.1f}]"
    )

    return mesh


def preprocess_hnl_csv(
    csv_file: str,
    mesh: trimesh.Trimesh,
    origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
    show_progress: bool | None = None,
) -> pd.DataFrame:
    """
    Read a Pythia HNL CSV and precompute geometry quantities
    (independent of lifetime and couplings).

    Expected input columns (at minimum):
        - event
        - parent_id
        - eta
        - phi
        - momentum
        - mass
        - weight (optional; if absent, set to 1.0)

    Output DataFrame contains original columns plus:
        - beta_gamma        : p / m (dimensionless)
        - hits_tube         : bool
        - entry_distance    : distance from origin to entry point [m]
        - path_length       : path length inside tube [m]

    The event-level logic (per-event decay probability, etc.) is *not*
    handled here; this is purely per-particle geometry.

    Note: Pythia events can contain multiple HNLs from different parent
    mesons (e.g., event 44 might have D0→N, B0→N, Ds→N). We compute
    geometry for EACH HNL individually (one row per HNL), not per event.
    This enables per-parent counting in the analysis layer.

    IMPORTANT - Weight Column Semantics
    ------------------------------------
    The 'weight' column is interpreted as a RELATIVE MC event weight
    (e.g., for phase-space reweighting), NOT an absolute cross-section.

    Absolute cross-sections come from production_xsecs.get_parent_sigma_pb().

    The analysis computes efficiency as:
        ε_parent = Σ(weight * P_decay) / Σ(weight)

    Then multiplies by external cross-section:
        N_sig = L × σ_parent × BR × ε_parent

    If 'weight' contains absolute cross-sections (e.g., σ_gen from Pythia),
    this will DOUBLE-COUNT the cross-section. Keep weights as relative!

    Typical valid weight values:
        - All 1.0 (unweighted MC)
        - 0.1 - 10.0 (phase-space reweighting factors)
        - pythia.info.weight() for weighted generation

    Invalid weight values:
        - pythia.info.sigmaGen() (absolute cross-section in pb)
        - Any value >> 1000 (likely an absolute cross-section mistake)

    Args:
        show_progress: If None (default), auto-detect based on whether stderr
            is a TTY. Set False to disable progress bars (batch/parallel mode).
    """
    if show_progress is None:
        show_progress = sys.stderr.isatty()

    df = pd.read_csv(csv_file)

    # Handle both old and new CSV formats
    # Old format: parent_id, momentum
    # New format: parent_pdg, p
    if "parent_pdg" in df.columns and "parent_id" not in df.columns:
        df["parent_id"] = df["parent_pdg"].abs()  # Use absolute value
    if "p" in df.columns and "momentum" not in df.columns:
        df["momentum"] = df["p"]

    required_cols = ["event", "parent_id", "eta", "phi", "momentum", "mass"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"CSV {csv_file} is missing required columns: {missing}"
        )

    if "weight" not in df.columns:
        df["weight"] = 1.0

    # Initialize geometry columns
    df["hits_tube"] = False
    df["entry_distance"] = np.nan
    df["path_length"] = np.nan

    # beta * gamma = p / m (in natural units)
    if "beta_gamma" not in df.columns:
        if "momentum" in df.columns and "mass" in df.columns:
            df["beta_gamma"] = df["momentum"] / df["mass"]
        else:
            raise ValueError("Cannot compute beta_gamma: missing 'beta_gamma' column or 'momentum'+'mass' columns")

    origin_arr = np.array(origin, dtype=float)

    print(f"Precomputing geometry for {len(df)} HNLs from {csv_file} ...")

    hits_tube = np.zeros(len(df), dtype=bool)
    entry_distance = np.full(len(df), np.nan, dtype=float)
    path_length = np.full(len(df), np.nan, dtype=float)

    eta = df["eta"].to_numpy(dtype=float)
    phi = df["phi"].to_numpy(dtype=float)

    valid_mask = np.isfinite(eta) & np.isfinite(phi)
    directions = eta_phi_to_directions(eta[valid_mask], phi[valid_mask])
    valid_dir_mask = np.all(np.isfinite(directions), axis=1)

    valid_indices = np.flatnonzero(valid_mask)[valid_dir_mask]
    directions = directions[valid_dir_mask]

    bad_dir_count = int(len(df) - len(valid_indices))

    def _intersects_location_safe(
        ray_origins: np.ndarray,
        ray_directions: np.ndarray,
        offset: int = 0,
    ) -> tuple[np.ndarray, np.ndarray, int]:
        try:
            locations, index_ray, _ = mesh.ray.intersects_location(
                ray_origins=ray_origins,
                ray_directions=ray_directions,
            )
            return locations, index_ray + offset, 0
        except RTreeError:
            n = len(ray_directions)
            if n <= 1:
                return np.empty((0, 3), dtype=float), np.empty((0,), dtype=int), 1
            mid = n // 2
            loc_a, idx_a, err_a = _intersects_location_safe(ray_origins[:mid], ray_directions[:mid], offset=offset)
            loc_b, idx_b, err_b = _intersects_location_safe(
                ray_origins[mid:], ray_directions[mid:], offset=offset + mid
            )
            if len(loc_a) == 0:
                loc = loc_b
            elif len(loc_b) == 0:
                loc = loc_a
            else:
                loc = np.concatenate([loc_a, loc_b], axis=0)
            if len(idx_a) == 0:
                idx = idx_b
            elif len(idx_b) == 0:
                idx = idx_a
            else:
                idx = np.concatenate([idx_a, idx_b], axis=0)
            return loc, idx, err_a + err_b

    rtree_errors = 0
    batch_size = 10_000
    n_batches = (len(valid_indices) + batch_size - 1) // batch_size
    iterator = range(0, len(valid_indices), batch_size)
    if show_progress:
        iterator = tqdm(iterator, total=n_batches, desc="Geometry rays", unit="batch")

    for start in iterator:
        end = min(start + batch_size, len(valid_indices))
        idx_chunk = valid_indices[start:end]
        dir_chunk = directions[start:end]
        origin_chunk = np.repeat(origin_arr[None, :], len(dir_chunk), axis=0)

        locations, index_ray, n_err = _intersects_location_safe(origin_chunk, dir_chunk, offset=0)
        rtree_errors += n_err

        if len(locations) == 0:
            continue

        deltas = locations - origin_arr[None, :]
        t = np.einsum("ij,ij->i", deltas, dir_chunk[index_ray])
        mask_t = np.isfinite(t) & (t > 0.0)
        if not np.any(mask_t):
            continue

        t = t[mask_t]
        index_ray = index_ray[mask_t]

        order = np.lexsort((t, index_ray))
        t = t[order]
        index_ray = index_ray[order]

        rays, counts = np.unique(index_ray, return_counts=True)
        pos = 0
        for ray, count in zip(rays, counts):
            ts = t[pos : pos + count]
            pos += count
            if len(ts) < 2:
                continue
            n_pairs = len(ts) // 2
            entry_val = float(ts[0])
            path_val = float(np.sum(ts[1 : 2 * n_pairs : 2] - ts[0 : 2 * n_pairs : 2]))
            df_idx = int(idx_chunk[int(ray)])
            hits_tube[df_idx] = True
            entry_distance[df_idx] = entry_val
            path_length[df_idx] = path_val

    if bad_dir_count > 0:
        print(f"[WARN] Skipped {bad_dir_count} HNL(s) with non-finite or degenerate directions.")

    if rtree_errors > 0:
        print(f"[WARN] Skipped {rtree_errors} ray batch(es)/ray(s) due to RTreeError in ray-mesh intersection.")

    df["hits_tube"] = hits_tube
    df["entry_distance"] = entry_distance
    df["path_length"] = path_length
    return df
