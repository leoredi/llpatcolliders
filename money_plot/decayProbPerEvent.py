import numpy as np
import pandas as pd
import trimesh
from tqdm import tqdm
import argparse
import os
import re

# ------------------ Constants ------------------

SPEED_OF_LIGHT = 299792458.0  # m/s

# ------------------ Kinematics helpers ------------------

def eta_phi_to_direction(eta, phi):
    """
    Convert pseudorapidity (eta) and azimuthal angle (phi) to a 3D unit vector.
    """
    theta = 2 * np.arctan(np.exp(-eta))
    dx = np.sin(theta) * np.cos(phi)
    dy = np.sin(theta) * np.sin(phi)
    dz = np.cos(theta)
    v = np.array([dx, dy, dz])
    return v / np.linalg.norm(v)


def calculate_decay_length(momentum, mass, lifetime_seconds):
    """
    Lab-frame decay length in meters: L = gamma * beta * c * tau
    momentum, mass in GeV; lifetime in seconds.
    """
    energy = np.sqrt(momentum**2 + mass**2)
    beta = momentum / energy
    gamma = energy / mass
    return beta * gamma * SPEED_OF_LIGHT * lifetime_seconds

# ------------------ Geometry: tube mesh ------------------

def create_tube_mesh(path_points, radius=1.0, n_segments=16):
    """Create a tube mesh along a path with circular cross-section."""
    vertices = []
    faces = []
    
    for i in range(len(path_points)):
        if i == 0:
            tangent = path_points[1] - path_points[0]
        elif i == len(path_points) - 1:
            tangent = path_points[i] - path_points[i-1]
        else:
            tangent = path_points[i+1] - path_points[i-1]
        
        tangent = tangent / np.linalg.norm(tangent)
        
        if abs(tangent[2]) < 0.9:
            up = np.array([0, 0, 1])
        else:
            up = np.array([1, 0, 0])
        
        right = np.cross(tangent, up)
        right = right / np.linalg.norm(right)
        up = np.cross(right, tangent)
        up = up / np.linalg.norm(up)
        
        for j in range(n_segments):
            angle = 2 * np.pi * j / n_segments
            offset = radius * (np.cos(angle) * right + np.sin(angle) * up)
            vertex = path_points[i] + offset
            vertices.append(vertex)
        
        if i > 0:
            for j in range(n_segments):
                v1 = (i-1) * n_segments + j
                v2 = (i-1) * n_segments + (j + 1) % n_segments
                v3 = i * n_segments + (j + 1) % n_segments
                v4 = i * n_segments + j
                
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
    last_ring_start = (len(path_points) - 1) * n_segments
    for j in range(n_segments):
        v1 = last_ring_start + j
        v2 = last_ring_start + (j + 1) % n_segments
        faces.append([center_end, v2, v1])
    
    return np.array(vertices), np.array(faces)


def build_tube_mesh():
    # === your correctedVert, as in the Higgs script ===
    correctedVert = [
        (-86.57954338701529, 0.1882163986665546  ),
        (-1731.590867740335, 3.764327973349282   ),
        (-3549.761278867689, 7.716872345365118   ),
        (-5887.408950317142, 12.798715109387558  ),
        (-8053.403266181902, -504.23173203003535 ),
        (-10046.991360867298, -1282.5065405198511),
        (-11783.350377373874, -2930.9057600491833),
        (-12913.652590171332, -4580.622494369192 ),
        (-13095.344153684957, -7536.749251839814 ),
        (-13099.610392054752, -9015.000846973791 ),
        (-13278.792403586143, -11101.567842600896),
        (-13372.39869252341, -13536.146959364076 ),
        (-13292.093029091975, -15710.234580371536),
        (-12779.140603923677, -17972.21925955668   ), 
        (-11659.12755425337, -19887.69754879509    ),
        (-10105.714877251532, -21630.204967658145  ),
        (-7512.845769209047, -23201.0590309365     ),
        (-5262.530506741277, -23466.820585854904   ),
        (-2751.72374851779, -23472.278861416264    ),
        (-241.41890069074725, -23651.64908934632   ),
        (1749.6596420124115, -23742.93404270002    ),
        (3827.568683300815, -23747.45123626804     ),
        (6078.6368113632525, -23752.344862633392   ),
        (8502.613071001502, -23844.570897980426    ),
        (11446.568501358292, -23764.01427935077    ),
        (13438.399909656131, -23594.431304151418   ),
        (15777.051401898476, -23251.689242178036   ),
        (18289.614846509525, -22648.455684448927   ),
        (20889.761655300477, -21697.58643838109    ),
        (23143.841245741598, -20659.00835053422    ),
        (25486.006110759066, -19098.88262197991    ),
        (27742.09334278597, -17364.656724658227    ),
        (28871.391734790544, -16062.763895075637   ),
        (30781.662703665817, -14153.873179790575   ),
        (32518.021720172394, -12505.473960261239   ),
        (34513.49197884447, -11075.029330388788    ),
        (36636.57295581305, -10427.47081077351     ),
        (38759.40297758341, -9866.868267342572     ),
        (41357.416667189485, -9655.12481884172     ),
        (43694.93886103982, -9703.684649697909     ),
        (46379.03018363646, -9666.041369964427     ),
        (49409.43967978114, -9629.150955825604     ),
        (51660.88424064092, -9503.610617914434     ),
        (54258.0195870532, -9596.213086058811      ),
        (57028.564975437745, -9602.236010816167    ),
        (59539.87364405768, -9433.782334008818     ),
        (62050.42944708294, -9526.196585754526     )
    ]

    correctedVertWithShift = []
    for x, y in correctedVert:
        correctedVertWithShift.append(((x - 11908.8279764855)/1000.0,
                                       (y + 13591.106147774964)/1000.0))

    Z_POSITION = 22.0
    path_3d = np.array([[x, y, Z_POSITION] for x, y in correctedVertWithShift])

    tube_radius = 1.4 * 1.1  # same as before
    vertices, faces = create_tube_mesh(path_3d, radius=tube_radius, n_segments=32)
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
    if mesh.volume < 0:
        mesh = mesh.copy()
        mesh.invert()
    return mesh

# ------------------ Geometry caching ------------------

def precompute_geometry(csv_file, mesh, origin=(0.0, 0.0, 0.0), out_file=None):
    """
    Compute entry/exit distances of each HNL through the tube.
    This is independent of lifetime, so we can reuse it for all τ.
    """
    origin = np.array(origin, dtype=float)

    df = pd.read_csv(csv_file)
    df.columns = [c.strip() for c in df.columns]

    required = ['event', 'eta', 'phi', 'momentum', 'mass']
    for col in required:
        if col not in df.columns:
            raise ValueError(f"CSV {csv_file} missing required column '{col}'")

    df['hits_tube'] = False
    df['entry_distance'] = np.nan
    df['exit_distance'] = np.nan
    df['path_length_in_tube'] = np.nan

    print(f"Precomputing geometry for {len(df)} particles "
          f"from {df['event'].nunique()} events...")

    # Use ray intersections
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Geometry (rays)"):
        eta = row['eta']
        phi = row['phi']
        direction = eta_phi_to_direction(eta, phi)

        locations, _, _ = mesh.ray.intersects_location(
            ray_origins=[origin],
            ray_directions=[direction]
        )

        if len(locations) >= 2:
            df.at[idx, 'hits_tube'] = True
            distances = [np.linalg.norm(loc - origin) for loc in locations]
            distances = sorted(distances)
            entry = distances[0]
            exit_ = distances[1]
            df.at[idx, 'entry_distance'] = entry
            df.at[idx, 'exit_distance'] = exit_
            df.at[idx, 'path_length_in_tube'] = exit_ - entry

    if out_file is None:
        base, ext = os.path.splitext(csv_file)
        out_file = base + "_geomcached.csv"

    df.to_csv(out_file, index=False)
    print(f"Saved geometry-cached CSV: {out_file}")
    return out_file

# ------------------ Lifetime scan & BR limits ------------------

def scan_lifetimes_from_cached(geom_csv, sigma_gen_pb, lumi_fb=3000.0,
                               lifetimes_s=None, zero_background=True):
    """
    Given a geometry-cached CSV and Pythia σ, scan lifetimes and compute BR limits.
    """
    df = pd.read_csv(geom_csv)
    df.columns = [c.strip() for c in df.columns]

    required = ['event', 'momentum', 'mass', 'hits_tube',
                'entry_distance', 'path_length_in_tube']
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Geometry CSV {geom_csv} missing '{col}'")

    if lifetimes_s is None:
        # e.g. from 10^-10 s to 10^-5 s
        lifetimes_s = np.logspace(-10, -5, 25)

    n_events = df['event'].nunique()
    sigma_fb = sigma_gen_pb * 1e3  # pb -> fb

    print(f"\nTotal events in file: {n_events}")
    print(f"Pythia σ_gen = {sigma_gen_pb:.3e} pb = {sigma_fb:.3e} fb")
    print(f"Luminosity    = {lumi_fb:.1f} fb^-1\n")

    results = {
        'lifetime_s': [],
        'ctau_m': [],
        'mean_event_decay_prob': [],
        'BR_limit': [],
    }

    # how many signal events we require (3 for 95% CL, zero background)
    N_req = 3.0 if zero_background else 5.0  # adjust if needed

    for tau in tqdm(lifetimes_s, desc="Scanning lifetimes"):
        # particle-level: decay lengths
        decay_lengths = calculate_decay_length(df['momentum'].values,
                                               df['mass'].values,
                                               tau)
        df['decay_length'] = decay_lengths

        # per-particle decay probability in tube (only if it hits)
        p_decay = np.zeros(len(df))
        mask_hit = df['hits_tube'].astype(bool) & df['path_length_in_tube'].notna()
        entry = df.loc[mask_hit, 'entry_distance'].values
        path = df.loc[mask_hit, 'path_length_in_tube'].values
        L = decay_lengths[mask_hit]

        p_survive_to_entry = np.exp(- entry / L)
        p_in_tube = p_survive_to_entry * (1.0 - np.exp(- path / L))
        p_decay[mask_hit] = p_in_tube

        df['decay_probability'] = p_decay

        # event-level probability: 1 - prod(1 - p_i)
        event_probs = []
        for evt, grp in df.groupby('event'):
            ps = grp['decay_probability'].values
            p_all_survive = np.prod(1.0 - ps)
            p_evt = 1.0 - p_all_survive
            event_probs.append(p_evt)
        event_probs = np.array(event_probs)

        mean_event_prob = event_probs.mean()

        # For BR(W/B→ℓN) = 1 in the MC, the expected # of signal events is:
        #   N_sig(BR=1) = σ * L * ε   with ε = mean_event_prob
        # For generic BR, N_sig = BR * σ * L * ε
        # => BR_limit = N_req / (σ * L * ε)
        if mean_event_prob > 0:
            br_lim = N_req / (sigma_fb * lumi_fb * mean_event_prob)
        else:
            br_lim = np.inf

        results['lifetime_s'].append(tau)
        results['ctau_m'].append(tau * SPEED_OF_LIGHT)
        results['mean_event_decay_prob'].append(mean_event_prob)
        results['BR_limit'].append(br_lim)

    return results

# ------------------ Main CLI ------------------

def read_sigma_from_meta(meta_file):
    with open(meta_file) as f:
        for line in f:
            if line.strip().startswith("sigma_gen_pb"):
                return float(line.split()[1])
    raise RuntimeError(f"Could not find sigma_gen_pb in {meta_file}")


def main():
    parser = argparse.ArgumentParser(description="HNL-from-B: BR vs cτ scan")
    parser.add_argument("csv_file", help="HNL CSV from Pythia")
    parser.add_argument("--meta", required=True, help=".meta file with sigma_gen_pb")
    parser.add_argument("--lumi", type=float, default=3000.0,
                        help="Integrated luminosity in fb^-1 (default 3000)")
    args = parser.parse_args()

    csv_file = args.csv_file
    meta_file = args.meta
    lumi_fb = args.lumi

    base = os.path.splitext(os.path.basename(csv_file))[0]

    # 1) Build tube mesh
    mesh = build_tube_mesh()
    print("Tube mesh created.")
    print(f"Mesh bounds: {mesh.bounds}")

    # 2) Geometry caching - check if already exists
    base_path, ext = os.path.splitext(csv_file)
    geom_csv = base_path + "_geomcached.csv"

    if os.path.exists(geom_csv):
        print(f"Using existing geometry cache: {geom_csv}")
    else:
        geom_csv = precompute_geometry(csv_file, mesh, origin=(0.0, 0.0, 0.0))

    # 3) Read σ from meta
    sigma_gen_pb = read_sigma_from_meta(meta_file)

    # 4) Lifetime scan
    lifetimes_s = np.logspace(-10, -5, 25)
    results = scan_lifetimes_from_cached(geom_csv,
                                         sigma_gen_pb=sigma_gen_pb,
                                         lumi_fb=lumi_fb,
                                         lifetimes_s=lifetimes_s)

    # 5) Save results for coupling-mapping step
    out_csv = f"../output/csv/analysis/{base}_BR_vs_ctau.csv"
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    df_out = pd.DataFrame(results)
    df_out.to_csv(out_csv, index=False)
    print(f"\nSaved BR vs cτ data to: {out_csv}")


if __name__ == "__main__":
    main()