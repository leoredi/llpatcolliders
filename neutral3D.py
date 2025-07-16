import numpy as np
import trimesh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm
from scipy.integrate import quad

# First, recreate the tube mesh at z=22.5
def create_tube_mesh(path_points, radius=1.0, n_segments=16):
    """Create a tube mesh along a path with circular cross-section"""
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
                
                faces.append([v1, v2, v3])
                faces.append([v1, v3, v4])
    
    # Cap the ends
    center_start = len(vertices)
    vertices.append(path_points[0])
    for j in range(n_segments):
        v1 = j
        v2 = (j + 1) % n_segments
        faces.append([center_start, v2, v1])
    
    center_end = len(vertices)
    vertices.append(path_points[-1])
    last_ring_start = (len(path_points) - 1) * n_segments
    for j in range(n_segments):
        v1 = last_ring_start + j
        v2 = last_ring_start + (j + 1) % n_segments
        faces.append([center_end, v1, v2])
    
    return np.array(vertices), np.array(faces)

# Your path vertices
path_vertices = [
    (-12, 12.1), (-13, 12.1), (-14, 12.1), (-15, 12.1), (-16, 12.1), (-17, 12.1), (-18, 12.1), (-18.0, 12.1),
    (-18.627001066178522, 12.06478837540216), (-19.24611723015536, 11.959596308218213), 
    (-19.849562746948934, 11.785746649726857), (-20.429748939058324, 11.545425660253546),
    (-20.979379628485884, 11.241655515678392), (-21.49154289040891, 10.878256301820967),
    (-21.959797974644665, 10.459797974644665), (-22.378256301820965, 9.991542890408908),
    (-22.74165551567839, 9.479379628485885), (-23.045425660253546, 8.929748939058326),
    (-23.285746649726857, 8.349562746948937), (-23.459596308218213, 7.746117230155361),
    (-23.56478837540216, 7.127001066178525), (-23.6, 6.500000000000001), (-23.6, 5.5),
    (-23.6, 4.5), (-23.6, 3.5), (-23.6, 2.5), (-23.6, 1.5), (-23.6, 0.5), (-23.6, -0.5),
    (-23.6, -1.5), (-23.6, -2.5), (-23.6, -3.5), (-23.6, -3.499999999999999),
    (-23.56478837540216, -4.127001066178523), (-23.459596308218213, -4.74611723015536),
    (-23.285746649726857, -5.349562746948934), (-23.045425660253546, -5.9297489390583245),
    (-22.741655515678392, -6.479379628485884), (-22.37825630182097, -6.9915428904089065),
    (-21.959797974644665, -7.459797974644665), (-21.49154289040891, -7.878256301820966),
    (-20.979379628485887, -8.24165551567839), (-20.429748939058328, -8.545425660253546),
    (-19.849562746948937, -8.785746649726857), (-19.246117230155363, -8.959596308218213),
    (-18.627001066178526, -9.06478837540216), (-18.0, -9.1), (-18, -9.1), (5, -9.1),
    (6.640000000000001, -8.5), (10.505298148418134, -7.7631578947368425),
    (13.66077853914565, -7.026315789473684), (16.17845166933955, -6.2894736842105265),
    (18.130328036156875, -5.552631578947368), (19.58841813675463, -4.815789473684211),
    (20.62473246828984, -4.078947368421053), (21.311281527919522, -3.3421052631578947),
    (21.7200758128007, -2.605263157894737), (21.92312582009039, -1.8684210526315796),
    (21.99244204694562, -1.1315789473684212), (22.0000349905234, -0.3947368421052637),
    (22.017915147980755, 0.3421052631578938), (22.118093016474706, 1.0789473684210513),
    (22.37257909316227, 1.8157894736842106), (22.853383875200468, 2.552631578947368),
    (23.63251785974632, 3.2894736842105257), (24.781991543956842, 4.026315789473683),
    (26.37381542498906, 4.763157894736841), (28.48, 5.5)
]
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
(62050.42944708294, -9526.196585754526     )]

correctedVertWithShift = []
for x,y in correctedVert:
    correctedVertWithShift.append(((x- 11908.8279764855)/1000,(y+13591.106147774964)/1000))
# Create 3D path at z=22
path_vertices = correctedVertWithShift
Z_POSITION = 22
path_3d = np.array([[x, y, Z_POSITION] for x, y in path_vertices])

# Create the tube mesh
tube_radius = 1.4*1.1
vertices, faces = create_tube_mesh(path_3d, radius=tube_radius, n_segments=32)
mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
if mesh.volume < 0:
    mesh.invert()

print("Tube mesh created at z=22")
print(f"Mesh bounds: X:[{mesh.bounds[0][0]:.1f}, {mesh.bounds[1][0]:.1f}], "
      f"Y:[{mesh.bounds[0][1]:.1f}, {mesh.bounds[1][1]:.1f}], "
      f"Z:[{mesh.bounds[0][2]:.1f}, {mesh.bounds[1][2]:.1f}]")

# Particle decay analysis functions
def exponential_decay_probability(distance, decay_length):
    """
    Calculate probability of decay within a given distance
    
    Args:
        distance: Distance traveled
        decay_length: Characteristic decay length (lifetime * speed)
    
    Returns:
        Probability of decay before reaching this distance
    """
    return 1 - np.exp(-distance / decay_length)

def survival_probability(distance, decay_length):
    """
    Calculate probability of surviving to a given distance
    
    Args:
        distance: Distance traveled
        decay_length: Characteristic decay length
    
    Returns:
        Probability of surviving to this distance
    """
    return np.exp(-distance / decay_length)

def calculate_decay_fraction_in_tube(mesh, origin, decay_length, 
                                   azimuth_samples=180, elevation_samples=90,
                                   show_progress=True):
    """
    Calculate fraction of particles that decay within the tube volume
    
    Args:
        mesh: Trimesh object representing the tube
        origin: Starting point of particles [x, y, z]
        decay_length: Characteristic decay length (e.g., lifetime * speed of light)
        azimuth_samples: Number of azimuth angle samples
        elevation_samples: Number of elevation angle samples
        show_progress: Whether to show progress bar
    
    Returns:
        Dictionary with detailed results
    """
    origin = np.array(origin)
    
    # Create angle grids
    azimuths = np.linspace(0, 2*np.pi, azimuth_samples, endpoint=False)
    elevations = np.linspace(-np.pi/2, np.pi/2, elevation_samples)
    
    total_rays = len(azimuths) * len(elevations)
    
    # Results storage
    decay_fractions = []
    azimuth_all = []
    elevation_all = []
    weighted_decay_fraction = 0
    total_solid_angle = 0
    rays_hitting_tube = 0
    decay_in_tube_map = np.zeros((elevation_samples, azimuth_samples))
    
    if show_progress:
        pbar = tqdm(total=total_rays, desc="Analyzing particle decay")
    
    # Analyze each direction
    for i, elevation in enumerate(elevations):
        for j, azimuth in enumerate(azimuths):
            # Convert to direction vector
            dx = np.cos(elevation) * np.cos(azimuth)
            dy = np.cos(elevation) * np.sin(azimuth)
            dz = np.sin(elevation)
            direction = np.array([dx, dy, dz])
            
            # Find intersections with tube
            locations, _, _ = mesh.ray.intersects_location(
                ray_origins=[origin],
                ray_directions=[direction]
            )
            
            if len(locations) > 0:
                rays_hitting_tube += 1
                
                # Sort intersections by distance
                distances = [np.linalg.norm(loc - origin) for loc in locations]
                sorted_indices = np.argsort(distances)
                sorted_distances = [distances[idx] for idx in sorted_indices]
                
                # For a tube, we typically have entry and exit points
                if len(sorted_distances) >= 2:
                    entry_distance = sorted_distances[0]
                    exit_distance = sorted_distances[1]
                    
                    # Probability of surviving to entry point
                    p_survive_to_entry = survival_probability(entry_distance, decay_length)
                    
                    # Probability of decaying between entry and exit
                    # P(decay in tube) = P(survive to entry) * P(decay before exit | survived to entry)
                    p_decay_in_tube = p_survive_to_entry * (1 - survival_probability(exit_distance - entry_distance, decay_length))
                    
                    decay_fractions.append(p_decay_in_tube)
                    azimuth_all.append(azimuth)
                    elevation_all.append(elevation)
                    decay_in_tube_map[i, j] = p_decay_in_tube
                    
                    # Weight by solid angle element
                    solid_angle_element = np.cos(elevation) * (2*np.pi/azimuth_samples) * (np.pi/elevation_samples)
                    weighted_decay_fraction += p_decay_in_tube * solid_angle_element
                    total_solid_angle += solid_angle_element
            
            if show_progress:
                pbar.update(1)
    
    if show_progress:
        pbar.close()
    
    # Calculate overall statistics
    if total_solid_angle > 0:
        average_decay_fraction = weighted_decay_fraction / total_solid_angle
    else:
        average_decay_fraction = 0
    
    # Calculate what fraction of ALL particles (in all directions) decay in tube
    total_decay_in_tube = weighted_decay_fraction / (4 * np.pi)
    
    return {
        'rays_hitting_tube': rays_hitting_tube,
        'total_rays': total_rays,
        'decay_fractions': decay_fractions,
        'average_decay_fraction': average_decay_fraction,
        'total_decay_fraction': total_decay_in_tube,
        'decay_map': decay_in_tube_map,
        'azimuths': azimuths,
        'elevations': elevations,
        'azimuth_all': azimuth_all,
        'elevation_all': elevation_all
    }

# Run the analysis
print("\n" + "="*50)
print("PARTICLE DECAY ANALYSIS")
print("="*50)

# Parameters
origin = [0, 0, 0]
lifetime_distance = 10.0  # meters (or whatever unit your mesh uses)

print(f"\nParameters:")
print(f"Origin: {origin}")
print(f"Particle lifetime (decay length): {lifetime_distance} m")
print(f"Decay probability at distance d: 1 - exp(-d/{lifetime_distance})")

# Run analysis
results = calculate_decay_fraction_in_tube(mesh, origin, lifetime_distance,
                                         azimuth_samples=120, elevation_samples=60)

print(f"\nResults:")
print(f"Rays hitting tube: {results['rays_hitting_tube']} out of {results['total_rays']}")
print(f"Average decay fraction (for particles that hit tube): {results['average_decay_fraction']:.4f} ({results['average_decay_fraction']*100:.2f}%)")
print(f"Total fraction of all particles decaying in tube: {results['total_decay_fraction']:.6f} ({results['total_decay_fraction']*100:.4f}%)")

# Visualization
fig = plt.figure(figsize=(15, 10))

# 1. Decay probability map
ax1 = fig.add_subplot(231)
im = ax1.imshow(results['decay_map'], extent=[0, 360, -90, 90], 
                origin='lower', aspect='auto', cmap='hot')
ax1.set_xlabel('Azimuth (degrees)')
ax1.set_ylabel('Elevation (degrees)')
ax1.set_title('Decay Probability in Tube by Direction')
plt.colorbar(im, ax=ax1, label='P(decay in tube)')

# 2. Histogram of decay probabilities
ax2 = fig.add_subplot(232)
if results['decay_fractions']:
    ax2.hist(results['decay_fractions'], bins=50, edgecolor='black', alpha=0.7)
    ax2.set_xlabel('Decay probability in tube')
    ax2.set_ylabel('Count (directions)')
    ax2.set_title('Distribution of Decay Probabilities')
    ax2.grid(True, alpha=0.3)

# 3. Decay probability vs lifetime
ax3 = fig.add_subplot(233)
lifetimes = np.logspace(0, 3, 50)  # 1 to 1000 meters
total_decay_fractions = []

print("\nCalculating decay fraction vs lifetime...")
for lt in tqdm(lifetimes, desc="Lifetime scan"):
    quick_results = calculate_decay_fraction_in_tube(mesh, origin, lt,
                                                   azimuth_samples=60, elevation_samples=30,
                                                   show_progress=False)
    total_decay_fractions.append(quick_results['total_decay_fraction'])

ax3.loglog(lifetimes, total_decay_fractions, 'b-', linewidth=2)
ax3.axvline(lifetime_distance, color='red', linestyle='--', label=f'Current: {lifetime_distance}m')
ax3.set_xlabel('Decay length (m)')
ax3.set_ylabel('Fraction decaying in tube')
ax3.set_title('Decay Fraction vs Particle Lifetime')
ax3.grid(True, which="both", ls="-", alpha=0.2)
ax3.legend()

# 4. Example decay curve along a specific ray
ax4 = fig.add_subplot(234)
# Find a ray that hits the tube
example_azimuth = 300*np.pi/180#np.pi/4
example_elevation = 62*np.pi/180
dx = np.cos(example_elevation) * np.cos(example_azimuth)
dy = np.cos(example_elevation) * np.sin(example_azimuth)
dz = np.sin(example_elevation)
direction = np.array([dx, dy, dz])

locations, _, _ = mesh.ray.intersects_location(
    ray_origins=[origin],
    ray_directions=[direction]
)

if len(locations) >= 2:
    distances = sorted([np.linalg.norm(loc - origin) for loc in locations])
    entry_dist = distances[0]
    exit_dist = distances[1]
    
    # Plot survival probability along ray
    d = np.linspace(0, exit_dist * 1.5, 1000)
    survival = survival_probability(d, lifetime_distance)
    
    ax4.plot(d, survival, 'b-', linewidth=2, label='Survival probability')
    ax4.axvline(entry_dist, color='green', linestyle='--', label='Tube entry')
    ax4.axvline(exit_dist, color='red', linestyle='--', label='Tube exit')
    ax4.fill_between([entry_dist, exit_dist], [0, 0], [1, 1], alpha=0.2, color='gray', label='Inside tube')
    
    ax4.set_xlabel('Distance (m)')
    ax4.set_ylabel('Survival probability')
    ax4.set_title(f'Example: Ray at Az={np.degrees(example_azimuth):.0f}°, El={np.degrees(example_elevation):.0f}°')
    ax4.grid(True, alpha=0.3)
    ax4.legend()

# 5. 3D visualization
ax5 = fig.add_subplot(235, projection='3d')

# Plot tube outline
sample_indices = np.linspace(0, len(path_3d)-1, 15, dtype=int)
for idx in sample_indices:
    circle_idx = idx * 32
    if circle_idx + 32 <= len(vertices):
        circle_verts = vertices[circle_idx:circle_idx+32:4]
        circle_verts = np.vstack([circle_verts, circle_verts[0]])
        ax5.plot(circle_verts[:, 0], circle_verts[:, 1], circle_verts[:, 2], 
                'b-', alpha=0.3, linewidth=1)

# Plot some example rays with decay probability coloring
n_example_rays = 30
for i in range(n_example_rays):
    idx = np.random.randint(len(results['decay_fractions']))
    if results['decay_fractions'][idx] > 0:
        # Get the corresponding angles
        i_elev = idx // len(results['azimuths'])
        j_azim = idx % len(results['azimuths'])
        azimuth = results['azimuth_all'][idx]
        elevation = results['elevation_all'][idx]
        
        # Create ray
        dx = np.cos(elevation) * np.cos(azimuth)
        dy = np.cos(elevation) * np.sin(azimuth)
        dz = np.sin(elevation)
        
        endpoint = origin + 40 * np.array([dx, dy, dz])
        
        # Color based on decay probability
        color_intensity = results['decay_fractions'][idx]*10
        ax5.plot([origin[0], endpoint[0]], 
                [origin[1], endpoint[1]], 
                [origin[2], endpoint[2]], 
                color=(color_intensity, 0, 1-color_intensity),
                alpha=0.5, linewidth=2)

ax5.scatter(*origin, color='green', s=200, marker='o')
ax5.set_xlabel('X')
ax5.set_ylabel('Y')
ax5.set_zlabel('Z')
ax5.set_title('Sample Rays (colored by decay probability)')

# 6. Summary statistics
ax6 = fig.add_subplot(236, frameon=False)
ax6.axis('off')

summary_text = f"""Summary Statistics:

Particle lifetime: {lifetime_distance:.1f} m

Tube position: z = {Z_POSITION} m
Tube radius: {tube_radius} m

Fraction of particles hitting tube:
  {results['rays_hitting_tube']/results['total_rays']:.4f}

For particles that hit the tube:
  Average decay probability: {results['average_decay_fraction']:.4f}

Overall fraction decaying in tube:
  {results['total_decay_fraction']:.6f}
  ({results['total_decay_fraction']*100:.4f}%)

This means {results['total_decay_fraction']*1e6:.1f} per million
particles decay within the tube volume."""

ax6.text(0.1, 0.9, summary_text, transform=ax6.transAxes, 
         fontsize=10, verticalalignment='top', fontfamily='monospace')

plt.tight_layout()
plt.show()

# Function for easy parameter scanning
def scan_lifetime_values(mesh, origin, lifetime_values):
    """
    Scan multiple lifetime values and return decay fractions
    
    Args:
        mesh: Trimesh object
        origin: Particle origin
        lifetime_values: List or array of lifetime values to test
    
    Returns:
        Array of decay fractions
    """
    decay_fractions = []
    
    print(f"\nScanning {len(lifetime_values)} lifetime values...")
    for lifetime in tqdm(lifetime_values):
        results = calculate_decay_fraction_in_tube(mesh, origin, lifetime,
                                                 azimuth_samples=60, elevation_samples=30,
                                                 show_progress=False)
        decay_fractions.append(results['total_decay_fraction'])
    
    return np.array(decay_fractions)

# Example: Find lifetime for specific decay fraction
def find_lifetime_for_decay_fraction(mesh, origin, target_fraction, 
                                    initial_guess=15.0, tolerance=1e-6):
    """
    Find the lifetime that gives a specific decay fraction in the tube
    
    Args:
        mesh: Trimesh object
        origin: Particle origin
        target_fraction: Desired decay fraction
        initial_guess: Initial lifetime guess
        tolerance: Convergence tolerance
    
    Returns:
        Lifetime value
    """
    from scipy.optimize import brentq
    
    def objective(lifetime):
        results = calculate_decay_fraction_in_tube(mesh, origin, lifetime,
                                                 azimuth_samples=30, elevation_samples=15,
                                                 show_progress=False)
        return results['total_decay_fraction'] - target_fraction
    
    # Find bounds
    lower = initial_guess / 100
    upper = initial_guess * 100
    
    # Use root finding
    optimal_lifetime = brentq(objective, lower, upper, xtol=tolerance)
    
    return optimal_lifetime

# Example usage
# print("\n" + "="*50)
# print("EXAMPLE: Finding lifetime for 1% decay in tube")
# target_lifetime = find_lifetime_for_decay_fraction(mesh, origin, 0.01)
# print(f"Lifetime needed for 1% decay in tube: {target_lifetime:.2f} m")
print("\n" + "="*50)
print("TUBE SHAPE VISUALIZATION")
print("="*50)

def plot_tube_shape(mesh, path_3d, tube_radius, origin=[0,0,0]):
    """
    Create a detailed 3D visualization of the tube shape
    
    Args:
        mesh: Trimesh object
        path_3d: 3D centerline path
        tube_radius: Radius of the tube
        origin: Origin point to show for reference
    """
    fig = plt.figure(figsize=(16, 12))
    
    # 1. Full 3D mesh view
    ax1 = fig.add_subplot(221, projection='3d')
    
    # Plot the mesh faces (subsample for performance)
    if len(mesh.faces) > 1000:
        face_indices = np.random.choice(len(mesh.faces), 1000, replace=False)
        sampled_faces = mesh.faces[face_indices]
    else:
        sampled_faces = mesh.faces
    
    # Plot mesh triangles
    for face in sampled_faces:
        triangle = mesh.vertices[face]
        triangle = np.vstack([triangle, triangle[0]])  # Close the triangle
        ax1.plot(triangle[:, 0], triangle[:, 1], triangle[:, 2], 
                'b-', alpha=0.1, linewidth=0.5)
    
    # Plot centerline
    ax1.plot(path_3d[:, 0], path_3d[:, 1], path_3d[:, 2], 
            'r-', linewidth=3, label='Centerline')
    
    # Plot origin
    ax1.scatter(*origin, color='green', s=200, marker='o', label='Origin (0,0,0)')
    
    # Add coordinate arrows at origin
    arrow_length = 10
    ax1.quiver(origin[0], origin[1], origin[2], arrow_length, 0, 0, 
              color='red', arrow_length_ratio=0.1, label='X')
    ax1.quiver(origin[0], origin[1], origin[2], 0, arrow_length, 0, 
              color='green', arrow_length_ratio=0.1, label='Y')
    ax1.quiver(origin[0], origin[1], origin[2], 0, 0, arrow_length, 
              color='blue', arrow_length_ratio=0.1, label='Z')
    
    ax1.set_xlabel('X (m)')
    ax1.set_ylabel('Y (m)')
    ax1.set_zlabel('Z (m)')
    ax1.set_title('3D Tube Mesh')
    ax1.legend()
    
    # Set equal aspect ratio
    max_range = np.array([
        mesh.vertices[:, 0].max() - mesh.vertices[:, 0].min(),
        mesh.vertices[:, 1].max() - mesh.vertices[:, 1].min(),
        mesh.vertices[:, 2].max() - mesh.vertices[:, 2].min()
    ]).max() / 2.0
    
    mid_x = (mesh.vertices[:, 0].max() + mesh.vertices[:, 0].min()) * 0.5
    mid_y = (mesh.vertices[:, 1].max() + mesh.vertices[:, 1].min()) * 0.5
    mid_z = (mesh.vertices[:, 2].max() + mesh.vertices[:, 2].min()) * 0.5
    
    ax1.set_xlim(mid_x - max_range, mid_x + max_range)
    ax1.set_ylim(mid_y - max_range, mid_y + max_range)
    ax1.set_zlim(0, mesh.vertices[:, 2].max() + 5)
    
    # 2. Top view (XY plane)
    ax2 = fig.add_subplot(222)
    
    # Plot tube outline
    theta = np.linspace(0, 2*np.pi, 50)
    for i in range(0, len(path_3d), 5):
        x_circle = path_3d[i, 0] + tube_radius * np.cos(theta)
        y_circle = path_3d[i, 1] + tube_radius * np.sin(theta)
        ax2.plot(x_circle, y_circle, 'b-', alpha=0.3, linewidth=0.5)
    
    # Plot centerline
    ax2.plot(path_3d[:, 0], path_3d[:, 1], 'r-', linewidth=2, label='Centerline')
    
    # Plot origin
    ax2.scatter(origin[0], origin[1], color='green', s=200, marker='o', label='Origin')
    
    # Draw lines from origin to tube extremes
    tube_points = np.array([[p[0], p[1]] for p in path_3d])
    distances = np.sqrt((tube_points[:, 0] - origin[0])**2 + 
                       (tube_points[:, 1] - origin[1])**2)
    nearest_idx = np.argmin(distances)
    farthest_idx = np.argmax(distances)
    
    ax2.plot([origin[0], path_3d[nearest_idx, 0]], 
            [origin[1], path_3d[nearest_idx, 1]], 
            'g--', alpha=0.5, label=f'Nearest: {distances[nearest_idx]:.1f}m')
    ax2.plot([origin[0], path_3d[farthest_idx, 0]], 
            [origin[1], path_3d[farthest_idx, 1]], 
            'r--', alpha=0.5, label=f'Farthest: {distances[farthest_idx]:.1f}m')
    
    ax2.set_xlabel('X (m)')
    ax2.set_ylabel('Y (m)')
    ax2.set_title('Top View (XY Plane)')
    ax2.grid(True, alpha=0.3)
    ax2.axis('equal')
    ax2.legend()
    
    # 3. Side view (XZ plane)
    ax3 = fig.add_subplot(223)
    
    # Plot tube profile
    ax3.fill_between(path_3d[:, 0], 
                    path_3d[:, 2] - tube_radius,
                    path_3d[:, 2] + tube_radius, 
                    alpha=0.3, color='blue', label='Tube volume')
    
    # Plot centerline
    ax3.plot(path_3d[:, 0], path_3d[:, 2], 'r-', linewidth=2, label='Centerline')
    
    # Plot origin and height reference
    ax3.scatter(origin[0], origin[2], color='green', s=200, marker='o', label='Origin')
    ax3.axhline(y=origin[2], color='green', linestyle=':', alpha=0.5)
    ax3.axhline(y=Z_POSITION, color='red', linestyle=':', alpha=0.5, label=f'z={Z_POSITION}m')
    
    # Add distance annotations
    ax3.annotate(f'{Z_POSITION}m', xy=(origin[0]-5, Z_POSITION), 
                xytext=(origin[0]-10, Z_POSITION), fontsize=10)
    
    ax3.set_xlabel('X (m)')
    ax3.set_ylabel('Z (m)')
    ax3.set_title('Side View (XZ Plane)')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    ax3.set_ylim(-5, Z_POSITION + 10)
    
    # 4. End view (YZ plane) - looking along X axis
    ax4 = fig.add_subplot(224)
    
    # Find a few cross-sections along the tube
    x_positions = np.linspace(mesh.bounds[0][0], mesh.bounds[1][0], 5)
    
    for x_pos in x_positions:
        # Find closest path point to this x position
        idx = np.argmin(np.abs(path_3d[:, 0] - x_pos))
        
        # Draw circle at this position
        theta = np.linspace(0, 2*np.pi, 50)
        y_circle = path_3d[idx, 1] + tube_radius * np.cos(theta)
        z_circle = path_3d[idx, 2] + tube_radius * np.sin(theta)
        
        alpha = 0.3 + 0.5 * (x_pos - mesh.bounds[0][0]) / (mesh.bounds[1][0] - mesh.bounds[0][0])
        ax4.plot(y_circle, z_circle, 'b-', alpha=alpha, linewidth=1)
    
    # Plot origin
    ax4.scatter(origin[1], origin[2], color='green', s=200, marker='o', label='Origin')
    
    # Draw sight lines from origin
    angles = np.linspace(0, 2*np.pi, 8)
    for angle in angles:
        y_end = origin[1] + 50 * np.cos(angle)
        z_end = origin[2] + 50 * np.sin(angle)
        ax4.plot([origin[1], y_end], [origin[2], z_end], 'g-', alpha=0.1, linewidth=0.5)
    
    ax4.set_xlabel('Y (m)')
    ax4.set_ylabel('Z (m)')
    ax4.set_title('End View (YZ Plane) - Looking Along X')
    ax4.grid(True, alpha=0.3)
    ax4.axis('equal')
    ax4.legend()
    
    plt.tight_layout()
    plt.show()
    
    # Print tube statistics
    print("\nTube Statistics:")
    print(f"Centerline length: {np.sum(np.linalg.norm(np.diff(path_3d, axis=0), axis=1)):.1f} m")
    print(f"Tube radius: {tube_radius} m")
    print(f"Tube volume: {mesh.volume:.1f} m³")
    print(f"Tube surface area: {mesh.area:.1f} m²")
    print(f"Height above origin: {Z_POSITION} m")
    print(f"Horizontal distance range: {distances.min():.1f} - {distances.max():.1f} m")

# Plot the tube
plot_tube_shape(mesh, path_3d, tube_radius, origin)

# Create an interactive 3D plot using matplotlib
def plot_tube_interactive(mesh, path_3d, origin=[0,0,0]):
    """
    Create an interactive 3D plot of the tube that can be rotated
    """
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot tube surface using mesh vertices
    # Create rings at regular intervals
    n_rings = 30
    for i in np.linspace(0, len(path_3d)-1, n_rings, dtype=int):
        idx = i * 32
        if idx + 32 <= len(vertices):
            ring = vertices[idx:idx+32]
            ring = np.vstack([ring, ring[0]])  # Close the ring
            ax.plot(ring[:, 0], ring[:, 1], ring[:, 2], 'b-', alpha=0.5, linewidth=1)
    
    # Plot longitudinal lines
    for j in range(0, 32, 4):  # Every 4th vertex for clarity
        long_line = []
        for i in range(len(path_3d)):
            idx = i * 32 + j
            if idx < len(vertices):
                long_line.append(vertices[idx])
        if long_line:
            long_line = np.array(long_line)
            ax.plot(long_line[:, 0], long_line[:, 1], long_line[:, 2], 
                   'b-', alpha=0.3, linewidth=0.5)
    
    # Plot centerline
    ax.plot(path_3d[:, 0], path_3d[:, 1], path_3d[:, 2], 
           'r-', linewidth=3, label='Centerline')
    
    # Plot origin and axes
    ax.scatter(*origin, color='green', s=300, marker='o', 
              edgecolors='black', linewidth=2, label='Origin')
    
    # Coordinate axes
    axis_length = 15
    ax.quiver(origin[0], origin[1], origin[2], axis_length, 0, 0, 
             color='red', arrow_length_ratio=0.1, linewidth=2)
    ax.quiver(origin[0], origin[1], origin[2], 0, axis_length, 0, 
             color='green', arrow_length_ratio=0.1, linewidth=2)
    ax.quiver(origin[0], origin[1], origin[2], 0, 0, axis_length, 
             color='blue', arrow_length_ratio=0.1, linewidth=2)
    
    ax.text(origin[0] + axis_length, origin[1], origin[2], 'X', fontsize=12)
    ax.text(origin[0], origin[1] + axis_length, origin[2], 'Y', fontsize=12)
    ax.text(origin[0], origin[1], origin[2] + axis_length, 'Z', fontsize=12)
    
    # Set labels and title
    ax.set_xlabel('X (m)', fontsize=10)
    ax.set_ylabel('Y (m)', fontsize=10)
    ax.set_zlabel('Z (m)', fontsize=10)
    ax.set_title('3D Tube Visualization (Rotate with mouse)', fontsize=14)
    
    # Set viewing angle
    ax.view_init(elev=20, azim=45)
    
    # Set aspect ratio
    ax.set_box_aspect([1,1,0.5])
    
    ax.legend(fontsize=10)
    plt.show()

# Create the interactive plot
print("\nCreating interactive 3D visualization...")
plot_tube_interactive(mesh, path_3d, origin)
