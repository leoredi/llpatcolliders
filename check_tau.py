import pandas as pd
import numpy as np

# Read the CSV
df = pd.read_csv('output/csv/hnlLL_m39GeVLLP.csv', sep=',\s*', engine='python')

# tau is in mm/c units in PYTHIA
# Convert to proper lifetime in seconds: tau_s = tau_mm / c_mm/s
c_mm_per_s = 299792458.0 * 1000.0  # speed of light in mm/s

tau_mm_c = df['tau'].values
tau_seconds = tau_mm_c / c_mm_per_s

print(f"Number of particles: {len(df)}")
print(f"\nProper decay time τ (mm/c):")
print(f"  Mean: {tau_mm_c.mean():.2f} mm/c")
print(f"  Median: {np.median(tau_mm_c):.2f} mm/c")
print(f"  Std: {tau_mm_c.std():.2f} mm/c")
print(f"  Min: {tau_mm_c.min():.2f} mm/c")
print(f"  Max: {tau_mm_c.max():.2f} mm/c")

print(f"\nProper lifetime τ (seconds):")
print(f"  Mean: {tau_seconds.mean():.3e} s  ({tau_seconds.mean()*1e9:.2f} ns)")
print(f"  Median: {np.median(tau_seconds):.3e} s  ({np.median(tau_seconds)*1e9:.2f} ns)")
print(f"  Expected from config: 1e4 mm/c = {1e4/c_mm_per_s:.3e} s = {1e4/c_mm_per_s*1e9:.2f} ns")

print(f"\nProper decay length c*τ (meters):")
ctau_meters = tau_mm_c / 1000.0  # convert mm to m
print(f"  Mean: {ctau_meters.mean():.2f} m")
print(f"  Median: {np.median(ctau_meters):.2f} m")
print(f"  Expected from config: 10 m")

# Calculate actual decay length in lab frame
print(f"\nLab frame decay distance:")
decay_dist = np.sqrt((df['xDec'] - df['xProd'])**2 +
                     (df['yDec'] - df['yProd'])**2 +
                     (df['zDec'] - df['zProd'])**2)
print(f"  Mean: {decay_dist.mean()/1000:.2f} m")
print(f"  Median: {np.median(decay_dist)/1000:.2f} m")

# Calculate gamma factor
gamma = df['momentum'] / df['mass']
print(f"\nGamma factor (boost):")
print(f"  Mean: {gamma.mean():.2f}")
print(f"  Median: {np.median(gamma):.2f}")

# Expected lab-frame decay length = gamma * c * tau
expected_lab_decay = gamma * tau_mm_c / 1000.0  # in meters
print(f"\nExpected lab-frame decay length γ*c*τ:")
print(f"  Mean: {expected_lab_decay.mean():.2f} m")
print(f"  Median: {np.median(expected_lab_decay):.2f} m")
