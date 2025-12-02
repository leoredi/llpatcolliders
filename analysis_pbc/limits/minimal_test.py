#!/usr/bin/env python
"""Minimal test to process one mass point"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Setup paths
sys.path.insert(0, str(Path(__file__).parent.parent))

print("="*60)
print("MINIMAL TEST: Processing 2.6 GeV muon")
print("="*60)

# Import modules
print("\n[1/5] Importing modules...")
from geometry.per_parent_efficiency import build_drainage_gallery_mesh, preprocess_hnl_csv
from models.hnl_model_hnlcalc import HNLModel
from config.production_xsecs import get_parent_sigma_pb
from limits.u2_limit_calculator import expected_signal_events
print("✓ Imports successful")

# Load simulation
print("\n[2/5] Loading simulation CSV...")
sim_csv = Path("../output/csv/simulation_new/HNL_2p60GeV_muon_beauty.csv")
geom_csv = Path("../output/csv/geometry/HNL_2p60GeV_muon_geom.csv")

if not sim_csv.exists():
    print(f"ERROR: {sim_csv} not found!")
    sys.exit(1)

print(f"✓ Found {sim_csv}")

# Process geometry
print("\n[3/5] Processing geometry...")
if geom_csv.exists():
    print(f"✓ Loading cached geometry from {geom_csv}")
    geom_df = pd.read_csv(geom_csv)
else:
    print("✓ Computing geometry (this will take ~30 sec)...")
    mesh = build_drainage_gallery_mesh()
    geom_df = preprocess_hnl_csv(sim_csv, mesh)
    geom_df.to_csv(geom_csv, index=False)
    print(f"✓ Saved geometry to {geom_csv}")

print(f"✓ Loaded {len(geom_df)} HNLs, {geom_df['hits_tube'].sum()} hit detector")

# Test signal calculation
print("\n[4/5] Testing signal calculation...")
mass_GeV = 2.6
eps2_test = 1e-6
benchmark = "010"

N_sig = expected_signal_events(
    geom_df=geom_df,
    mass_GeV=mass_GeV,
    eps2=eps2_test,
    benchmark=benchmark,
    lumi_fb=3000.0
)

print(f"✓ N_sig = {N_sig:.2f} events at |Uμ|² = {eps2_test:.1e}")

# Scan |U|^2
print("\n[5/5] Scanning |U|² (10 points)...")
eps2_scan = np.logspace(-10, -4, 10)
results = []

for eps2 in eps2_scan:
    N = expected_signal_events(geom_df, mass_GeV, eps2, benchmark, 3000.0)
    results.append((eps2, N))
    if N >= 3:
        print(f"  |U|² = {eps2:.2e} → N = {N:.1f} {'✓ EXCLUDED' if N >= 3 else ''}")

print("\n" + "="*60)
print("MINIMAL TEST COMPLETE!")
print("="*60)
