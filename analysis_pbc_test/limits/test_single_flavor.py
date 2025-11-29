"""
Test script to run just muon flavor
"""
import sys
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from limits.u2_limit_calculator import run_reach_scan

if __name__ == "__main__":
    L_HL_LHC_FB = 3000.0
    N_CORES = 4

    print("Starting muon scan...")
    run_reach_scan("muon", "010", L_HL_LHC_FB, n_jobs=N_CORES)
    print("Muon scan complete!")
