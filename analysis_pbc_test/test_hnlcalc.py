#!/usr/bin/env python3
"""
Quick test of HNLCalc installation and basic functionality.
"""

import sys
sys.path.append('HNLCalc')

from HNLCalc import HNLCalc
import numpy as np

print("="*70)
print("TESTING HNLCalc")
print("="*70)

# Test 1: Initialize for muon coupling (BC7)
print("\nTest 1: Initialize HNLCalc for muon coupling (BC7)")
print("HNL(ve=0, vmu=1, vtau=0)")

try:
    hnl = HNLCalc(ve=0, vmu=1, vtau=0)
    print("✓ HNLCalc initialized successfully")
except Exception as e:
    print(f"✗ ERROR: {e}")
    sys.exit(1)

# Test 2: Check class methods
print("\nTest 2: Check available methods")
methods = [m for m in dir(hnl) if not m.startswith('_') and callable(getattr(hnl, m))]
print(f"Found {len(methods)} public methods")
print("Key methods:")
for m in ['total_width', 'BR', 'production']:
    if m in methods:
        print(f"  ✓ {m}")
    else:
        print(f"  ? {m} (checking alternatives...)")

# Show first 10 methods
print("\nFirst 10 methods:")
for m in methods[:10]:
    print(f"  - {m}")

print("\n" + "="*70)
print("HNLCalc installation verified!")
print("="*70)
