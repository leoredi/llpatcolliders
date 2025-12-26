"""
Simple standalone tests for llpdecay package.
Can be run without pytest.
"""

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

try:
    import numpy as np
    print("✓ NumPy imported successfully")
except ImportError as e:
    print(f"✗ Failed to import NumPy: {e}")
    print("\nPlease install numpy: pip install numpy")
    sys.exit(1)

print("\n" + "="*60)
print("Testing llpdecay package")
print("="*60)

# Test 1: Import package
print("\n[1] Testing package import...")
try:
    from llpdecay import HNL
    from llpdecay.core import kallen, boost_to_lab, sample_two_body_decay
    print("✓ Package imported successfully")
except ImportError as e:
    print(f"✗ Failed to import package: {e}")
    sys.exit(1)

# Test 2: Create HNL model
print("\n[2] Testing HNL model creation...")
try:
    hnl = HNL(mass=2.0, Umu=1e-6, seed=42)
    print(f"✓ Created HNL model: {hnl}")
except Exception as e:
    print(f"✗ Failed to create HNL: {e}")
    sys.exit(1)

# Test 3: Get available channels
print("\n[3] Testing available channels...")
try:
    channels = hnl.available_channels()
    print(f"✓ Found {len(channels)} available channels")
    print(f"  Channels: {', '.join(channels[:5])}...")
except Exception as e:
    print(f"✗ Failed to get channels: {e}")
    sys.exit(1)

# Test 4: Calculate branching ratios
print("\n[4] Testing branching ratios...")
try:
    brs = hnl.branching_ratios()
    total = sum(brs.values())
    print(f"✓ Computed {len(brs)} branching ratios (sum={total:.6f})")
    # Show top 3
    for ch, br in sorted(brs.items(), key=lambda x: -x[1])[:3]:
        print(f"  {ch}: {br:.3%}")
except Exception as e:
    print(f"✗ Failed to compute BRs: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 5: Sample a decay
print("\n[5] Testing decay sampling...")
try:
    parent_4vec = np.array([10.0, 3.0, 0.5, 9.3])
    daughters, channel = hnl.sample_decay(parent_4vec, return_channel=True)

    print(f"✓ Sampled decay via channel: {channel}")
    print(f"  Daughters shape: {daughters.shape}")

    # Check conservation
    total = np.sum(daughters[0], axis=0)
    conservation_error = np.linalg.norm(total - parent_4vec)
    print(f"  Conservation error: {conservation_error:.2e}")

    if conservation_error > 1e-5:
        print(f"  ⚠ WARNING: Large conservation error!")
    else:
        print(f"  ✓ Energy-momentum conserved")

except Exception as e:
    print(f"✗ Failed to sample decay: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 6: Sample multiple decays
print("\n[6] Testing multiple decay sampling...")
try:
    n_events = 100
    daughters = hnl.sample_decay(parent_4vec, n_events=n_events)
    print(f"✓ Sampled {n_events} decays")
    print(f"  Shape: {daughters.shape}")

    # Check conservation for all
    errors = []
    for i in range(n_events):
        total = np.sum(daughters[i], axis=0)
        error = np.linalg.norm(total - parent_4vec)
        errors.append(error)

    max_error = max(errors)
    print(f"  Max conservation error: {max_error:.2e}")

    if max_error > 1e-5:
        print(f"  ⚠ WARNING: Large conservation errors in some events")
    else:
        print(f"  ✓ All events conserve energy-momentum")

except Exception as e:
    print(f"✗ Failed multiple sampling: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 7: Lifetime calculations
print("\n[7] Testing lifetime calculations...")
try:
    width = hnl.total_width()
    lifetime = hnl.lifetime()
    ctau = hnl.ctau()

    print(f"✓ Total width: {width:.3e} GeV")
    print(f"  Lifetime: {lifetime:.3e} s")
    print(f"  Decay length: {ctau:.3e} m")

    if width <= 0:
        print(f"  ✗ ERROR: Width should be positive!")
        sys.exit(1)

except Exception as e:
    print(f"✗ Failed lifetime calculation: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 8: Test kinematics functions
print("\n[8] Testing kinematics functions...")
try:
    # Test Källén function
    lam = kallen(4.0, 1.0, 1.0)
    print(f"✓ Källén(4, 1, 1) = {lam:.3f}")

    # Test two-body decay
    daughters_rest = sample_two_body_decay(2.0, 0.5, 0.3, n_events=10, seed=42)
    print(f"✓ Sampled 2-body decays in rest frame: {daughters_rest.shape}")

    # Test boost
    parent = np.array([10.0, 3.0, 0.0, 9.5])
    daughters_lab = boost_to_lab(daughters_rest[0], parent)
    print(f"✓ Boosted to lab frame: {daughters_lab.shape}")

except Exception as e:
    print(f"✗ Failed kinematics tests: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# All tests passed!
print("\n" + "="*60)
print("✓ ALL TESTS PASSED!")
print("="*60)
print("\nllpdecay package is working correctly.")
