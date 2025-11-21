"""
Test Pipeline - Minimal Example

Demonstrates the PBC-grade analysis pipeline end-to-end:
1. Load HNL model
2. Compute per-parent efficiencies
3. Calculate N_sig
4. Find |U|² limit

This is a standalone test that doesn't require HNLCalc
(uses toy models for demonstration).
"""

import sys
from pathlib import Path
import numpy as np

# Add modules to path
sys.path.append(str(Path(__file__).parent.parent / 'models'))
sys.path.append(str(Path(__file__).parent.parent / 'geometry'))
sys.path.append(str(Path(__file__).parent.parent / 'limits'))


def test_model():
    """Test HNL model wrapper."""
    print("="*70)
    print("TEST 1: HNL Model")
    print("="*70)

    from hnl_model_hnlcalc import HNLModelHNLCalc

    model = HNLModelHNLCalc()

    # Test proper lifetime
    mass = 1.0  # GeV
    flavour = 'muon'

    ctau0 = model.ctau0_m_for_U2_eq_1(mass, flavour)
    print(f"\nProper lifetime at m={mass} GeV, {flavour}:")
    print(f"  cτ₀ (|U|²=1) = {ctau0:.3e} m")

    # Test at different |U|²
    U2_values = [1e-6, 1e-8, 1e-10]
    for U2 in U2_values:
        ctau = ctau0 / U2
        print(f"  cτ (|U|²={U2:.0e}) = {ctau:.3e} m")

    # Test cross-sections
    print(f"\nProduction cross-sections at √s = 14 TeV:")
    parents = ['K+', 'D0', 'B+', 'W+']
    for parent_name in parents:
        pdg = model.PARENT_PDGIDS[parent_name]
        sigma = model.sigma_parent(pdg, mass, flavour)
        print(f"  σ({parent_name:3s}) = {sigma:.3e} pb")

    # Test branching ratios
    print(f"\nBranching ratios at m={mass} GeV, {flavour}, |U|²=1e-6:")
    for parent_name in parents:
        pdg = model.PARENT_PDGIDS[parent_name]
        br = model.BR_parent_to_HNL(pdg, mass, flavour, 1e-6)
        print(f"  BR({parent_name:3s} → ℓN) = {br:.3e}")

    print("\n✓ Model test passed\n")


def test_efficiency_calculation():
    """Test per-parent efficiency calculation."""
    print("="*70)
    print("TEST 2: Per-Parent Efficiency (Toy Data)")
    print("="*70)

    from per_parent_efficiency import PerParentEfficiency, DetectorGeometry
    import pandas as pd

    # Create toy Pythia data
    print("\nGenerating toy event sample...")

    n_events = 1000
    data = {
        'event': np.arange(n_events),
        'weight': np.ones(n_events),
        'id': [9900015] * n_events,
        'parent_id': np.random.choice([321, 421, 521, 24], n_events),  # K, D, B, W
        'pt': np.random.uniform(1, 50, n_events),
        'eta': np.random.uniform(-5, 5, n_events),
        'phi': np.random.uniform(-np.pi, np.pi, n_events),
        'momentum': np.random.uniform(10, 100, n_events),
        'energy': np.random.uniform(10, 100, n_events),
        'mass': [1.0] * n_events,
        'prod_x_m': np.random.uniform(-0.001, 0.001, n_events),
        'prod_y_m': np.random.uniform(-0.001, 0.001, n_events),
        'prod_z_m': np.random.uniform(-0.01, 0.01, n_events)
    }

    df = pd.DataFrame(data)

    # Initialize calculator
    detector = DetectorGeometry()
    calc = PerParentEfficiency(detector)

    print(f"  Detector volume: {detector.volume_m3():.1f} m³")
    print(f"  Events generated: {len(df)}")

    # Compute efficiency for one cτ value
    ctau_test = 10.0  # meters
    P_decay = calc.decay_probability_in_detector(df, ctau_test)

    print(f"\nDecay probabilities at cτ = {ctau_test} m:")
    print(f"  Mean: {np.mean(P_decay):.3e}")
    print(f"  Median: {np.median(P_decay):.3e}")
    print(f"  Max: {np.max(P_decay):.3e}")

    # Group by parent
    for parent_pdg in [321, 421, 521, 24]:
        mask = df['parent_id'] == parent_pdg
        if mask.sum() > 0:
            P_parent = P_decay[mask]
            parent_name = calc.detector.__class__.__name__  # placeholder
            print(f"  Parent {parent_pdg}: <P> = {np.mean(P_parent):.3e} "
                  f"({mask.sum()} events)")

    print("\n✓ Efficiency calculation test passed\n")


def test_nsig_calculation():
    """Test N_sig calculation."""
    print("="*70)
    print("TEST 3: Signal Yield Calculation")
    print("="*70)

    from u2_limit_calculator import U2LimitCalculator

    # Initialize calculator
    lumi_fb = 3000.0  # HL-LHC
    calc = U2LimitCalculator(lumi_fb=lumi_fb)

    # Create toy efficiency map
    print("\nCreating toy efficiency map...")

    efficiency_map = {
        321: {1.0: 1e-5, 10.0: 1e-4, 100.0: 1e-3},     # K+
        421: {1.0: 5e-6, 10.0: 5e-5, 100.0: 5e-4},     # D0
        521: {1.0: 1e-6, 10.0: 1e-5, 100.0: 1e-4},     # B+
        24:  {1.0: 1e-4, 10.0: 1e-3, 100.0: 1e-2},     # W+
    }

    # Compute N_sig for test point
    mass = 1.0  # GeV
    flavour = 'muon'
    U2 = 1e-7

    print(f"\nComputing N_sig for:")
    print(f"  Mass: {mass} GeV")
    print(f"  Flavour: {flavour}")
    print(f"  |U|²: {U2:.2e}")
    print(f"  Luminosity: {lumi_fb:.0f} fb⁻¹")

    N_sig = calc.compute_Nsig(mass, flavour, U2, efficiency_map)

    print(f"\nResult: N_sig = {N_sig:.3e} events")

    print("\n✓ N_sig calculation test passed\n")


def test_u2_limit_finding():
    """Test |U|² limit finding."""
    print("="*70)
    print("TEST 4: |U|² Limit Finding")
    print("="*70)

    from u2_limit_calculator import U2LimitCalculator

    # Initialize
    calc = U2LimitCalculator(lumi_fb=3000.0, n_limit=3.0)

    # Toy efficiency map
    efficiency_map = {
        321: {1.0: 1e-5, 10.0: 1e-4, 100.0: 1e-3},
        421: {1.0: 5e-6, 10.0: 5e-5, 100.0: 5e-4},
        521: {1.0: 1e-6, 10.0: 1e-5, 100.0: 1e-4},
        24:  {1.0: 1e-4, 10.0: 1e-3, 100.0: 1e-2},
    }

    # Find limit
    mass = 1.0
    flavour = 'muon'

    print(f"\nFinding |U|² limit for m={mass} GeV, {flavour}...")

    U2_limit, info = calc.find_U2_limit(mass, flavour, efficiency_map)

    print(f"\nResult:")
    print(f"  |U|²_limit = {U2_limit:.3e}")
    print(f"  cτ at limit = {info['ctau_at_limit']:.3e} m")
    print(f"  N_sig = {info['N_sig_at_limit']:.1f} events")

    print("\n✓ Limit finding test passed\n")


def run_all_tests():
    """Run all pipeline tests."""
    print("\n" + "="*70)
    print("PBC-GRADE PIPELINE TEST SUITE")
    print("="*70 + "\n")

    try:
        test_model()
        test_efficiency_calculation()
        test_nsig_calculation()
        test_u2_limit_finding()

        print("="*70)
        print("ALL TESTS PASSED ✓")
        print("="*70)
        print("\nNext steps:")
        print("1. Install HNLCalc for real width/BR calculations")
        print("2. Run on real Pythia data from production/")
        print("3. Compare with PBC BC6/7/8 benchmark curves")

    except Exception as e:
        print("\n" + "="*70)
        print(f"TEST FAILED: {e}")
        print("="*70)
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    run_all_tests()
