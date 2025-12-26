"""
Demonstration of Phase 2 features.

Shows:
1. Three-body decay sampling
2. ALP model usage
3. Form factors and radiative corrections
4. HNLCalc validation
"""

import numpy as np
import sys
from pathlib import Path

# Add package to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from llpdecay import HNL, ALP
from llpdecay.decays.three_body import sample_three_body_decay, hnl_three_body_leptonic_me
from llpdecay.advanced import form_factor_pion, qed_correction_lepton_pair
from llpdecay.validation import print_comparison_table


def demo_three_body_decays():
    """Demonstrate 3-body decay physics."""
    print("=" * 70)
    print("DEMO 1: Three-Body Decay Sampling")
    print("=" * 70)

    # Create HNL with 3-body channels
    hnl = HNL(mass=2.0, Umu=1e-6, seed=42)

    print(f"\nHNL: {hnl}")
    print(f"Available channels: {len(hnl.available_channels())}")

    # Get 3-body channels
    channels = hnl.available_channels()
    three_body = [ch for ch in channels if 'nu_' in ch and '_pi0' not in ch]

    if three_body:
        print(f"\n3-body channels: {three_body}")

        # Sample a 3-body decay
        parent_4vec = np.array([10.0, 3.0, 0.0, 9.5])
        print(f"\nParent 4-vector: {parent_4vec}")

        for ch in three_body[:2]:  # Show first 2
            print(f"\n  Channel: {ch}")

            daughters, _ = hnl.sample_decay(
                parent_4vec,
                channel=ch,
                n_events=5,
                return_channel=True
            )

            print(f"  Generated {len(daughters)} events")
            print(f"  Daughter 1: {daughters[0, 0]}")
            print(f"  Daughter 2: {daughters[0, 1]}")
            print(f"  Daughter 3: {daughters[0, 2]}")

            # Check conservation
            total = np.sum(daughters[0], axis=0)
            error = np.linalg.norm(total - parent_4vec)
            print(f"  Conservation error: {error:.2e}")
    else:
        print("\nNo 3-body channels available at this mass")

    print()


def demo_alp_model():
    """Demonstrate ALP physics."""
    print("=" * 70)
    print("DEMO 2: Axion-Like Particle (ALP) Model")
    print("=" * 70)

    # Photophilic ALP
    print("\n[A] Photophilic ALP (couples to photons)")
    alp_photon = ALP(mass=0.5, g_agg=1e-4, seed=42)
    print(f"  {alp_photon}")

    brs = alp_photon.branching_ratios()
    print(f"\n  Branching ratios:")
    for ch, br in sorted(brs.items(), key=lambda x: -x[1]):
        print(f"    {ch}: {br:.3%}")

    print(f"\n  Lifetime: {alp_photon.lifetime():.3e} s")
    print(f"  Decay length: {alp_photon.ctau():.3e} m")

    # Sample decay
    parent = np.array([5.0, 2.0, 0.0, 4.5])
    daughters, ch = alp_photon.sample_decay(parent, n_events=1, return_channel=True)
    print(f"\n  Sampled decay: {ch[0]}")
    print(f"  Daughter PDGs: {alp_photon.get_daughter_pdgs(ch[0])}")

    # Leptophilic ALP
    print("\n[B] Leptophilic ALP (couples to leptons)")
    alp_lepton = ALP(mass=1.0, g_agg=0, f_a=1e8, c_e=1.0, c_mu=1.0, seed=42)
    print(f"  {alp_lepton}")

    brs = alp_lepton.branching_ratios()
    print(f"\n  Branching ratios:")
    for ch, br in sorted(brs.items(), key=lambda x: -x[1]):
        print(f"    {ch}: {br:.3%}")

    print()


def demo_form_factors():
    """Demonstrate form factors and corrections."""
    print("=" * 70)
    print("DEMO 3: Form Factors and Radiative Corrections")
    print("=" * 70)

    # Pion form factor
    print("\n[A] Pion Form Factor")
    q2_values = [0.0, 0.1, 0.3, 0.5]
    print("  q² (GeV²)  |  f₊(q²)")
    print("  " + "-" * 25)
    for q2 in q2_values:
        ff = form_factor_pion(q2)
        print(f"  {q2:8.2f}   |  {ff:6.3f}")

    # QED corrections
    print("\n[B] QED Corrections for e⁺e⁻ Production")
    m_e = 0.000511
    s_values = [0.1, 1.0, 10.0]
    print("  s (GeV²)   |  QED Correction")
    print("  " + "-" * 30)
    for s in s_values:
        corr = qed_correction_lepton_pair(s, m_e)
        delta = (corr - 1.0) * 100
        print(f"  {s:8.2f}   |  {corr:.4f}  (+{delta:.2f}%)")

    print()


def demo_hnlcalc_validation():
    """Demonstrate HNLCalc validation."""
    print("=" * 70)
    print("DEMO 4: HNLCalc Validation (if available)")
    print("=" * 70)

    print("\nComparing llpdecay with HNLCalc...")

    try:
        # This will attempt to import HNLCalc
        print_comparison_table(mass=2.0, Umu=1e-6)
    except Exception as e:
        print(f"\nHNLCalc not available: {e}")
        print("Skipping validation demo.")

        # Show internal calculation instead
        print("\nShowing llpdecay internal calculation:")
        hnl = HNL(mass=2.0, Umu=1e-6)
        brs = hnl.branching_ratios()

        print(f"\nBranching Ratios (m_N = {hnl.mass} GeV, Uμ² = {hnl.Umu:.1e}):")
        print("-" * 50)
        for ch, br in sorted(brs.items(), key=lambda x: -x[1])[:10]:
            print(f"  {ch:<15} {br:>10.3%}")

    print()


def demo_advanced_usage():
    """Demonstrate advanced features."""
    print("=" * 70)
    print("DEMO 5: Advanced Usage Examples")
    print("=" * 70)

    # Majorana vs Dirac comparison
    print("\n[A] Majorana vs Dirac HNL")

    hnl_maj = HNL(mass=2.0, Umu=1e-6, is_majorana=True)
    hnl_dir = HNL(mass=2.0, Umu=1e-6, is_majorana=False)

    width_maj = hnl_maj.total_width()
    width_dir = hnl_dir.total_width()

    print(f"  Majorana width: {width_maj:.3e} GeV")
    print(f"  Dirac width:    {width_dir:.3e} GeV")
    print(f"  Ratio:          {width_maj/width_dir:.2f}")

    # Flavor mixing comparison
    print("\n[B] Flavor Mixing Effects")

    mixings = [
        ("e-only", {'Ue': 1e-6, 'Umu': 0, 'Utau': 0}),
        ("μ-only", {'Ue': 0, 'Umu': 1e-6, 'Utau': 0}),
        ("equal", {'Ue': 1e-6/3, 'Umu': 1e-6/3, 'Utau': 1e-6/3}),
    ]

    for name, mixing in mixings:
        hnl = HNL(mass=2.0, **mixing)
        brs = hnl.branching_ratios()

        print(f"\n  {name} mixing:")
        # Show top 3 channels
        for ch, br in sorted(brs.items(), key=lambda x: -x[1])[:3]:
            print(f"    {ch:<12} {br:>6.1%}")

    print()


def main():
    """Run all demos."""
    print("\n")
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 68 + "║")
    print("║" + "  LLPDecay Phase 2 Features Demonstration".center(68) + "║")
    print("║" + " " * 68 + "║")
    print("╚" + "=" * 68 + "╝")
    print()

    demo_three_body_decays()
    demo_alp_model()
    demo_form_factors()
    demo_hnlcalc_validation()
    demo_advanced_usage()

    print("=" * 70)
    print("All demos completed successfully!")
    print("=" * 70)
    print()


if __name__ == '__main__':
    main()
