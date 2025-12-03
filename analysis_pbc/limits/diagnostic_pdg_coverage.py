"""
limits/diagnostic_pdg_coverage.py

Diagnostic tool to check for PDG coverage mismatches between:
1. Simulation CSV parent_id values (from Pythia)
2. HNLCalc production channels (from HNLModel.production_brs())
3. Cross-section lookup table (from production_xsecs.py)

This helps identify silent data loss when unknown PDG codes appear.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Set
import pandas as pd
import numpy as np

# Add analysis_pbc to path
THIS_FILE = Path(__file__).resolve()
ANALYSIS_DIR = THIS_FILE.parent
ANALYSIS_ROOT = ANALYSIS_DIR.parent
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from models.hnl_model_hnlcalc import HNLModel
from config.production_xsecs import get_parent_sigma_pb

REPO_ROOT = ANALYSIS_ROOT.parent
SIM_DIR = REPO_ROOT / "output" / "csv" / "simulation"


def get_pdg_from_csvs(csv_dir: Path = SIM_DIR, max_files: int = 10) -> Set[int]:
    """
    Extract unique parent PDG codes from simulation CSV files.

    Parameters
    ----------
    csv_dir : Path
        Directory containing simulation CSV files
    max_files : int
        Maximum number of files to scan (to avoid reading all 102 CSVs)

    Returns
    -------
    Set[int]
        Set of unique absolute parent PDG codes found in CSV files
    """
    csv_files = list(csv_dir.glob("HNL_mass_*_Meson.csv"))[:max_files]

    if not csv_files:
        print(f"[WARN] No CSV files found in {csv_dir}")
        return set()

    unique_pdgs = set()

    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file)
            if "parent_id" in df.columns:
                pdgs = df["parent_id"].dropna().to_numpy()
                pdgs = pdgs[np.isfinite(pdgs)].astype(int)
                unique_pdgs.update(abs(p) for p in pdgs)
        except Exception as e:
            print(f"[WARN] Failed to read {csv_file.name}: {e}")
            continue

    return unique_pdgs


def get_pdg_from_hnlcalc(mass_GeV: float = 1.0, benchmark: str = "010") -> Set[int]:
    """
    Extract parent PDG codes from HNLCalc production channels.

    Parameters
    ----------
    mass_GeV : float
        HNL mass to test (default 1.0 GeV, should have many channels)
    benchmark : str
        Coupling benchmark: "100" (electron), "010" (muon), "001" (tau)

    Returns
    -------
    Set[int]
        Set of parent PDG codes that HNLCalc knows about
    """
    # Use small coupling to avoid numerical issues
    if benchmark == "100":
        model = HNLModel(mass_GeV=mass_GeV, Ue2=1e-6, Umu2=0.0, Utau2=0.0)
    elif benchmark == "010":
        model = HNLModel(mass_GeV=mass_GeV, Ue2=0.0, Umu2=1e-6, Utau2=0.0)
    elif benchmark == "001":
        model = HNLModel(mass_GeV=mass_GeV, Ue2=0.0, Umu2=0.0, Utau2=1e-6)
    else:
        raise ValueError(f"Unknown benchmark: {benchmark}")

    br_dict = model.production_brs()

    return set(br_dict.keys())


def get_pdg_from_xsec_table() -> Set[int]:
    """
    Extract parent PDG codes that have cross-section entries.

    Returns
    -------
    Set[int]
        Set of parent PDG codes with known cross-sections
    """
    # Test common PDG codes used in HNL physics
    test_pdgs = [
        321,   # K+
        130,   # KL0
        310,   # KS0
        411,   # D+
        421,   # D0
        431,   # Ds+
        4122,  # Λc+
        511,   # B0
        521,   # B+
        531,   # Bs0
        541,   # Bc+
        5122,  # Λb0
        5232,  # Σb
        5332,  # Ωb
        15,    # τ
    ]

    known_pdgs = set()

    for pdg in test_pdgs:
        sigma = get_parent_sigma_pb(pdg)
        if sigma > 0.0:
            known_pdgs.add(pdg)

    return known_pdgs


def diagnose_coverage(verbose: bool = True):
    """
    Run full diagnostic to check for PDG coverage gaps.

    Parameters
    ----------
    verbose : bool
        If True, print detailed information about each PDG code
    """
    print("=" * 80)
    print("PDG COVERAGE DIAGNOSTIC")
    print("=" * 80)
    print()

    # 1. Get PDG codes from each source
    print("1. Scanning simulation CSV files...")
    csv_pdgs = get_pdg_from_csvs()
    print(f"   Found {len(csv_pdgs)} unique parent PDG codes in CSV files")

    print("\n2. Querying HNLCalc production channels (m=1.0 GeV, muon coupling)...")
    hnlcalc_pdgs = get_pdg_from_hnlcalc(mass_GeV=1.0, benchmark="010")
    print(f"   Found {len(hnlcalc_pdgs)} parent PDG codes in HNLCalc")

    print("\n3. Checking cross-section lookup table...")
    xsec_pdgs = get_pdg_from_xsec_table()
    print(f"   Found {len(xsec_pdgs)} parent PDG codes with known cross-sections")

    # 2. Find gaps
    print("\n" + "=" * 80)
    print("COVERAGE ANALYSIS")
    print("=" * 80)

    # PDGs in CSV but not in HNLCalc
    csv_only = csv_pdgs - hnlcalc_pdgs
    if csv_only:
        print(f"\n⚠️  WARNING: {len(csv_only)} PDG codes in CSV but NOT in HNLCalc:")
        print(f"   {sorted(csv_only)}")
        print("   → These events will have BR=0.0 and contribute nothing to signal")
    else:
        print("\n✓ All CSV parent PDGs are covered by HNLCalc")

    # PDGs in CSV but not in cross-section table
    csv_no_xsec = csv_pdgs - xsec_pdgs
    if csv_no_xsec:
        print(f"\n⚠️  WARNING: {len(csv_no_xsec)} PDG codes in CSV but NO cross-section:")
        print(f"   {sorted(csv_no_xsec)}")
        print("   → These events will have σ=0.0 and contribute nothing to signal")
    else:
        print("\n✓ All CSV parent PDGs have cross-sections")

    # PDGs in HNLCalc but not in CSV
    hnlcalc_only = hnlcalc_pdgs - csv_pdgs
    if hnlcalc_only:
        print(f"\n⚠️  INFO: {len(hnlcalc_only)} PDG codes in HNLCalc but NOT in CSV:")
        print(f"   {sorted(hnlcalc_only)}")
        print("   → These channels exist but are not simulated")
        print("   → Consider adding them to Pythia production if important")

    # PDGs in cross-section table but not in HNLCalc
    xsec_only = xsec_pdgs - hnlcalc_pdgs
    if xsec_only:
        print(f"\n⚠️  INFO: {len(xsec_only)} PDG codes with σ but NOT in HNLCalc:")
        print(f"   {sorted(xsec_only)}")
        print("   → These mesons have cross-sections but HNLCalc doesn't model them")

    # 3. Detailed breakdown (if verbose)
    if verbose and csv_pdgs:
        print("\n" + "=" * 80)
        print("DETAILED PDG BREAKDOWN")
        print("=" * 80)
        print(f"{'PDG':<8} {'Name':<10} {'In CSV':<8} {'HNLCalc':<10} {'Has σ':<8} {'Status'}")
        print("-" * 80)

        pdg_names = {
            321: "K+",
            130: "KL0",
            310: "KS0",
            411: "D+",
            421: "D0",
            431: "Ds+",
            4122: "Λc+",
            511: "B0",
            521: "B+",
            531: "Bs0",
            541: "Bc+",
            5122: "Λb0",
            5232: "Σb",
            5332: "Ωb",
            15: "τ",
        }

        all_pdgs = csv_pdgs | hnlcalc_pdgs | xsec_pdgs

        for pdg in sorted(all_pdgs):
            name = pdg_names.get(pdg, "Unknown")
            in_csv = "✓" if pdg in csv_pdgs else "✗"
            in_hnl = "✓" if pdg in hnlcalc_pdgs else "✗"
            has_xsec = "✓" if pdg in xsec_pdgs else "✗"

            # Determine status
            if pdg in csv_pdgs and pdg in hnlcalc_pdgs and pdg in xsec_pdgs:
                status = "OK"
            elif pdg in csv_pdgs and (pdg not in hnlcalc_pdgs or pdg not in xsec_pdgs):
                status = "⚠️  LOST"
            elif pdg not in csv_pdgs and pdg in hnlcalc_pdgs:
                status = "Not simulated"
            else:
                status = "Incomplete"

            print(f"{pdg:<8} {name:<10} {in_csv:<8} {in_hnl:<10} {has_xsec:<8} {status}")

    print("\n" + "=" * 80)
    print("RECOMMENDATION")
    print("=" * 80)

    if csv_only or csv_no_xsec:
        print("⚠️  POTENTIAL SILENT DATA LOSS DETECTED")
        print()
        print("Some parent PDG codes in your simulation CSV files are missing from:")
        if csv_only:
            print("  - HNLCalc production channels (BR will be 0.0)")
        if csv_no_xsec:
            print("  - Cross-section lookup table (σ will be 0.0)")
        print()
        print("These events are being silently discarded during analysis.")
        print()
        print("Actions:")
        print("  1. Check if these PDG codes are physically relevant")
        print("  2. Add missing entries to config/production_xsecs.py")
        print("  3. Verify HNLCalc supports these production channels")
        print("  4. Consider adding diagnostic warnings in expected_signal_events()")
    else:
        print("✓ No coverage gaps detected - all CSV parents are properly handled")

    print("=" * 80)


if __name__ == "__main__":
    diagnose_coverage(verbose=True)
