#!/usr/bin/env python3
"""
Mass scan script for HNL simulation
Generates PYTHIA configurations for different HNL masses and runs simulations
"""

import os
import subprocess
import sys

# Configuration
START_MASS = 15  # GeV
MASS_STEP = 8    # GeV
MAX_MASS = 79    # GeV (kinematic limit from W decay)
N_EVENTS = 1000  # Test run to verify fix

# Base configuration file
BASE_CMND = "hnlLL.cmnd"

# Output directory for mass scan
OUTPUT_DIR = "mass_scan_hnl"

def create_config_file(mass, base_file=BASE_CMND):
    """Create a modified configuration file for a specific mass"""

    # Read the base configuration
    with open(base_file, 'r') as f:
        lines = f.readlines()

    # Modify the mass line
    modified_lines = []
    for line in lines:
        if line.startswith('9900015:m0'):
            modified_lines.append(f'9900015:m0 = {mass:.1f}                    ! HNL mass = {mass} GeV\n')
        elif line.startswith('Main:numberOfEvents'):
            modified_lines.append(f'Main:numberOfEvents = {N_EVENTS}\n')
        else:
            modified_lines.append(line)

    # Create output filename in the mass scan directory
    output_file = os.path.join(OUTPUT_DIR, f"hnlLL_m{mass}GeV.cmnd")

    # Write modified configuration
    with open(output_file, 'w') as f:
        f.writelines(modified_lines)

    return output_file

def run_simulation(config_file, mass):
    """Run PYTHIA simulation with the given configuration"""

    print(f"\n{'='*70}")
    print(f"Running simulation for HNL mass = {mass} GeV")
    print(f"Config file: {config_file}")
    print(f"Events: {N_EVENTS}")
    print(f"{'='*70}\n")

    # Run PYTHIA from the mass_scan_hnl directory
    cmd = ['../main144', '-c', os.path.basename(config_file)]

    try:
        # Run in the OUTPUT_DIR so CSV files are created there
        result = subprocess.run(cmd, check=True, capture_output=False, cwd=OUTPUT_DIR)
        print(f"\n✓ Simulation completed successfully for m = {mass} GeV")
        return True
    except subprocess.CalledProcessError as e:
        print(f"\n✗ Simulation failed for m = {mass} GeV")
        print(f"Error: {e}")
        return False

def main():
    """Main execution function"""

    # Create output directory if it doesn't exist
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print(f"Created output directory: {OUTPUT_DIR}")
    else:
        print(f"Using existing output directory: {OUTPUT_DIR}")

    # Calculate mass points
    masses = []
    mass = START_MASS
    while mass <= MAX_MASS:
        masses.append(mass)
        mass += MASS_STEP

    print(f"\n{'='*70}")
    print(f"HNL MASS SCAN CONFIGURATION")
    print(f"{'='*70}")
    print(f"Mass range: {START_MASS} - {MAX_MASS} GeV")
    print(f"Mass step: {MASS_STEP} GeV")
    print(f"Mass points: {masses}")
    print(f"Events per mass: {N_EVENTS}")
    print(f"Total mass points: {len(masses)}")
    print(f"Output directory: {OUTPUT_DIR}/")
    print(f"{'='*70}\n")

    # Run simulations for each mass
    successful = []
    failed = []

    for mass in masses:
        # Create configuration file
        config_file = create_config_file(mass)
        print(f"Created config file: {config_file}")

        # Run simulation
        success = run_simulation(config_file, mass)

        if success:
            successful.append(mass)
            # Expected output CSV file
            csv_file = os.path.join(OUTPUT_DIR, f"hnlLL_m{mass}GeVLLP.csv")
            if os.path.exists(csv_file):
                print(f"✓ Output CSV: {csv_file}")
        else:
            failed.append(mass)

    # Summary
    print(f"\n{'='*70}")
    print(f"MASS SCAN SUMMARY")
    print(f"{'='*70}")
    print(f"Successful: {len(successful)}/{len(masses)} mass points")
    if successful:
        print(f"  Masses: {successful}")
    if failed:
        print(f"Failed: {len(failed)}/{len(masses)} mass points")
        print(f"  Masses: {failed}")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
