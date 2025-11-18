#!/usr/bin/env python3
"""
Mass scan script for HNL simulation
Generates PYTHIA configurations for different HNL masses and runs simulations
Supports both muon-coupled (mu) and tau-coupled (tau) scenarios
"""

import os
import subprocess
import sys
import argparse
from multiprocessing import Pool
from datetime import datetime

# Configuration
START_MASS = 15  # GeV
MASS_STEP = 8    # GeV
MAX_MASS = 71    # GeV (exclude 79 GeV point)
N_EVENTS = 1_000_000  # 1 million events per mass point

# Scenario configurations
SCENARIO_CONFIG = {
    'mu': {
        'base_cmnd': 'hnlLL.cmnd',
        'output_prefix': 'hnlLL',
        'description': 'muon-coupled HNL (W → μ N)'
    },
    'tau': {
        'base_cmnd': 'hnlTauLL.cmnd',
        'output_prefix': 'hnlTauLL',
        'description': 'tau-coupled HNL (W → τ N)'
    }
}

def create_config_file(mass, base_file, output_prefix):
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

    # Create output filename in the pythiaStuff directory
    output_file = f"{output_prefix}_m{mass}GeV.cmnd"

    # Write modified configuration
    with open(output_file, 'w') as f:
        f.writelines(modified_lines)

    return output_file

def run_simulation(args):
    """Run PYTHIA simulation with the given configuration"""
    mass, scenario_cfg = args

    base_file = scenario_cfg['base_cmnd']
    output_prefix = scenario_cfg['output_prefix']

    # Create config file
    config_file = create_config_file(mass, base_file, output_prefix)

    start_time = datetime.now()
    print(f"[{start_time.strftime('%H:%M:%S')}] STARTING: {scenario_cfg['description']} mass = {mass} GeV")
    print(f"  Config: {config_file}")
    print(f"  Events: {N_EVENTS:,}")
    sys.stdout.flush()

    # Run PYTHIA from current directory (pythiaStuff/)
    cmd = ['./main144', '-c', config_file]

    try:
        # Run from current directory (pythiaStuff/)
        # CSV files will be created in ../output/csv/ by main144
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()

        print(f"[{end_time.strftime('%H:%M:%S')}] COMPLETED: {scenario_cfg['description']} mass = {mass} GeV ({duration:.1f}s)")

        # Check for output CSV
        csv_file = f"../output/csv/{output_prefix}_m{mass}GeVLLP.csv"
        if os.path.exists(csv_file):
            print(f"  Output: ../output/csv/{output_prefix}_m{mass}GeVLLP.csv")

        sys.stdout.flush()
        return (mass, True, duration)
    except subprocess.CalledProcessError as e:
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        print(f"[{end_time.strftime('%H:%M:%S')}] FAILED: {scenario_cfg['description']} mass = {mass} GeV ({duration:.1f}s)")
        print(f"  Error: {e}")
        sys.stdout.flush()
        return (mass, False, duration)

def main():
    """Main execution function"""

    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Run HNL mass scan for muon-coupled or tau-coupled scenarios',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run muon-coupled scan (default)
  python run_mass_scan.py
  python run_mass_scan.py --scenario mu

  # Run tau-coupled scan
  python run_mass_scan.py --scenario tau
        """
    )
    parser.add_argument(
        '--scenario',
        type=str,
        choices=['mu', 'tau'],
        default='mu',
        help='HNL coupling scenario: mu (W→μN) or tau (W→τN). Default: mu'
    )
    args = parser.parse_args()

    # Get scenario configuration
    scenario_cfg = SCENARIO_CONFIG[args.scenario]

    # Calculate mass points
    masses = []
    mass = START_MASS
    while mass <= MAX_MASS:
        masses.append(mass)
        mass += MASS_STEP

    print(f"\n{'='*70}")
    print(f"HNL MASS SCAN CONFIGURATION")
    print(f"{'='*70}")
    print(f"Scenario: {scenario_cfg['description']}")
    print(f"Base config: {scenario_cfg['base_cmnd']}")
    print(f"Output prefix: {scenario_cfg['output_prefix']}")
    print(f"Mass range: {START_MASS} - {MAX_MASS} GeV")
    print(f"Mass step: {MASS_STEP} GeV")
    print(f"Mass points: {masses}")
    print(f"Events per mass: {N_EVENTS:,}")
    print(f"Total mass points: {len(masses)}")
    print(f"Parallel processes: 2")
    print(f"CSV output directory: ../output/csv/")
    print(f"{'='*70}\n")

    # Run simulations in parallel (2 at a time)
    scan_start = datetime.now()
    print(f"Scan started at {scan_start.strftime('%H:%M:%S')}\n")

    # Create argument tuples for run_simulation
    sim_args = [(mass, scenario_cfg) for mass in masses]

    with Pool(processes=2) as pool:
        results = pool.map(run_simulation, sim_args)

    scan_end = datetime.now()
    total_duration = (scan_end - scan_start).total_seconds()

    # Process results
    successful = []
    failed = []
    total_sim_time = 0

    for mass, success, duration in results:
        total_sim_time += duration
        if success:
            successful.append(mass)
        else:
            failed.append(mass)

    # Summary
    print(f"\n{'='*70}")
    print(f"MASS SCAN SUMMARY")
    print(f"{'='*70}")
    print(f"Scan completed at {scan_end.strftime('%H:%M:%S')}")
    print(f"Wall-clock time: {total_duration/60:.1f} minutes")
    print(f"Total simulation time: {total_sim_time/60:.1f} minutes")
    print(f"Speedup from parallelization: {total_sim_time/total_duration:.2f}x")
    print(f"\nSuccessful: {len(successful)}/{len(masses)} mass points")
    if successful:
        print(f"  Masses: {successful} GeV")
    if failed:
        print(f"Failed: {len(failed)}/{len(masses)} mass points")
        print(f"  Masses: {failed} GeV")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
