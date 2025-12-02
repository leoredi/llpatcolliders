#!/usr/bin/env python3
"""
MadGraph HNL Production Driver (Docker-friendly)

Generates HNL events via electroweak production (W/Z → ℓ N) at LHC 14 TeV.

This script:
1. Generates MadGraph process directories from proc_card_*.dat templates
2. Writes run_card.dat and param_card.dat into each process Cards/ directory
3. Runs MadGraph's generate_events to produce LHE files
4. Converts LHE → CSV using lhe_to_csv.LHEParser
5. Extracts cross-sections and appends a summary CSV

Usage (inside Docker container, with repo mounted at /work):

    # From /work/production/madgraph:
    python3 scripts/run_hnl_scan.py               # full scan
    python3 scripts/run_hnl_scan.py --test        # quick test (muon, 15 GeV)
    python3 scripts/run_hnl_scan.py --flavour muon
    python3 scripts/run_hnl_scan.py --masses 10 15 20 --flavour electron
    python3 scripts/run_hnl_scan.py --nevents 10000

Environment:

    MG5_PATH  (optional)  Absolute path to mg5_aMC binary.
                          Default inside Docker: /opt/MG5_aMC_v3_6_6/bin/mg5_aMC
"""

import sys
import os
import subprocess
import shutil
import re
import argparse
from pathlib import Path
from datetime import datetime

# Import mass grid from central configuration
# Add project root to path to import config_mass_grid
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent.parent  # production/madgraph_production/scripts/ -> repo root
sys.path.insert(0, str(PROJECT_ROOT))

try:
    from config_mass_grid import (
        ELECTRON_MASSES_EW,
        MUON_MASSES_EW,
        TAU_MASSES_EW,
    )
except ImportError as e:
    print(f"ERROR: Could not import mass grid from config_mass_grid.py")
    print(f"  Project root: {PROJECT_ROOT}")
    print(f"  Error: {e}")
    sys.exit(1)


# ============================================================================
# CONFIGURATION
# ============================================================================

# Default MG5 path – matches your Dockerfile setup
DEFAULT_MG5_EXE = os.environ.get("MG5_PATH", "/opt/MG5_aMC_v3_6_6/bin/mg5_aMC")

# Mass grid for EW regime (5-80 GeV) - loaded from config_mass_grid.py
# Combine all flavours for backward compatibility (when --masses not specified)
MASS_GRID_FULL = sorted(set(ELECTRON_MASSES_EW + MUON_MASSES_EW + TAU_MASSES_EW))

# Flavours
FLAVOURS = ['electron', 'muon', 'tau']

# Number of events per mass point
N_EVENTS_DEFAULT = 50000

# K-factor for NLO correction (approximate)
K_FACTOR = 1.3

# Mixing parameter configurations
# At generation time, we set |U_ℓ|² = 1 to get σ_max
# Analysis will rescale by actual |U_ℓ|²
MIXING_CONFIGS = {
    'electron': {'ve1': 1.0, 'vmu1': 0.0, 'vtau1': 0.0},
    'muon':    {'ve1': 0.0, 'vmu1': 1.0, 'vtau1': 0.0},
    'tau':     {'ve1': 0.0, 'vmu1': 0.0, 'vtau1': 1.0},
}


# ============================================================================
# DIRECTORY STRUCTURE
# ============================================================================

class ProjectPaths:
    """Manage project directory structure"""

    def __init__(self, base_dir=None):
        """
        Initialize paths relative to production/madgraph_production/

        Args:
            base_dir: Base directory (default: auto-detect from script location)
        """
        if base_dir is None:
            # Script is in production/madgraph_production/scripts/
            script_dir = Path(__file__).parent
            base_dir = script_dir.parent

        self.base = Path(base_dir).resolve()
        self.cards = self.base / 'cards'
        self.scripts = self.base / 'scripts'
        self.lhe_dir = self.base / 'lhe'

        # CSV output directory: use central output/csv/simulation_new/
        # This matches the Pythia production output location
        self.csv_dir = self.base.parent.parent / 'output' / 'csv' / 'simulation_new'

        self.work_dir = self.base / 'work'
        self.mg5_dir = self.base / 'mg5'  # unused in Docker but kept for completeness

        # MadGraph executable
        if Path(DEFAULT_MG5_EXE).is_absolute():
            self.mg5_exe = Path(DEFAULT_MG5_EXE)
        else:
            # Treat as relative to base (for non-Docker users)
            self.mg5_exe = self.base / DEFAULT_MG5_EXE

    def lhe_path(self, flavour, mass):
        """Get LHE output directory for (flavour, mass) (not strictly used)"""
        return self.lhe_dir / flavour / f"m_{mass}GeV"

    def csv_path(self, flavour, mass):
        """Get CSV output path for (flavour, mass) - matches Pythia format"""
        # Use same directory structure as Pythia production for compatibility
        # Format: HNL_{mass}GeV_{flavour}_ew.csv
        mass_str = f"{mass:.1f}".replace('.', 'p')  # e.g., 15.0 → 15p0
        return self.csv_dir / f"HNL_{mass_str}GeV_{flavour}_ew.csv"

    def summary_csv_path(self):
        """Get summary CSV path"""
        return self.csv_dir / "summary_HNL_ew_production.csv"

    def work_subdir(self, flavour, mass):
        """Get working directory for this run"""
        return self.work_dir / f"hnl_{flavour}_{mass}GeV"


# ============================================================================
# CARD GENERATION & PROCESS SETUP
# ============================================================================

def generate_process(paths, flavour, mass):
    """
    Step 1: Generate MadGraph process directory

    Creates the process directory with bin/, Cards/, SubProcesses/, etc.
    Does NOT yet fill in the Cards/ - that comes next.

    Args:
        paths: ProjectPaths object
        flavour: Lepton flavour
        mass: HNL mass in GeV

    Returns:
        Path: work directory
    """
    print(f"  Step 1: Generating process for {flavour}, m = {mass} GeV")

    # Read process card template
    proc_template_path = paths.cards / f"proc_card_{flavour}.dat"
    if not proc_template_path.exists():
        raise FileNotFoundError(f"Process card not found: {proc_template_path}")

    # Create temporary command file
    work_dir = paths.work_subdir(flavour, mass)
    work_dir.parent.mkdir(parents=True, exist_ok=True)

    mg5_cmd_file = work_dir.parent / f"mg5_gen_{flavour}_{mass}.txt"

    # Read process definition from template
    with open(proc_template_path) as f:
        proc_lines = f.readlines()

    # Write MadGraph command file
    with open(mg5_cmd_file, 'w') as f:
        f.write("# MadGraph process generation\n")
        f.write("import model SM_HeavyN_CKM_AllMasses_LO\n\n")

        # Copy process definitions (generate/add process lines)
        for line in proc_lines:
            if ('generate' in line or 'add process' in line) and not line.strip().startswith('#'):
                f.write(line)

        # Output to work directory
        f.write(f"\noutput {work_dir} -nojpeg\n")
        f.write("quit\n")

    print(f"    → Running MadGraph to create process directory...")

    # Run MadGraph directly (Docker: mg5_aMC is installed and in path /opt)
    try:
        cmd = [
            str(paths.mg5_exe),
            str(mg5_cmd_file)
        ]

        log_file = work_dir.parent / f"mg5_gen_{flavour}_{mass}.log"
        with open(log_file, 'w') as log:
            result = subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=300
            )

        if result.returncode != 0:
            print(f"    ✗ Process generation failed")
            print(f"      See: {log_file}")
            return None

        # Verify process directory was created
        if not (work_dir / 'bin' / 'generate_events').exists():
            print(f"    ✗ Process directory incomplete")
            return None

        print(f"    ✓ Process directory created: {work_dir}")

        # Clean up temp file
        try:
            mg5_cmd_file.unlink()
        except OSError:
            pass

        return work_dir

    except subprocess.TimeoutExpired:
        print(f"    ✗ Process generation timed out")
        return None
    except Exception as e:
        print(f"    ✗ Error: {e}")
        return None


def write_cards_to_process(paths, work_dir, flavour, mass, n_events):
    """
    Step 2: Write run_card and param_card into the process Cards/ directory

    Args:
        paths: ProjectPaths object
        work_dir: Process directory from generate_process()
        flavour: Lepton flavour
        mass: HNL mass in GeV
        n_events: Number of events

    Returns:
        bool: True if successful
    """
    print(f"  Step 2: Writing cards to {work_dir}/Cards/")

    cards_dir = work_dir / 'Cards'
    cards_dir.mkdir(exist_ok=True)

    # Read templates
    run_template_path = paths.cards / "run_card_template.dat"
    param_template_path = paths.cards / "param_card_template.dat"

    if not run_template_path.exists():
        raise FileNotFoundError(f"Run card template not found: {run_template_path}")
    if not param_template_path.exists():
        raise FileNotFoundError(f"Param card template not found: {param_template_path}")

    # 1. Run card
    run_content = run_template_path.read_text()
    run_content = run_content.replace('N_EVENTS_PLACEHOLDER', str(n_events))
    (cards_dir / 'run_card.dat').write_text(run_content)

    # 2. Param card (mass + mixing)
    param_content = param_template_path.read_text()

    mixing = MIXING_CONFIGS[flavour]
    param_content = param_content.replace('MASS_N1_PLACEHOLDER', f'{mass:.6e}')
    param_content = param_content.replace('VE1_PLACEHOLDER', f'{mixing["ve1"]:.6e}')
    param_content = param_content.replace('VMU1_PLACEHOLDER', f'{mixing["vmu1"]:.6e}')
    param_content = param_content.replace('VTAU1_PLACEHOLDER', f'{mixing["vtau1"]:.6e}')

    (cards_dir / 'param_card.dat').write_text(param_content)

    print(f"    ✓ Cards written to {cards_dir}")

    return True


# ============================================================================
# EVENT GENERATION
# ============================================================================

def generate_events(paths, work_dir):
    """
    Step 3: Run event generation using the process directory's generate_events script

    Args:
        paths: ProjectPaths object
        work_dir: Process directory from generate_process()

    Returns:
        Path: LHE file path if successful, None otherwise
    """
    print(f"  Step 3: Generating events...")

    # Verify generate_events script exists
    generate_events_script = work_dir / 'bin' / 'generate_events'

    if not generate_events_script.exists():
        print(f"    ✗ generate_events script not found: {generate_events_script}")
        return None

    try:
        # In the Docker image we have python3 in PATH
        cmd = [
            'python3',
            'bin/generate_events',
            '-f',                  # Force (non-interactive)
            '--laststep=parton',   # Stop at parton level (no shower)
        ]

        log_file = work_dir / 'generate_events.log'

        print(f"    → Running: {' '.join(cmd)}")
        print(f"    → Working directory: {work_dir}")
        print(f"    → Log: {log_file}")

        env = os.environ.copy()

        with open(log_file, 'w') as log:
            result = subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                cwd=work_dir,
                env=env,
                timeout=3600  # 1 hour timeout
            )

        if result.returncode != 0:
            print(f"    ✗ Event generation failed (exit code {result.returncode})")
            print(f"      See: {log_file}")
            return None

        # Find LHE file
        events_dir = work_dir / 'Events'
        if not events_dir.exists():
            print(f"    ✗ Events directory not created")
            return None

        # Look for first run_*/unweighted_events.lhe(.gz)
        for run_dir in sorted(events_dir.glob("run_*")):
            for lhe_file in [run_dir / "unweighted_events.lhe.gz",
                             run_dir / "unweighted_events.lhe"]:
                if lhe_file.exists():
                    print(f"    ✓ Events generated: {lhe_file}")
                    return lhe_file

        print(f"    ✗ No LHE file found in {events_dir}")
        return None

    except subprocess.TimeoutExpired:
        print(f"    ✗ Event generation timed out")
        return None
    except Exception as e:
        print(f"    ✗ Error: {e}")
        return None


# ============================================================================
# CROSS-SECTION EXTRACTION
# ============================================================================

def extract_cross_section(work_dir):
    """
    Extract cross-section from MadGraph output

    Args:
        work_dir: Working directory with MadGraph output

    Returns:
        dict: {'xsec_pb': float, 'xsec_error_pb': float, 'n_events': int}
              or None if extraction fails
    """
    # Look for cross-section in banner file or log
    banner_paths = list(Path(work_dir).rglob('*_banner.txt'))
    log_path = work_dir / 'madgraph.log'  # may not exist, that's fine

    xsec_pb = None
    xsec_error_pb = None
    n_events = None

    # Try banner file first
    if banner_paths:
        banner_path = banner_paths[0]
        with open(banner_path) as f:
            content = f.read()

        # Pattern 1: "Integrated weight (pb)  :  1.234e+02"
        match = re.search(r'Integrated weight.*:\s*([\d.eE+-]+)', content)
        if match:
            xsec_pb = float(match.group(1))

        # Pattern 2: Look for error
        match_err = re.search(r'error.*:\s*([\d.eE+-]+)', content, re.IGNORECASE)
        if match_err:
            xsec_error_pb = float(match_err.group(1))

        # Pattern 3: Number of events
        match_n = re.search(r'Number of Events.*:\s*(\d+)', content)
        if match_n:
            n_events = int(match_n.group(1))

    # Fallback: try log file
    if xsec_pb is None and log_path.exists():
        with open(log_path) as f:
            content = f.read()

        match = re.search(r'Cross-section.*:\s*([\d.eE+-]+)\s*\+/-\s*([\d.eE+-]+)', content)
        if match:
            xsec_pb = float(match.group(1))
            xsec_error_pb = float(match.group(2))

    if xsec_pb is None:
        print(f"    Warning: Could not extract cross-section from {work_dir}")
        return None

    return {
        'xsec_pb': xsec_pb,
        'xsec_error_pb': xsec_error_pb or 0.0,
        'n_events': n_events or 0
    }


# ============================================================================
# LHE → CSV CONVERSION
# ============================================================================

def convert_lhe_to_csv(paths, flavour, mass, lhe_file):
    """
    Convert LHE output to CSV format

    Args:
        paths: ProjectPaths object
        flavour: Lepton flavour
        mass: HNL mass in GeV
        lhe_file: Path to LHE file

    Returns:
        int: Number of events converted, or None if failed
    """
    print(f"  Step 4: Converting LHE → CSV...")
    print(f"    → LHE file: {lhe_file}")

    # Output CSV path
    csv_output = paths.csv_path(flavour, mass)
    csv_output.parent.mkdir(parents=True, exist_ok=True)

    # Import and run LHE parser
    sys.path.insert(0, str(paths.scripts))
    from lhe_to_csv import LHEParser

    try:
        parser = LHEParser(lhe_file, mass, flavour)
        n_events = parser.write_csv(csv_output)
        print(f"    ✓ Converted {n_events} events")
        print(f"    ✓ CSV: {csv_output}")
        return n_events

    except Exception as e:
        print(f"    ✗ LHE conversion failed: {e}")
        return None


# ============================================================================
# SUMMARY CSV
# ============================================================================

def initialize_summary_csv(paths):
    """Create summary CSV with header if it doesn't exist"""
    summary_path = paths.summary_csv_path()
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    if not summary_path.exists():
        header = "mass_hnl_GeV,flavour,xsec_pb,xsec_error_pb,k_factor,n_events_generated,csv_path,timestamp"
        with open(summary_path, 'w') as f:
            f.write(header + '\n')


def append_to_summary(paths, flavour, mass, xsec_data, n_events_csv):
    """Append result to summary CSV"""
    summary_path = paths.summary_csv_path()
    # Use just filename since all CSVs are in same directory
    csv_rel_path = paths.csv_path(flavour, mass).name
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    row = (
        f"{mass:.1f},"
        f"{flavour},"
        f"{xsec_data['xsec_pb']:.6e},"
        f"{xsec_data['xsec_error_pb']:.6e},"
        f"{K_FACTOR:.2f},"
        f"{n_events_csv},"
        f"{csv_rel_path},"
        f"{timestamp}"
    )

    with open(summary_path, 'a') as f:
        f.write(row + '\n')


# ============================================================================
# MAIN PIPELINE
# ============================================================================

def run_single_point(paths, flavour, mass, n_events):
    """
    Run pipeline for single (flavour, mass) point

    Three-step workflow:
    1. Generate MadGraph process directory
    2. Write run_card and param_card into Cards/
    3. Run generate_events to produce LHE
    4. Extract cross-section
    5. Convert LHE → CSV
    6. Append summary row

    Args:
        paths: ProjectPaths object
        flavour: Lepton flavour
        mass: HNL mass in GeV
        n_events: Number of events to generate

    Returns:
        bool: True if successful
    """
    print(f"\n{'='*70}")
    print(f"Processing: {flavour} coupling, m_HNL = {mass} GeV")
    print(f"{'='*70}")

    try:
        # Step 1: Generate process directory
        work_dir = generate_process(paths, flavour, mass)
        if work_dir is None:
            return False

        # Step 2: Write cards into Cards/ directory
        success = write_cards_to_process(paths, work_dir, flavour, mass, n_events)
        if not success:
            return False

        # Step 3: Run event generation
        lhe_file = generate_events(paths, work_dir)
        if lhe_file is None:
            return False

        # Step 4: Extract cross-section
        xsec_data = extract_cross_section(work_dir)
        if xsec_data is None:
            print("  ✗ Failed to extract cross-section")
            return False

        print(f"  Cross-section: {xsec_data['xsec_pb']:.3e} ± {xsec_data['xsec_error_pb']:.3e} pb")

        # Step 5: Convert LHE → CSV
        n_events_csv = convert_lhe_to_csv(paths, flavour, mass, lhe_file)
        if n_events_csv is None:
            return False

        # Step 6: Update summary
        append_to_summary(paths, flavour, mass, xsec_data, n_events_csv)
        print(f"  ✓ Added to summary CSV")

        return True

    except Exception as e:
        print(f"  ✗ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description='MadGraph HNL Production Driver'
    )
    parser.add_argument(
        '--test',
        action='store_true',
        help='Test mode: single mass point (15 GeV muon, 1000 events)'
    )
    parser.add_argument(
        '--flavour',
        choices=FLAVOURS,
        help='Run single flavour only'
    )
    parser.add_argument(
        '--masses',
        type=float,
        nargs='+',
        help='Custom mass list (GeV)'
    )
    parser.add_argument(
        '--nevents',
        type=int,
        default=N_EVENTS_DEFAULT,
        help=f'Number of events per point (default: {N_EVENTS_DEFAULT})'
    )

    args = parser.parse_args()

    # Setup paths
    paths = ProjectPaths()

    print("="*70)
    print("MadGraph HNL Production Pipeline")
    print("="*70)
    print(f"Base directory: {paths.base}")
    print(f"MadGraph:       {paths.mg5_exe}")
    print()

    # Verify MadGraph exists (only if absolute path; if not, we trust PATH)
    if paths.mg5_exe.is_absolute() and not paths.mg5_exe.exists():
        print(f"ERROR: MadGraph not found at {paths.mg5_exe}")
        print(f"Set MG5_PATH env var or adjust DEFAULT_MG5_EXE in run_hnl_scan.py")
        return 1

    # Determine run configuration
    if args.test:
        masses = [15.0]
        flavours = ['muon']
        n_events = 1000
        print("TEST MODE: Single point (15 GeV muon, 1000 events)")
    else:
        masses = args.masses if args.masses else MASS_GRID_FULL
        flavours = [args.flavour] if args.flavour else FLAVOURS
        n_events = args.nevents

    print(f"Masses:          {masses}")
    print(f"Flavours:        {flavours}")
    print(f"Events per point:{n_events}")
    print()

    # Initialize summary CSV
    initialize_summary_csv(paths)

    # Run pipeline
    n_total = len(masses) * len(flavours)
    n_success = 0
    n_failed = 0

    for flavour in flavours:
        for mass in masses:
            success = run_single_point(paths, flavour, mass, n_events)
            if success:
                n_success += 1
            else:
                n_failed += 1

    # Final summary
    print("\n" + "="*70)
    print("PIPELINE COMPLETE")
    print("="*70)
    print(f"Total points: {n_total}")
    print(f"Successful:   {n_success}")
    print(f"Failed:       {n_failed}")
    print(f"\nSummary CSV: {paths.summary_csv_path()}")
    print(f"Event CSVs:   {paths.csv_dir}/")
    print("="*70)

    return 0 if n_failed == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
