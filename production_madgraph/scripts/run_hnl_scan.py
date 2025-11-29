import os
import sys
import subprocess
import shutil
import glob
import logging
import argparse
import re
import gzip
from pathlib import Path

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger("HNL_Scan")

# Constants
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
PRODUCTION_DIR = PROJECT_ROOT / "production_madgraph"
CARDS_DIR = PRODUCTION_DIR / "cards"
LHE_DIR = PRODUCTION_DIR / "lhe"
CSV_DIR = PRODUCTION_DIR / "csv"
SCRIPTS_DIR = PRODUCTION_DIR / "scripts"
WORK_DIR = PRODUCTION_DIR / "work"

SUMMARY_CSV = CSV_DIR / "summary_HNL_EW_production.csv"

# Configuration Defaults
DEFAULT_MG5_PATH = os.environ.get("MG5_PATH", "mg5_aMC") # Default assumes in PATH
DEFAULT_NEVENTS = 10000
K_FACTOR = 1.3

# Mass Grid (as per instructions)
MASS_GRID = [
    5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
    25, 28, 30, 32, 35, 38, 40, 45, 50, 55, 60, 65, 70, 75, 80
]
FLAVOURS = ["electron", "muon", "tau"]

# Flavor mappings
FLAVOUR_CONFIG = {
    "electron": {
        "l_plus": "e+", "l_minus": "e-", "vl": "ve", "vl_tilde": "ve~",
        "V1": 1.0, "V2": 0.0, "V3": 0.0
    },
    "muon": {
        "l_plus": "mu+", "l_minus": "mu-", "vl": "vm", "vl_tilde": "vm~",
        "V1": 0.0, "V2": 1.0, "V3": 0.0
    },
    "tau": {
        "l_plus": "ta+", "l_minus": "ta-", "vl": "vt", "vl_tilde": "vt~",
        "V1": 0.0, "V2": 0.0, "V3": 1.0
    }
}

def ensure_directories():
    for d in [LHE_DIR, CSV_DIR, WORK_DIR]:
        d.mkdir(parents=True, exist_ok=True)

    # Initialize summary CSV if not exists
    if not SUMMARY_CSV.exists():
        with open(SUMMARY_CSV, 'w') as f:
            f.write("mass_hnl_GeV,flavour,xsec_pb,xsec_error_pb,k_factor,n_events_generated,csv_path\n")

def get_cross_section(process_dir):
    """
    Parses the cross section from the run log or HTML.
    MadGraph usually puts it in `crossx.html` or the log file.
    The easiest is often reading `Events/run_01/run_01_tag_1_banner.txt` or similar,
    but we might just parse the stdout of the MG run if we captured it,
    or look for the `crossx.html` which contains the summary.

    However, robust way is reading the html/log.
    Let's check `Events/run_01/run_01_tag_1_banner.txt` which contains Integrated weight.
    Or `Events/run_01/unweighted_events.lhe.gz` header (it has <init> block with xsec).

    <init>
    2212 2212 7.000000e+03 7.000000e+03 0 0 247000 247000 -4 1
    1.23456e+00 4.56789e-03 1.00000e+00 1
    </init>
    Line 2: xsec(pb), error(pb), weight, process_id
    """

    # Find the LHE file in the process dir
    events_dir = Path(process_dir) / "Events"
    # Find latest run
    runs = sorted([d for d in events_dir.iterdir() if d.name.startswith("run_")], key=lambda x: x.name)
    if not runs:
        logger.error(f"No run directory found in {events_dir}")
        return 0.0, 0.0, None

    latest_run = runs[-1]
    lhe_gz = latest_run / "unweighted_events.lhe.gz"
    lhe_file = latest_run / "unweighted_events.lhe"

    target_lhe = None
    if lhe_gz.exists():
        target_lhe = lhe_gz
        open_func = gzip.open
        mode = 'rt'
    elif lhe_file.exists():
        target_lhe = lhe_file
        open_func = open
        mode = 'r'
    else:
        logger.error(f"No LHE file found in {latest_run}")
        return 0.0, 0.0, None

    xsec = 0.0
    err = 0.0

    try:
        with open_func(target_lhe, mode) as f:
            in_init = False
            for line in f:
                if '<init>' in line:
                    in_init = True
                    # Skip the beam line
                    next(f)
                    # Next line has xsec
                    xsec_line = next(f)
                    parts = xsec_line.strip().split()
                    xsec = float(parts[0])
                    err = float(parts[1])
                    break
    except Exception as e:
        logger.error(f"Error reading xsec from LHE: {e}")

    return xsec, err, target_lhe

def run_point(mg5_path, mass, flavour, nevents, dry_run=False):
    logger.info(f"--- Running Mass={mass} GeV, Flavour={flavour} ---")

    # Setup paths
    proc_name = f"hnl_{flavour}_m{mass}".replace(".", "p")
    proc_dir = WORK_DIR / proc_name

    # Cleanup previous run if exists
    if proc_dir.exists():
        shutil.rmtree(proc_dir)

    # 1. Generate Process Card
    # We substitute placeholders in ew_process_template.dat
    # Note: MadGraph needs absolute paths usually to be safe

    params = FLAVOUR_CONFIG[flavour]
    with open(CARDS_DIR / "ew_process_template.dat", 'r') as f:
        proc_card_template = f.read()

    proc_card_content = proc_card_template.format(
        l_plus=params["l_plus"],
        l_minus=params["l_minus"],
        vl=params["vl"],
        vl_tilde=params["vl_tilde"],
        output_dir=proc_dir.absolute()
    )

    proc_card_path = WORK_DIR / f"proc_card_{proc_name}.dat"
    with open(proc_card_path, 'w') as f:
        f.write(proc_card_content)

    # 2. Run MadGraph to generate process directory
    if not dry_run:
        cmd = [mg5_path, str(proc_card_path)]
        logger.info(f"Generating process: {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL) # Capture output to keep log clean?
        except subprocess.CalledProcessError as e:
            logger.error(f"MadGraph process generation failed: {e}")
            return False

    if not proc_dir.exists() and not dry_run:
        logger.error(f"Process directory {proc_dir} was not created.")
        return False

    # 3. Prepare Run and Param cards
    # Copy/Template param card
    with open(CARDS_DIR / "param_card_hnl_template.dat", 'r') as f:
        param_card_template = f.read()

    param_card_content = param_card_template.format(
        mass_n1=f"{mass:.6e}",
        V1=f"{params['V1']:.6e}",
        V2=f"{params['V2']:.6e}",
        V3=f"{params['V3']:.6e}"
    )

    # Overwrite the default param_card in the process directory
    cards_out_dir = proc_dir / "Cards"
    if not dry_run:
        with open(cards_out_dir / "param_card.dat", 'w') as f:
            f.write(param_card_content)

    # Template Run Card
    with open(CARDS_DIR / "run_card_ew_template.dat", 'r') as f:
        run_card_template = f.read()

    run_card_content = run_card_template.format(
        nevents=nevents
    )

    if not dry_run:
        with open(cards_out_dir / "run_card.dat", 'w') as f:
            f.write(run_card_content)

    # 4. Launch Run
    # We can run ./bin/generate_events from the process dir
    if not dry_run:
        bin_script = proc_dir / "bin" / "generate_events"
        # we need to pass arguments to answer prompts: 0 (no shower), 0 (no analysis) - usually handled by 'force' or input file
        # Using input file method for generate_events
        # generate_events asks:
        # 1. Run name? (Enter)
        # 2. Edit cards? (0 = done)
        # we can feed this via stdin

        # Or easier: ./bin/madevent -> launch run_01

        logger.info("Launching event generation...")
        # Input: 0 (keep defaults/done editing)
        try:
            # The 'generate_events' script usually takes interaction.
            # providing "0\n" tells it we are done editing cards.
            # If it asks for run name, it might be tricky.
            # Using 'madevent' command file is safer.

            # Create a command file for madevent
            me_cmd_path = proc_dir / "run_cmd.txt"
            with open(me_cmd_path, 'w') as f:
                f.write("launch run_01\n") # launch
                f.write("0\n")             # done editing
                f.write("set nevents {}\n".format(nevents)) # Ensure nevents is set (redundant but safe)
                f.write("0\n")             # done

            bin_madevent = proc_dir / "bin" / "madevent"
            subprocess.run([str(bin_madevent), str(me_cmd_path)], check=True, stdout=subprocess.DEVNULL)

        except subprocess.CalledProcessError as e:
            logger.error(f"Event generation failed: {e}")
            return False

    # 5. Collect Output
    xsec, xsec_err, lhe_path = get_cross_section(proc_dir)
    logger.info(f"Cross-section: {xsec} +/- {xsec_err} pb")

    if not lhe_path:
        logger.error("Could not find output LHE.")
        return False

    # Copy LHE to final destination
    dest_dir = LHE_DIR / flavour / f"m_{mass}GeV"
    dest_dir.mkdir(parents=True, exist_ok=True)
    final_lhe_path = dest_dir / f"HNL_mass_{mass}GeV_{flavour}_EW.lhe.gz"

    # Check if we need to compress
    if str(lhe_path).endswith('.lhe') and str(final_lhe_path).endswith('.gz'):
        logger.info(f"Compressing {lhe_path} to {final_lhe_path}")
        with open(lhe_path, 'rb') as f_in:
            with gzip.open(final_lhe_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copy(lhe_path, final_lhe_path)

    # 6. Convert to CSV
    csv_filename = f"HNL_mass_{mass}GeV_{flavour}_EW.csv"
    csv_path = CSV_DIR / flavour / csv_filename

    # Import conversion logic
    sys.path.append(str(SCRIPTS_DIR))
    from lhe_to_csv import process_lhe_file

    # We pass the final_lhe_path
    if not process_lhe_file(str(final_lhe_path), str(csv_path), mass):
        logger.error("LHE to CSV conversion failed.")
        return False

    # 7. Update Summary
    # Check if entry exists? No, just append.
    with open(SUMMARY_CSV, 'a') as f:
        # mass_hnl_GeV,flavour,xsec_pb,xsec_error_pb,k_factor,n_events_generated,csv_path
        # relative path for csv
        rel_csv_path = os.path.relpath(csv_path, PROJECT_ROOT)
        f.write(f"{mass},{flavour},{xsec},{xsec_err},{K_FACTOR},{nevents},{rel_csv_path}\n")

    # Cleanup work dir to save space?
    # shutil.rmtree(proc_dir) # Optional

    return True

def main():
    parser = argparse.ArgumentParser(description="Run HNL MadGraph Production Scan")
    parser.add_argument("--test", action="store_true", help="Run a small test (1 point, few events)")
    parser.add_argument("--mg5", default=DEFAULT_MG5_PATH, help="Path to mg5_aMC executable")
    args = parser.parse_args()

    ensure_directories()

    if args.test:
        logger.info("Running in TEST mode")
        grid = [20] # Single mass point
        flavs = ["electron"] # Single flavour
        nevents = 100
    else:
        grid = MASS_GRID
        flavs = FLAVOURS
        nevents = DEFAULT_NEVENTS

    # Check MG5
    if shutil.which(args.mg5) is None and not os.path.exists(args.mg5):
        logger.error(f"MadGraph executable '{args.mg5}' not found. Please set MG5_PATH or passed --mg5.")
        sys.exit(1)

    for flav in flavs:
        for mass in grid:
            success = run_point(args.mg5, mass, flav, nevents)
            if not success:
                logger.error(f"Failed at Mass={mass} Flavour={flav}")
                # Continue or exit? Continue.

    logger.info("Scan complete.")

if __name__ == "__main__":
    main()
