#!/usr/bin/env python3
"""
production/madgraph/run_wz_production.py

Generate HNL events via electroweak production (W/Z → ℓ N) at LHC 14 TeV.

Ported from the upstream llpatcolliders_FONLL run_wz_production.py:
  - Uses a configurable vendored MG5 (see _resolve_mg5_exe)
  - No Docker assumptions
  - Uses config_mass_grid.py mass loop
  - Parton level (no shower)
  - Output: weight,E,px,py,pz headerless CSV at
    output/llp_4vectors/{flavor}/WZ/mN_{mass}.csv

Weight: the MG5 unweighted-event weight (XWGTUP, = σ_LO / N per event) is read
directly from the LHE and then scaled by K_FACTOR_EW so the summed weight is
σ_NLO. The per-row weight is therefore σ_NLO / N in pb at U²=1, matching the
meson-channel convention.

Usage:
    python run_wz_production.py                          # full scan
    python run_wz_production.py --flavor Umu             # single flavor
    python run_wz_production.py --test                   # quick test
    python run_wz_production.py --masses 1.0 2.0 5.0
"""

import os
import sys
import subprocess
import shutil
import argparse
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config_mass_grid import MASS_GRID, N_EVENTS_DEFAULT, format_mass_for_filename
from production.constants import K_FACTOR_EW, FLAVOR_TO_MG5

PYTHON_EXE = Path(sys.executable)

# Vendored HeavyN UFO model (loaded by absolute path so MG5 picks up this copy
# regardless of which MG5 install resolves below).
MODEL_DIR = PROJECT_ROOT / "vendored" / "SM_HeavyN_CKM_AllMasses_LO"

# Directories
CARDS_DIR = Path(__file__).parent / "cards"
WORK_DIR = Path(__file__).parent / "work"
OUTPUT_BASE = PROJECT_ROOT / "output" / "llp_4vectors"

# Mixing configurations (U²=1 for the active flavor)
MIXING_CONFIGS = {
    "Ue":   {"ve1": 1.0, "vmu1": 0.0, "vtau1": 0.0},
    "Umu":  {"ve1": 0.0, "vmu1": 1.0, "vtau1": 0.0},
    "Utau": {"ve1": 0.0, "vmu1": 0.0, "vtau1": 1.0},
}


def _resolve_mg5_exe():
    """Resolve the mg5_aMC executable.

    Resolution order:
      1. $HNL_MG5_EXE (explicit path)
      2. hnl/vendored/MG5_aMC_v3_6_6/bin/mg5_aMC (drop-in vendoring)
      3. sibling llpatcolliders_FONLL/vendored/MG5_aMC_v3_6_6/bin/mg5_aMC

    The 148 MB MG5 install is intentionally not committed; see
    hnl/vendored/PROVENANCE.md.
    """
    env_path = os.environ.get("HNL_MG5_EXE")
    if env_path:
        return Path(env_path)

    vendored = PROJECT_ROOT / "vendored" / "MG5_aMC_v3_6_6" / "bin" / "mg5_aMC"
    if vendored.exists():
        return vendored

    # Sibling llpatcolliders_FONLL checkout. PROJECT_ROOT is the hnl/ dir, so
    # its parent is the repo root and its grandparent is the projects dir that
    # also holds llpatcolliders_FONLL.
    rel = Path("llpatcolliders_FONLL") / "vendored" / "MG5_aMC_v3_6_6" / "bin" / "mg5_aMC"
    for base in (PROJECT_ROOT.parent, PROJECT_ROOT.parent.parent):
        candidate = base / rel
        if candidate.exists():
            return candidate
    return PROJECT_ROOT.parent.parent / rel


MG5_EXE = _resolve_mg5_exe()


def _disable_mg5_html_opening(cards_dir):
    """Disable MG5 browser auto-open for this process directory."""
    cfg = cards_dir / "me5_configuration.txt"
    if not cfg.exists():
        return
    text = cfg.read_text()
    new_text = text.replace("# automatic_html_opening = True", "automatic_html_opening = False")
    new_text = new_text.replace("# automatic_html_opening = False", "automatic_html_opening = False")
    new_text = new_text.replace("automatic_html_opening = True", "automatic_html_opening = False")
    new_text = new_text.replace("# web_browser = None", "web_browser = None")
    if "automatic_html_opening" not in new_text:
        new_text += "\nautomatic_html_opening = False\n"
    if new_text != text:
        cfg.write_text(new_text)


def generate_process(flavor, mass):
    """Generate MadGraph process directory from proc_card template."""
    mg5_flavor = FLAVOR_TO_MG5[flavor]
    mass_label = format_mass_for_filename(mass)
    work_subdir = WORK_DIR / f"hnl_{mg5_flavor}_{mass_label}GeV"

    proc_card = CARDS_DIR / f"proc_card_{mg5_flavor}.dat"
    if not proc_card.exists():
        raise FileNotFoundError(f"Process card not found: {proc_card}")

    WORK_DIR.mkdir(parents=True, exist_ok=True)

    # Build MG5 command file
    cmd_file = WORK_DIR / f"mg5_gen_{mg5_flavor}_{mass_label}.txt"
    with open(proc_card) as f:
        proc_lines = f.readlines()

    with open(cmd_file, 'w') as f:
        # Load the vendored UFO model by absolute path.
        f.write(f"import model {MODEL_DIR}\n\n")
        f.write("set automatic_html_opening False\n")
        for line in proc_lines:
            if ('generate' in line or 'add process' in line) and not line.strip().startswith('#'):
                f.write(line)
        f.write(f"\noutput {work_subdir} -nojpeg\n")
        f.write("quit\n")

    log_file = WORK_DIR / f"mg5_gen_{mg5_flavor}_{mass_label}.log"
    with open(log_file, 'w') as log:
        result = subprocess.run(
            [str(PYTHON_EXE), str(MG5_EXE), str(cmd_file)],
            stdout=log, stderr=subprocess.STDOUT, timeout=300,
        )

    if result.returncode != 0 or not (work_subdir / 'bin' / 'generate_events').exists():
        print(f"    FAILED: process generation (see {log_file})")
        return None

    cmd_file.unlink(missing_ok=True)
    return work_subdir


def write_cards(work_subdir, flavor, mass, n_events):
    """Write run_card.dat and param_card.dat into process Cards/ directory."""
    cards_dir = work_subdir / "Cards"
    cards_dir.mkdir(exist_ok=True)

    # Run card
    run_content = (CARDS_DIR / "run_card_template.dat").read_text()
    run_content = run_content.replace("N_EVENTS_PLACEHOLDER", str(n_events))
    (cards_dir / "run_card.dat").write_text(run_content)
    _disable_mg5_html_opening(cards_dir)

    # Param card
    param_content = (CARDS_DIR / "param_card_template.dat").read_text()
    mixing = MIXING_CONFIGS[flavor]
    param_content = param_content.replace("MASS_N1_PLACEHOLDER", f"{mass:.6e}")
    param_content = param_content.replace("VE1_PLACEHOLDER", f"{mixing['ve1']:.6e}")
    param_content = param_content.replace("VMU1_PLACEHOLDER", f"{mixing['vmu1']:.6e}")
    param_content = param_content.replace("VTAU1_PLACEHOLDER", f"{mixing['vtau1']:.6e}")
    (cards_dir / "param_card.dat").write_text(param_content)


def run_events(work_subdir, nb_core=1):
    """Run MadGraph generate_events → return LHE file path or None."""
    log_file = work_subdir / "generate_events.log"
    cmd = [
        str(PYTHON_EXE), "bin/generate_events",
        "-f", "--laststep=parton",
    ]
    if nb_core > 1:
        cmd += ["--multicore", f"--nb_core={nb_core}"]
    else:
        cmd += ["--nb_core=1"]

    with open(log_file, 'w') as log:
        result = subprocess.run(
            cmd, stdout=log, stderr=subprocess.STDOUT,
            cwd=work_subdir, timeout=3600,
        )

    if result.returncode != 0:
        print(f"    FAILED: event generation (see {log_file})")
        return None

    events_dir = work_subdir / "Events"
    if events_dir.exists():
        for run_dir in sorted(events_dir.glob("run_*")):
            for lhe in [run_dir / "unweighted_events.lhe.gz",
                        run_dir / "unweighted_events.lhe"]:
                if lhe.exists():
                    return lhe
    return None


def convert_lhe(lhe_path, csv_path):
    """Convert LHE → CSV using our simplified parser."""
    from production.madgraph.lhe_to_csv import LHEParser
    parser = LHEParser(lhe_path)
    return parser.write_hnl_csv(csv_path)


def run_single_point(flavor, mass, n_events, nb_core=1):
    """Full pipeline for one (flavor, mass) point."""
    mass_label = format_mass_for_filename(mass)
    csv_dir = OUTPUT_BASE / flavor / "WZ"
    csv_dir.mkdir(parents=True, exist_ok=True)
    csv_path = csv_dir / f"mN_{mass_label}.csv"

    print(f"\n  [{flavor}] m_N = {mass} GeV")

    # Step 1: Generate process
    work_subdir = generate_process(flavor, mass)
    if work_subdir is None:
        return False

    # Step 2: Write cards
    write_cards(work_subdir, flavor, mass, n_events)

    # Step 3: Generate events
    lhe_path = run_events(work_subdir, nb_core=nb_core)
    if lhe_path is None:
        return False

    # Step 4: Convert LHE → CSV
    n_ev = convert_lhe(lhe_path, csv_path)
    if n_ev is None or n_ev == 0:
        print(f"    FAILED: no HNL events extracted")
        return False

    print(f"    OK: {n_ev} events → {csv_path}")

    # Step 5: Apply the EW K-factor by rescaling weights.
    # The parser writes the LHE per-event weight (XWGTUP = σ_LO / N) into
    # column 0; multiplying by K_FACTOR_EW promotes the summed weight from
    # σ_LO to σ_NLO without changing the kinematics.
    import numpy as np
    try:
        data = np.loadtxt(csv_path, delimiter=",")
    except ValueError:
        data = np.empty((0, 5))
    if data.ndim == 1 and data.size > 0:
        data = data.reshape(1, -1)
    if data.ndim == 2 and len(data) > 0:
        data[:, 0] *= K_FACTOR_EW
        np.savetxt(csv_path, data, delimiter=",", fmt="%.8e")

    # Cleanup
    shutil.rmtree(work_subdir, ignore_errors=True)
    return True


def main():
    parser = argparse.ArgumentParser(description="W/Z → ℓ N production via MadGraph")
    parser.add_argument("--flavor", choices=["Ue", "Umu", "Utau"], nargs="+",
                        default=["Ue", "Umu", "Utau"])
    parser.add_argument("--masses", type=float, nargs="+", default=None)
    parser.add_argument("--nevents", type=int, default=N_EVENTS_DEFAULT)
    parser.add_argument("--nb-core", type=int, default=1,
                        help="Number of CPU cores per MG5 generate_events job")
    parser.add_argument("--test", action="store_true",
                        help="Quick test: single point (Umu, 1.0 GeV, 1000 events)")
    parser.add_argument("--min-mass", type=float, default=None)
    args = parser.parse_args()

    if args.test:
        flavors = ["Umu"]
        masses = [1.0]
        n_events = 1000
    else:
        flavors = args.flavor
        masses = args.masses if args.masses else MASS_GRID
        n_events = args.nevents

    nb_core = args.nb_core
    if nb_core < 1:
        print("ERROR: --nb-core must be >= 1")
        return 1

    if args.min_mass is not None:
        masses = [m for m in masses if m >= args.min_mass]

    # Verify MG5 exists
    if not MG5_EXE.exists():
        print(f"ERROR: MadGraph not found at {MG5_EXE}")
        print("Set $HNL_MG5_EXE or vendor MG5 under hnl/vendored/ "
              "(see hnl/vendored/PROVENANCE.md).")
        return 1
    if not MODEL_DIR.exists():
        print(f"ERROR: HeavyN UFO model not found at {MODEL_DIR}")
        return 1

    print(f"W/Z → ℓ N Production")
    print(f"  MG5: {MG5_EXE}")
    print(f"  Model: {MODEL_DIR}")
    print(f"  Flavors: {flavors}")
    print(f"  Masses: {len(masses)} points")
    print(f"  Events/point: {n_events}")
    print(f"  Cores/point: {nb_core}")

    n_ok = 0
    n_fail = 0
    for flavor in flavors:
        for mass in masses:
            if run_single_point(flavor, mass, n_events, nb_core=nb_core):
                n_ok += 1
            else:
                n_fail += 1

    print(f"\nDone: {n_ok} OK, {n_fail} failed")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
