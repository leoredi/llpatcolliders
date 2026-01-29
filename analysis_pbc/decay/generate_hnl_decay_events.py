#!/usr/bin/env python3
"""
Generate HNL decay events (5-20 GeV) with MadGraph5 + Pythia8 hadronization.

Outputs stable final-state particles in the HNL rest frame for displaced-vertex
studies in transverse detectors.
"""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List

import numpy as np
import pandas as pd

import sys

THIS_FILE = Path(__file__).resolve()
ANALYSIS_ROOT = THIS_FILE.parents[1]
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from models.hnl_model_hnlcalc import HNLModel

HBARC_GEV_M = 1.973269804e-16  # GeV * m


@dataclass
class RunConfig:
    masses: List[float]
    nevents: int
    Ue2: float
    Umu2: float
    Utau2: float
    mg5_path: Path
    work_dir: Path
    out_dir: Path
    seed: int


def hnl_ctau_and_width(mass_GeV: float, Ue2: float, Umu2: float, Utau2: float) -> tuple[float, float]:
    """
    Return (ctau0_m, gamma_tot_GeV) using HNLCalc (per UCI-TR-2024-01 formulas).
    """
    model = HNLModel(mass_GeV=mass_GeV, Ue2=Ue2, Umu2=Umu2, Utau2=Utau2)
    ctau0_m = float(model.ctau0_m)
    if ctau0_m <= 0.0:
        raise ValueError(f"Invalid ctau0 for mass={mass_GeV} GeV")
    gamma_tot = HBARC_GEV_M / ctau0_m
    return ctau0_m, gamma_tot


def write_param_card(template_path: Path, out_path: Path, mass_GeV: float, Ue2: float, Umu2: float, Utau2: float) -> None:
    text = template_path.read_text()
    replacements = {
        "MASS_N1_PLACEHOLDER": f"{mass_GeV:.6e}",
        "VE1_PLACEHOLDER": f"{np.sqrt(Ue2):.6e}",
        "VMU1_PLACEHOLDER": f"{np.sqrt(Umu2):.6e}",
        "VTAU1_PLACEHOLDER": f"{np.sqrt(Utau2):.6e}",
    }
    for key, value in replacements.items():
        text = text.replace(key, value)
    out_path.write_text(text)


def build_mg5_process_command(proc_dir: Path) -> str:
    return f"""
import model SM_HeavyN_CKM_AllMasses_LO

define l- = e- mu- ta-
define l+ = e+ mu+ ta+
define vl = ve vm vt
define vl~ = ve~ vm~ vt~
define utype = u c
define dtype = d s b
define utype~ = u~ c~
define dtype~ = d~ s~ b~
define q = u c d s b
define q~ = u~ c~ d~ s~ b~

generate n1 > l- utype dtype~
add process n1 > l+ utype~ dtype
add process n1 > vl q q~
add process n1 > vl~ q q~
add process n1 > l- l+ vl
add process n1 > l- l+ vl~

output {proc_dir} -nojpeg
"""


def build_mg5_launch_command(proc_dir: Path, nevents: int, seed: int) -> str:
    return f"""
launch {proc_dir}
set nevents {nevents}
set iseed {seed}
set lpp1 0
set lpp2 0
set ebeam1 0
set ebeam2 0
set pdlabel none
set lhaid 0
"""


def run_madgraph(mg5_path: Path, cmd_text: str, work_dir: Path) -> None:
    cmd_file = work_dir / "mg5_decay_cmd.dat"
    cmd_file.write_text(cmd_text.strip() + "\n")
    subprocess.run([str(mg5_path), str(cmd_file)], check=True, cwd=work_dir)


def find_lhe_file(proc_dir: Path) -> Path:
    for candidate in proc_dir.glob("Events/*/unweighted_events.lhe*"):
        return candidate
    raise FileNotFoundError(f"No LHE file found under {proc_dir}")


def ensure_lhe_uncompressed(lhe_path: Path, work_dir: Path) -> Path:
    if lhe_path.suffix != ".gz":
        return lhe_path
    out_path = work_dir / "unweighted_events.lhe"
    with gzip.open(lhe_path, "rb") as src, out_path.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    return out_path


def load_pythia() -> "pythia8.Pythia":
    _ensure_pythia8_pythonpath()
    try:
        import pythia8
    except ImportError as exc:
        raise ImportError(
            "Pythia8 Python bindings not found. "
            "Set PYTHONPATH to the directory containing pythia8*.so, "
            "or install a build with bindings enabled."
        ) from exc
    return pythia8.Pythia()


def _ensure_pythia8_pythonpath() -> None:
    try:
        import pythia8  # noqa: F401
        return
    except ImportError:
        pass

    search_roots = [
        Path("/opt"),
        Path("/usr/local/lib"),
    ]
    pattern = "pythia8*.so"
    candidates = []
    for root in search_roots:
        if not root.exists():
            continue
        candidates.extend(root.glob("pythia*/lib*/python*/site-packages/" + pattern))
        candidates.extend(root.glob("pythia*/lib*/site-packages/" + pattern))
        candidates.extend(root.glob("pythia*/lib/" + pattern))
        candidates.extend(root.glob("python*/dist-packages/" + pattern))
        candidates.extend(root.glob("python*/site-packages/" + pattern))

    if not candidates:
        return

    module_dir = candidates[0].parent
    module_dir_str = str(module_dir)
    if module_dir_str not in sys.path:
        sys.path.insert(0, module_dir_str)
    existing = os.environ.get("PYTHONPATH", "")
    if module_dir_str not in existing.split(os.pathsep):
        os.environ["PYTHONPATH"] = (
            module_dir_str + (os.pathsep + existing if existing else "")
        )


def hadronize_lhe(
    lhe_path: Path,
    max_events: int,
    stable_pids: Iterable[int],
    seed: int,
) -> pd.DataFrame:
    sanitized = sanitize_lhe_beams(lhe_path, lhe_path.parent)
    pythia = load_pythia()
    pythia.readString("Beams:frameType = 4")
    pythia.readString(f"Beams:LHEF = {sanitized}")
    pythia.readString("ProcessLevel:all = off")
    pythia.readString("PartonLevel:all = off")
    pythia.readString("HadronLevel:all = on")
    pythia.readString("HadronLevel:Decay = on")
    pythia.readString("ParticleDecays:limitTau0 = on")
    pythia.readString("ParticleDecays:tau0Max = 10.0")
    pythia.readString(f"Random:seed = {seed}")
    pythia.readString("Random:setSeed = on")
    pythia.init()

    rows = []
    stable_set = {int(pid) for pid in stable_pids}
    event_id = 0

    while event_id < max_events and pythia.next():
        for p in pythia.event:
            if not p.isFinal():
                continue
            pid = int(p.id())
            if pid not in stable_set:
                continue
            rows.append(
                {
                    "event": event_id,
                    "pid": pid,
                    "E": float(p.e()),
                    "px": float(p.px()),
                    "py": float(p.py()),
                    "pz": float(p.pz()),
                    "mass": float(p.m()),
                }
            )
        event_id += 1

    if event_id == 0:
        raise RuntimeError("No events processed by Pythia.")
    return pd.DataFrame(rows)


def build_stable_pid_list() -> List[int]:
    # Charged tracks + photons for DV response.
    return [
        11, -11,  # e
        13, -13,  # mu
        211, -211,  # pi+/-
        321, -321,  # K+/-
        2212, -2212,  # p
        22,  # gamma
    ]


def sanitize_lhe_beams(lhe_path: Path, work_dir: Path) -> Path:
    """
    Ensure the <init> beam IDs/energies are valid for Pythia.
    MG5 decay LHE files often set beam2 id=0, which Pythia rejects.
    """
    text = lhe_path.read_text()
    lines = text.splitlines()
    init_idx = None
    for i, line in enumerate(lines):
        if line.strip() == "<init>":
            init_idx = i
            break
    if init_idx is None or init_idx + 1 >= len(lines):
        return lhe_path

    fields = lines[init_idx + 1].split()
    if len(fields) < 10:
        return lhe_path

    id_a = fields[0]
    id_b = fields[1]
    e_a = fields[2]
    e_b = fields[3]

    allowed_beams = {"11", "-11", "13", "-13", "22", "2212", "-2212", "2112", "-2112"}
    if id_a not in allowed_beams or id_b not in allowed_beams:
        # Use "no-beams" LHEF extension for hadron-level standalone.
        fields[0] = "0"
        fields[1] = "0"
        fields[2] = "0.000000e+00"
        fields[3] = "0.000000e+00"
    elif id_b == "0":
        fields[1] = id_a
        fields[3] = e_a if e_b == "0.000000e+00" else e_b
    lines[init_idx + 1] = " ".join(fields)

    out_path = work_dir / f"{lhe_path.stem}_sanitized{lhe_path.suffix}"
    out_path.write_text("\n".join(lines) + "\n")
    return out_path


def generate_for_mass(config: RunConfig, mass_GeV: float) -> None:
    mass_label = f"{mass_GeV:.2f}".replace(".", "p")
    proc_dir = config.work_dir / f"mg5_hnl_decay_{mass_label}"
    proc_dir.mkdir(parents=True, exist_ok=True)

    proc_text = build_mg5_process_command(proc_dir)
    run_madgraph(config.mg5_path, proc_text, config.work_dir)

    # Write param_card with mixings and mass into Cards/ before launch
    template = (Path(__file__).resolve().parents[2]
                / "production" / "madgraph_production" / "cards" / "param_card_template.dat")
    param_card = proc_dir / "Cards" / "param_card.dat"
    write_param_card(template, param_card, mass_GeV, config.Ue2, config.Umu2, config.Utau2)

    launch_text = build_mg5_launch_command(proc_dir, config.nevents, config.seed)
    run_madgraph(config.mg5_path, launch_text, config.work_dir)

    lhe_path = find_lhe_file(proc_dir)
    lhe_path = ensure_lhe_uncompressed(lhe_path, config.work_dir)

    stable_pids = build_stable_pid_list()
    df = hadronize_lhe(lhe_path, config.nevents, stable_pids, config.seed)

    ctau0_m, gamma_tot = hnl_ctau_and_width(mass_GeV, config.Ue2, config.Umu2, config.Utau2)
    df["mass_GeV"] = mass_GeV
    df["ctau0_m"] = ctau0_m
    df["gamma_tot_GeV"] = gamma_tot

    out_path = config.out_dir / f"HNL_decay_rest_{mass_label}.csv"
    df.to_csv(out_path, index=False)
    print(f"Saved {len(df)} particles to {out_path}")


def parse_masses(value: str) -> List[float]:
    masses = []
    for token in value.split(","):
        token = token.strip()
        if not token:
            continue
        masses.append(float(token))
    return masses


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate HNL decay events with MG5 + Pythia8.")
    parser.add_argument("--masses", type=str, default="5.0,10.0,15.0,20.0")
    parser.add_argument("--nevents", type=int, default=20000)
    parser.add_argument("--Ue2", type=float, default=0.0)
    parser.add_argument("--Umu2", type=float, default=1e-6)
    parser.add_argument("--Utau2", type=float, default=0.0)
    parser.add_argument("--mg5-path", type=str, default="/opt/MG5_aMC_v3_6_6/bin/mg5_aMC")
    parser.add_argument("--work-dir", type=str, default="analysis_pbc/decay/work")
    parser.add_argument("--out-dir", type=str, default="analysis_pbc/decay/output")
    parser.add_argument("--seed", type=int, default=12345)
    args = parser.parse_args()

    masses = parse_masses(args.masses)
    work_dir = Path(args.work_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    config = RunConfig(
        masses=masses,
        nevents=args.nevents,
        Ue2=args.Ue2,
        Umu2=args.Umu2,
        Utau2=args.Utau2,
        mg5_path=Path(args.mg5_path),
        work_dir=work_dir,
        out_dir=out_dir,
        seed=args.seed,
    )

    if not config.mg5_path.exists():
        raise FileNotFoundError(f"MG5 executable not found: {config.mg5_path}")

    for mass in config.masses:
        print(f"=== Generating decay sample for mN={mass:.2f} GeV ===")
        generate_for_mass(config, mass)


if __name__ == "__main__":
    main()
