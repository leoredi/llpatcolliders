#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List

import numpy as np
import pandas as pd

import sys

THIS_FILE = Path(__file__).resolve()
REPO_ROOT = THIS_FILE.parents[2]
ANALYSIS_ROOT = REPO_ROOT / "analysis_pbc"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from hnl_models.hnl_model_hnlcalc import HNLModel
from tools.decay.decay_library_io import (
    dataframe_to_event_daughters,
    format_mass_suffix,
    write_decay_library_txt,
)

HBARC_GEV_M = 1.973269804e-16


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
    output_format: str
    final_state_mode: str
    txt_prefix: str


def hnl_ctau_and_width(mass_GeV: float, Ue2: float, Umu2: float, Utau2: float) -> tuple[float, float]:
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
add process n1 > vl vl vl~
add process n1 > vl~ vl~ vl

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
        import pythia8
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


def _parse_lhe_events(lhe_path: Path) -> List[List[dict]]:
    """Parse LHE file into a list of events, each a list of particle dicts."""
    events: List[List[dict]] = []
    in_event = False
    current: List[dict] = []

    with open(lhe_path) as f:
        for line in f:
            stripped = line.strip()
            if stripped == "<event>":
                in_event = True
                current = []
                continue
            if stripped == "</event>":
                in_event = False
                if current:
                    events.append(current)
                continue
            if not in_event or stripped.startswith("#") or stripped.startswith("<"):
                continue
            parts = stripped.split()
            if len(parts) < 13:
                continue
            current.append({
                "pid": int(parts[0]),
                "lhe_status": int(parts[1]),
                "col": int(parts[4]),
                "acol": int(parts[5]),
                "px": float(parts[6]),
                "py": float(parts[7]),
                "pz": float(parts[8]),
                "e": float(parts[9]),
                "m": float(parts[10]),
            })
    return events


def _needs_hadronization(pid: int) -> bool:
    absid = abs(pid)
    return (1 <= absid <= 6) or absid == 21


_UNSTABLE_LEPTONS = {15}  # tau


def _needs_pythia(pid: int) -> bool:
    """True if this particle requires Pythia processing (hadronization or decay)."""
    return _needs_hadronization(pid) or abs(pid) in _UNSTABLE_LEPTONS


def hadronize_lhe(
    lhe_path: Path,
    max_events: int,
    seed: int,
    selected_pids: Iterable[int] | None = None,
) -> pd.DataFrame:
    # Parse LHE events directly instead of using Pythia's LHEF reader,
    # which fails on decay-at-rest events with BSM beam IDs.
    lhe_events = _parse_lhe_events(lhe_path)
    if max_events > 0:
        lhe_events = lhe_events[:max_events]

    pythia = load_pythia()
    # Skip process generation; we fill the event record manually.
    pythia.readString("ProcessLevel:all = off")
    pythia.readString("HadronLevel:all = on")
    pythia.readString("HadronLevel:Decay = on")
    pythia.readString("ParticleDecays:limitTau0 = on")
    pythia.readString("ParticleDecays:tau0Max = 10.0")
    pythia.readString("Random:setSeed = on")
    pythia.readString(f"Random:seed = {seed}")
    pythia.readString("Next:numberShowInfo = 0")
    pythia.readString("Next:numberShowProcess = 0")
    pythia.readString("Next:numberShowEvent = 0")
    pythia.init()

    rows: List[dict] = []
    selected_set = {int(pid) for pid in selected_pids} if selected_pids is not None else None
    n_failed = 0

    for event_id, particles in enumerate(lhe_events):
        final = [p for p in particles if p["lhe_status"] == 1]
        needs_pythia = any(_needs_pythia(p["pid"]) for p in final)

        if not needs_pythia:
            # Truly stable leptonic event (e/mu + neutrinos only) — safe to bypass.
            for p in final:
                pid = p["pid"]
                if selected_set is not None and pid not in selected_set:
                    continue
                rows.append({
                    "event": event_id, "pid": pid,
                    "E": p["e"], "px": p["px"], "py": p["py"],
                    "pz": p["pz"], "mass": p["m"],
                })
            continue

        # Fill Pythia event record — hadronize partons and/or decay unstable leptons.
        pythia.event.reset()
        for p in final:
            # Status 23 = outgoing from hardest subprocess.  Pythia
            # hadronizes colored particles (col/acol != 0) and leaves
            # color-singlet leptons untouched.
            pythia.event.append(
                p["pid"], 23, p["col"], p["acol"],
                p["px"], p["py"], p["pz"], p["e"], p["m"],
            )

        if not pythia.next():
            n_failed += 1
            continue

        for i in range(pythia.event.size()):
            if not pythia.event[i].isFinal():
                continue
            pid = int(pythia.event[i].id())
            if selected_set is not None and pid not in selected_set:
                continue
            rows.append({
                "event": event_id, "pid": pid,
                "E": float(pythia.event[i].e()),
                "px": float(pythia.event[i].px()),
                "py": float(pythia.event[i].py()),
                "pz": float(pythia.event[i].pz()),
                "mass": float(pythia.event[i].m()),
            })

    n_total = len(lhe_events)
    if n_total == 0:
        raise RuntimeError("No LHE events found in file.")

    if n_failed > 0:
        warnings.warn(
            f"hadronize_lhe: Pythia hadronization failed for "
            f"{n_failed}/{n_total} events.",
            UserWarning,
        )

    n_with_daughters = len({r["event"] for r in rows}) if rows else 0
    if max_events > 0 and n_with_daughters < 0.9 * max_events:
        warnings.warn(
            f"hadronize_lhe: only {n_with_daughters}/{max_events} events "
            f"produced usable daughters ({n_with_daughters / max_events:.1%}).",
            UserWarning,
        )

    return pd.DataFrame(rows)


def build_stable_pid_list() -> List[int]:
    # Legacy subset used by the old CSV workflow.
    return [
        11, -11,
        13, -13,
        211, -211,
        321, -321,
        2212, -2212,
        22,
    ]


def _selected_pids_for_mode(mode: str) -> List[int] | None:
    if mode == "all":
        return None
    if mode == "stable_subset":
        return build_stable_pid_list()
    raise ValueError(f"Unknown final-state mode '{mode}'. Use 'all' or 'stable_subset'.")


def generate_for_mass(config: RunConfig, mass_GeV: float) -> None:
    mass_label = f"{mass_GeV:.2f}".replace(".", "p")
    proc_dir = config.work_dir / f"mg5_hnl_decay_{mass_label}"
    proc_dir.mkdir(parents=True, exist_ok=True)

    proc_text = build_mg5_process_command(proc_dir)
    run_madgraph(config.mg5_path, proc_text, config.work_dir)

    template = REPO_ROOT / "production" / "madgraph_production" / "cards" / "param_card_template.dat"
    param_card = proc_dir / "Cards" / "param_card.dat"
    write_param_card(template, param_card, mass_GeV, config.Ue2, config.Umu2, config.Utau2)

    launch_text = build_mg5_launch_command(proc_dir, config.nevents, config.seed)
    run_madgraph(config.mg5_path, launch_text, config.work_dir)

    lhe_path = find_lhe_file(proc_dir)
    lhe_path = ensure_lhe_uncompressed(lhe_path, config.work_dir)

    selected_pids = _selected_pids_for_mode(config.final_state_mode)
    df = hadronize_lhe(lhe_path, config.nevents, config.seed, selected_pids=selected_pids)

    ctau0_m, gamma_tot = hnl_ctau_and_width(mass_GeV, config.Ue2, config.Umu2, config.Utau2)

    if config.output_format in {"csv", "both"}:
        csv_df = df.copy()
        csv_df["mass_GeV"] = mass_GeV
        csv_df["ctau0_m"] = ctau0_m
        csv_df["gamma_tot_GeV"] = gamma_tot
        out_csv = config.out_dir / f"HNL_decay_rest_{mass_label}.csv"
        csv_df.to_csv(out_csv, index=False)
        print(f"Saved {len(csv_df)} particles to {out_csv}")

    if config.output_format in {"txt", "both"}:
        event_daughters = dataframe_to_event_daughters(df)
        mass_suffix = format_mass_suffix(mass_GeV)
        out_txt = config.out_dir / f"{config.txt_prefix}_{mass_suffix}.txt"
        n_events = write_decay_library_txt(out_txt, mass_GeV=mass_GeV, event_daughters=event_daughters)
        print(f"Saved {n_events} decay events to {out_txt}")


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
    parser.add_argument("--work-dir", type=str, default=str(REPO_ROOT / "analysis_pbc" / "decay" / "work"))
    parser.add_argument("--out-dir", type=str, default=str(REPO_ROOT / "analysis_pbc" / "decay" / "output"))
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument(
        "--output-format",
        choices=["csv", "txt", "both"],
        default="csv",
        help="Output format: csv (legacy), txt (decay library), or both.",
    )
    parser.add_argument(
        "--final-state-mode",
        choices=["stable_subset", "all"],
        default="stable_subset",
        help="Particle selection for output rows; 'all' is recommended for decay-library txt output.",
    )
    parser.add_argument(
        "--txt-prefix",
        type=str,
        default="vN_Ntoall_generated",
        help="Filename prefix for txt mode (mass suffix is appended automatically).",
    )
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
        output_format=args.output_format,
        final_state_mode=args.final_state_mode,
        txt_prefix=args.txt_prefix,
    )

    if not config.mg5_path.exists():
        raise FileNotFoundError(f"MG5 executable not found: {config.mg5_path}")

    for mass in config.masses:
        print(f"=== Generating decay sample for mN={mass:.2f} GeV ===")
        generate_for_mass(config, mass)


if __name__ == "__main__":
    main()
