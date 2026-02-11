#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Tuple

THIS_FILE = Path(__file__).resolve()
REPO_ROOT = THIS_FILE.parents[2]
ANALYSIS_ROOT = REPO_ROOT / "analysis_pbc"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))

from config_mass_grid import MASS_GRID
from decay.rhn_decay_library import FLAVOUR_CONFIG
from tools.decay.decay_library_io import format_mass_suffix
from tools.decay.generate_hnl_decay_events import RunConfig, generate_for_mass

ALL_FLAVOURS = ("electron", "muon", "tau")
FLAVOUR_SEED_OFFSET = {
    "electron": 1000,
    "muon": 2000,
    "tau": 3000,
}


def _parse_csv_list(raw: str) -> List[str]:
    return [x.strip() for x in raw.split(",") if x.strip()]


def _parse_mass_list(raw: str) -> List[float]:
    return [float(x) for x in _parse_csv_list(raw)]


def couplings_for_flavour(flavour: str, u2_norm: float) -> Tuple[float, float, float]:
    if flavour == "electron":
        return float(u2_norm), 0.0, 0.0
    if flavour == "muon":
        return 0.0, float(u2_norm), 0.0
    if flavour == "tau":
        return 0.0, 0.0, float(u2_norm)
    raise ValueError(f"Unknown flavour '{flavour}'")


def _overlay_out_dir(output_root: Path, flavour: str) -> Path:
    cfg = FLAVOUR_CONFIG[flavour]
    return output_root / cfg["repo"] / cfg["decay_dir"]


def _hadronized_masses(masses: List[float], flavour: str) -> List[float]:
    threshold = float(FLAVOUR_CONFIG[flavour]["low_mass_threshold"])
    return [m for m in masses if float(m) > threshold]


def _resolve_masses(args: argparse.Namespace) -> List[float]:
    if args.masses:
        masses = _parse_mass_list(args.masses)
    elif args.from_mass_grid:
        masses = list(MASS_GRID)
    else:
        raise ValueError("No masses selected. Provide --masses or keep --from-mass-grid enabled.")
    return sorted({float(m) for m in masses})


def _resolve_flavours(raw: str) -> List[str]:
    flavours = _parse_csv_list(raw)
    if not flavours:
        return list(ALL_FLAVOURS)
    unknown = [f for f in flavours if f not in ALL_FLAVOURS]
    if unknown:
        raise ValueError(f"Unknown flavour(s): {unknown}. Allowed: {ALL_FLAVOURS}")
    return flavours


def _seed_for_point(base_seed: int, flavour: str, mass_GeV: float) -> int:
    offset = FLAVOUR_SEED_OFFSET.get(flavour)
    if offset is None:
        raise ValueError(f"Unknown flavour '{flavour}' for seed derivation.")
    raw = int(base_seed) + int(round(float(mass_GeV) * 1000.0)) + int(offset)
    # Keep positive, bounded seed values that are accepted by MG5/Pythia.
    return 1 + (raw % 900_000_000)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Precompute overlay decay-library txt files for hadronized mass region."
    )
    parser.add_argument("--flavours", type=str, default="electron,muon,tau")
    parser.add_argument("--masses", type=str, default=None, help="Comma-separated mass list in GeV.")
    parser.add_argument(
        "--from-mass-grid",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use MASS_GRID when --masses is not provided (default: true).",
    )
    parser.add_argument("--nevents", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--mg5-path", type=str, default="/opt/MG5_aMC_v3_6_6/bin/mg5_aMC")
    parser.add_argument("--work-dir", type=str, default=str(REPO_ROOT / "analysis_pbc" / "decay" / "work"))
    parser.add_argument("--output-root", type=str, default=str(REPO_ROOT / "output" / "decay" / "generated"))
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--u2-norm", type=float, default=1e-6)
    args = parser.parse_args()

    flavours = _resolve_flavours(args.flavours)
    masses = _resolve_masses(args)

    mg5_path = Path(args.mg5_path).resolve()
    if not mg5_path.exists():
        raise FileNotFoundError(f"MG5 executable not found: {mg5_path}")
    if args.u2_norm <= 0.0:
        raise ValueError("--u2-norm must be positive.")
    if args.nevents <= 0:
        raise ValueError("--nevents must be positive.")

    output_root = Path(args.output_root).resolve()
    work_root = Path(args.work_dir).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    work_root.mkdir(parents=True, exist_ok=True)

    total_jobs = 0
    completed = 0
    skipped = 0

    for flavour in flavours:
        mass_list = _hadronized_masses(masses, flavour)
        if not mass_list:
            print(f"[{flavour}] No masses above hadronized threshold; skipping.")
            continue

        Ue2, Umu2, Utau2 = couplings_for_flavour(flavour, args.u2_norm)
        out_dir = _overlay_out_dir(output_root, flavour)
        work_dir = work_root / flavour
        out_dir.mkdir(parents=True, exist_ok=True)
        work_dir.mkdir(parents=True, exist_ok=True)

        print(
            f"[{flavour}] masses={len(mass_list)} u2_norm={args.u2_norm:g} "
            f"(Ue2,Umu2,Utau2)=({Ue2:g},{Umu2:g},{Utau2:g})"
        )
        for mass in mass_list:
            total_jobs += 1
            txt_prefix = "vN_Ntoall_generated"
            out_txt = out_dir / f"{txt_prefix}_{format_mass_suffix(mass)}.txt"
            if out_txt.exists() and not args.overwrite:
                print(f"  [skip] exists: {out_txt.name}")
                skipped += 1
                continue
            point_seed = _seed_for_point(args.seed, flavour, mass)
            config = RunConfig(
                masses=[mass],
                nevents=args.nevents,
                Ue2=Ue2,
                Umu2=Umu2,
                Utau2=Utau2,
                mg5_path=mg5_path,
                work_dir=work_dir,
                out_dir=out_dir,
                seed=point_seed,
                output_format="txt",
                final_state_mode="all",
                txt_prefix=txt_prefix,
            )
            print(f"  [run ] mN={mass:.6g} GeV seed={point_seed} -> {out_txt.name}")
            generate_for_mass(config, mass)
            completed += 1

    print(
        f"Done. total={total_jobs} completed={completed} skipped={skipped} "
        f"output_root={output_root}"
    )


if __name__ == "__main__":
    main()
