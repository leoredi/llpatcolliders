#!/usr/bin/env python3
"""
production/combine_channels.py

Concatenate all available channel CSVs for each (flavor, mass) into combined files.

Input directories:
  output/llp_4vectors/{flavor}/{Bmeson,Dmeson,Bc,Kmeson,
                                tau,induced_tau,WZ}/mN_{mass}.csv

Output:
  output/llp_4vectors/{flavor}/combined/mN_{mass}.csv

Format: headerless, 5 columns: weight,E,px,py,pz
"""

import sys
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from config_mass_grid import MASS_GRID, format_mass_for_filename

OUTPUT_BASE = PROJECT_ROOT / "output" / "llp_4vectors"
# Channel inventory semantics:
#   - existing file, zero bytes  -> channel closed at this m_N (kinematics).
#   - existing file, >0 bytes    -> real contribution, read and concatenate.
#   - file missing               -> channel never ran for this point.
# combine_for_point distinguishes the latter two when strict=True (used by
# run_all so a silent driver failure cannot pass off an incomplete combined
# file as success). Tau is split into two physical sources written to
# different folders by their respective drivers:
#   tau/         — prompt: pp -> W -> tau nu  /  pp -> Z -> tau tau (MG5)
#   induced_tau/ — cascade: Ds, B+ -> tau nu -> tau -> N (FONLL meson source)
CHANNELS = ["Bmeson", "Dmeson", "Bc", "tau", "induced_tau", "Kmeson", "WZ"]
FLAVORS = ["Ue", "Umu", "Utau"]


def combine_for_point(flavor, mass, channels=None, strict=False):
    """
    Combine all channel CSVs for a single (flavor, mass) point.

    Parameters
    ----------
    flavor, mass
        Flavor label and HNL mass.
    channels : list of str, optional
        Per-call channel inventory (default: the module-level CHANNELS).
        Callers that intentionally skip a channel (``--no-wz``,
        ``--no-prompt-tau``) should pass the reduced list so a legitimately
        absent channel is not flagged as missing in strict mode.
    strict : bool
        If True, raise RuntimeError when any requested channel CSV is missing
        or unreadable. ``empty`` (zero-byte) files are always treated as
        "channel closed" and never trigger strict failure.

    Returns
    -------
    int
        Total number of events in combined file.
    """
    mass_label = format_mass_for_filename(mass)
    combined_dir = OUTPUT_BASE / flavor / "combined"
    combined_dir.mkdir(parents=True, exist_ok=True)
    combined_path = combined_dir / f"mN_{mass_label}.csv"

    channels = CHANNELS if channels is None else channels
    all_data, missing, bad = [], [], []
    for ch in channels:
        csv_path = OUTPUT_BASE / flavor / ch / f"mN_{mass_label}.csv"
        if not csv_path.exists():
            missing.append(ch)
            continue
        if csv_path.stat().st_size == 0:
            # Empty sentinel: channel closed at this m_N. Legitimate zero.
            continue
        try:
            data = np.loadtxt(csv_path, delimiter=",")
            if data.ndim == 1:
                data = data.reshape(1, -1)
            if data.shape[1] >= 5:
                all_data.append(data[:, :5])
            else:
                bad.append(f"{ch}: only {data.shape[1]} columns")
        except Exception as e:
            bad.append(f"{ch}: {e}")

    if strict and (missing or bad):
        raise RuntimeError(
            f"incomplete combine for {flavor} mN={mass_label}: "
            f"missing={missing}, bad={bad}"
        )

    if all_data:
        combined = np.vstack(all_data)
        np.savetxt(combined_path, combined, delimiter=",", fmt="%.8e")
        return len(combined)
    else:
        combined_path.write_text("")
        return 0


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Combine HNL production channels")
    parser.add_argument("--flavor", choices=FLAVORS, nargs="+", default=FLAVORS)
    parser.add_argument("--masses", type=float, nargs="+", default=None)
    args = parser.parse_args()

    masses = args.masses if args.masses else MASS_GRID

    total_files = 0
    total_events = 0
    n_empty = 0

    for flavor in args.flavor:
        print(f"\n{flavor}:")
        for mass in masses:
            n_ev = combine_for_point(flavor, mass)
            total_files += 1
            total_events += n_ev
            if n_ev == 0:
                n_empty += 1

        # Summary for this flavor
        combined_dir = OUTPUT_BASE / flavor / "combined"
        n_files = len(list(combined_dir.glob("mN_*.csv")))
        print(f"  {n_files} files in {combined_dir}")

    print(f"\nTotal: {total_files} combined files, {total_events} events, {n_empty} empty")


if __name__ == "__main__":
    main()
