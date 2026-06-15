#!/usr/bin/env python3
from __future__ import annotations

import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ENV = ROOT / "env"
LINUX = ROOT / "src" / "fonll" / "Linux"
RUNDIR = ROOT / "run" / "validate_cteq66_charm_y0"
LOGDIR = ROOT / "logs"

PUBLIC_D0_PT5_Y0 = 3.5221e7


def env() -> dict[str, str]:
    out = os.environ.copy()
    out["PATH"] = f"{ENV / 'bin'}:{out.get('PATH', '')}"
    out["LHAPDF_DATA_PATH"] = str(ENV / "share" / "LHAPDF")
    return out


def run(cmd: list[str], stdin: str, log_name: str) -> None:
    LOGDIR.mkdir(parents=True, exist_ok=True)
    RUNDIR.mkdir(parents=True, exist_ok=True)
    with (LOGDIR / log_name).open("w") as log:
        proc = subprocess.run(
            cmd,
            cwd=RUNDIR,
            input=stdin,
            text=True,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env(),
            check=False,
        )
    if proc.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed with exit {proc.returncode}; see {LOGDIR / log_name}")


def grid_input() -> str:
    return "\n".join(
        [
            "val_c_",
            " 1 7000. 0 0 10550",
            " 1 7000. 0 0 10550",
            " 1.5",
            " -1.",
            " 1. 1.",
            " -0.5",
            " 0",
            " 0.5",
            " 10000",
            " 0",
            " 200 80",
            " 1",
        ]
    ) + "\n"


def frag_input() -> str:
    return "\n".join(
        [
            "0",
            "0",
            "5",
            "val_c_.out",
            "",
            "0",
            "5",
            "1",
            "-1",
            "0.1",
            "9",
            "5 0",
            "-1 0",
        ]
    ) + "\n"


def read_fragmented_point() -> float:
    path = RUNDIR / "fragmfonll.dat"
    for line in path.read_text().splitlines():
        if line.startswith("pt") or not line.strip():
            continue
        cols = line.split()
        if len(cols) >= 3:
            return float(cols[2].replace("D", "E"))
    raise ValueError(f"no data point found in {path}")


def main() -> None:
    for stale in [
        RUNDIR / "val_c_.out",
        RUNDIR / "val_c_.outlog",
        RUNDIR / "val_c_fonll.log",
        RUNDIR / "fragmfonll.dat",
        RUNDIR / "fragfonll.log",
    ]:
        if stale.exists():
            stale.unlink()

    run([str(LINUX / "fonllgridlha")], grid_input(), "validate_cteq66_charm_y0_fonllgrid.log")
    run([str(LINUX / "fragmfonll")], frag_input(), "validate_cteq66_charm_y0_fragmfonll.log")

    local = read_fragmented_point()
    rel = (local - PUBLIC_D0_PT5_Y0) / PUBLIC_D0_PT5_Y0
    print(f"local={local:.12e}")
    print(f"public_D0={PUBLIC_D0_PT5_Y0:.12e}")
    print(f"rel_diff={rel:+.6%}")


if __name__ == "__main__":
    main()
