#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import pandas as pd

LIBRARY_HEADER = (
    "Format is groups of {llp fv (in pythia frame), mass, pid}, "
    "{daughter1 fv (in llp rest frame), pid, name}, "
    "{daughter2 fv (in llp rest frame), pid, name}..."
)

# Each daughter row is:
#   E, px, py, pz, mass, pid[, name]
DecayDaughter = Tuple[float, float, float, float, float, int]


def format_mass_suffix(mass_GeV: float, max_decimals: int = 6) -> str:
    """Return mass suffix compatible with existing decay-library filenames."""
    text = f"{float(mass_GeV):.{max_decimals}f}"
    text = text.rstrip("0").rstrip(".")
    return text if text else "0"


def dataframe_to_event_daughters(df: pd.DataFrame) -> List[List[DecayDaughter]]:
    """Convert event-level dataframe into daughters grouped by event id."""
    required = {"event", "E", "px", "py", "pz", "mass", "pid"}
    missing = required - set(df.columns)
    if missing:
        missing_str = ", ".join(sorted(missing))
        raise ValueError(f"Missing required columns for decay library conversion: {missing_str}")

    out: List[List[DecayDaughter]] = []
    if len(df) == 0:
        return out

    for event_id in sorted(df["event"].unique()):
        part = df[df["event"] == event_id]
        daughters: List[DecayDaughter] = []
        for row in part.itertuples(index=False):
            daughters.append(
                (
                    float(getattr(row, "E")),
                    float(getattr(row, "px")),
                    float(getattr(row, "py")),
                    float(getattr(row, "pz")),
                    float(getattr(row, "mass")),
                    int(getattr(row, "pid")),
                )
            )
        if daughters:
            out.append(daughters)
    return out


def _format_number(value: float) -> str:
    return f"{float(value):.16g}"


def _format_parent_line(mass_GeV: float, parent_pid: int) -> str:
    # Placeholder LLP 4-vector in "pythia frame". The loader ignores this row.
    return ",".join(
        [
            _format_number(mass_GeV),
            _format_number(0.0),
            _format_number(0.0),
            _format_number(0.0),
            _format_number(mass_GeV),
            str(int(parent_pid)),
        ]
    )


def _format_daughter_line(d: Sequence[object]) -> str:
    if len(d) < 6:
        raise ValueError(f"Expected at least 6 daughter fields, got {len(d)}")
    pieces = [
        _format_number(float(d[0])),
        _format_number(float(d[1])),
        _format_number(float(d[2])),
        _format_number(float(d[3])),
        _format_number(float(d[4])),
        str(int(d[5])),
    ]
    if len(d) >= 7 and d[6] not in (None, ""):
        pieces.append(str(d[6]))
    return ",".join(pieces)


def write_decay_library_txt(
    output_path: Path,
    mass_GeV: float,
    event_daughters: Iterable[Sequence[Sequence[object]]],
    parent_pid: int = 9900012,
) -> int:
    """
    Write a decay library .txt in the format consumed by rhn_decay_library.py.

    Returns number of event blocks written.
    """
    blocks: List[str] = []
    for daughters in event_daughters:
        if not daughters:
            continue
        lines = [_format_parent_line(mass_GeV, parent_pid)]
        lines.extend(_format_daughter_line(d) for d in daughters)
        blocks.append("\n".join(lines))

    if not blocks:
        raise ValueError("No non-empty decay events to write.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = LIBRARY_HEADER + "\n" + ("\n\n\n".join(blocks)) + "\n"
    output_path.write_text(payload)
    return len(blocks)

