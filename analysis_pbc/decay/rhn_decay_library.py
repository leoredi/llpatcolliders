from __future__ import annotations

import os
import re
import warnings
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, List, Tuple

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
GENERATED_ROOT = REPO_ROOT / "output" / "decay" / "generated"
EXTERNAL_ROOT = REPO_ROOT / "analysis_pbc" / "decay" / "external"

FLAVOUR_CONFIG = {
    "electron": {
        "repo": "MATHUSLA_LLPfiles_RHN_Ue",
        "decay_dir": "RHN_Ue_hadronic_decays_geant",
        "low_mass_threshold": 0.42,
    },
    "muon": {
        "repo": "MATHUSLA_LLPfiles_RHN_Umu",
        "decay_dir": "RHN_Umu_hadronic_decays_geant",
        "low_mass_threshold": 0.53,
    },
    "tau": {
        "repo": "MATHUSLA_LLPfiles_RHN_Utau",
        "decay_dir": "RHN_Utau_hadronic_decays_geant",
        "low_mass_threshold": 0.42,
    },
}

DECAY_PRIORITIES = {
    "electron": (
        "inclDs",
        "inclDD",
        "inclD",
        "nocharm",
        "nocharmnoss",
        "lightfonly",
        "analytical2and3bodydecays",
    ),
    "muon": (
        "inclDs",
        "inclDD",
        "inclD",
        "nocharm",
        "nocharmnoss",
        "lightfonly",
        "analytical2and3bodydecays",
    ),
    "tau": (
        "lightfonly",
        "lightfsonly",
        "lightfstau",
        "lightfstauK",
    ),
}

DECAY_CATEGORY_ORDER = (
    "lightfstauK",
    "lightfstau",
    "lightfsonly",
    "lightfonly",
    "inclDs",
    "inclDD",
    "inclD",
    "nocharmnoss",
    "nocharm",
    "analytical2and3bodydecays",
)

MAX_DECAY_FILE_DELTA_GEV = 0.5
OVERLAY_SWITCH_MASS_GEV = 5.0


@dataclass(frozen=True)
class DecayFileEntry:
    path: Path
    mass_GeV: float
    category: str
    source: str
    source_priority: int


def _decay_relpath(flavour: str) -> Path:
    config = FLAVOUR_CONFIG.get(flavour)
    if not config:
        raise ValueError(f"Unknown flavour '{flavour}' for decay library.")
    return Path(config["repo"]) / config["decay_dir"]


def _decay_dirs_with_priority(flavour: str) -> List[Tuple[Path, str, int]]:
    rel = _decay_relpath(flavour)
    return [
        (GENERATED_ROOT / rel, "generated", 0),
        (EXTERNAL_ROOT / rel, "external", 1),
    ]


def _parse_mass_from_name(path: Path) -> float | None:
    match = re.search(r"_([0-9]+(?:\.[0-9]+)?)\.txt$", path.name)
    if not match:
        return None
    return float(match.group(1))


def _categorize_name(path: Path) -> str:
    name = path.name
    for category in DECAY_CATEGORY_ORDER:
        if category in name:
            return category
    return "unknown"


def _category_rank(category: str) -> int:
    try:
        return DECAY_CATEGORY_ORDER.index(category)
    except ValueError:
        return len(DECAY_CATEGORY_ORDER)


def _entry_order_key(entry: DecayFileEntry, mass_GeV: float) -> Tuple[float, int, int, str]:
    return (
        abs(entry.mass_GeV - mass_GeV),
        entry.source_priority,
        _category_rank(entry.category),
        entry.path.name,
    )


@lru_cache(maxsize=6)
def list_decay_files(flavour: str) -> List[DecayFileEntry]:
    seen_names: set[str] = set()
    entries: List[DecayFileEntry] = []
    searched_dirs: List[Path] = []
    for decay_dir, source, source_priority in _decay_dirs_with_priority(flavour):
        searched_dirs.append(decay_dir)
        if not decay_dir.exists():
            continue
        for path in sorted(decay_dir.glob("*.txt")):
            # Overlay entries shadow external entries with identical filenames.
            if path.name in seen_names:
                continue
            mass = _parse_mass_from_name(path)
            if mass is None:
                continue
            seen_names.add(path.name)
            entries.append(
                DecayFileEntry(
                    path=path,
                    mass_GeV=mass,
                    category=_categorize_name(path),
                    source=source,
                    source_priority=source_priority,
                )
            )
    if not entries:
        searched_str = ", ".join(str(d) for d in searched_dirs)
        raise FileNotFoundError(f"No decay files found for flavour '{flavour}'. Searched: {searched_str}")
    return entries


def _nearest_entry(entries: Iterable[DecayFileEntry], mass_GeV: float) -> DecayFileEntry | None:
    entries = list(entries)
    if not entries:
        return None
    return min(entries, key=lambda e: _entry_order_key(e, mass_GeV))


def _allow_large_mismatch_from_env() -> bool:
    raw = os.environ.get("HNL_ALLOW_DECAY_MASS_MISMATCH", "").strip().lower()
    return raw in {"1", "true", "yes", "on"}


def _enforce_mass_mismatch_policy(chosen: DecayFileEntry, mass_GeV: float, max_delta: float) -> None:
    delta = abs(chosen.mass_GeV - mass_GeV)
    if delta <= max_delta:
        return

    message = (
        f"Large decay-file mass mismatch for {mass_GeV:.3f} GeV: "
        f"selected {chosen.mass_GeV:.3f} GeV (Δ={delta:.3f} GeV) "
        f"category={chosen.category} source={chosen.source} file={chosen.path}"
    )
    if _allow_large_mismatch_from_env():
        warnings.warn(f"{message} (allowed by HNL_ALLOW_DECAY_MASS_MISMATCH)", UserWarning)
        return
    raise ValueError(
        f"{message}. Refusing to extrapolate beyond {max_delta:.3f} GeV. "
        "Set HNL_ALLOW_DECAY_MASS_MISMATCH=1 to override for diagnostics."
    )


def _parse_pid_token(token: str, tol: float = 1e-6) -> int:
    """
    Parse PID tokens encoded as either integer text ("16") or integral float
    text ("16.0"). Reject non-integral values (e.g. "16.5").
    """
    pid_float = float(token)
    pid_int = int(round(pid_float))
    if abs(pid_float - pid_int) > tol:
        raise ValueError(f"Non-integral PID token: {token}")
    return pid_int


def select_decay_file(flavour: str, mass_GeV: float) -> DecayFileEntry:
    entries = list_decay_files(flavour)
    low_mass_threshold = FLAVOUR_CONFIG[flavour]["low_mass_threshold"]
    priorities = DECAY_PRIORITIES.get(flavour)
    if not priorities:
        chosen = _nearest_entry(entries, mass_GeV)
        if chosen is None:
            raise FileNotFoundError(f"No decay files available for flavour '{flavour}'.")
        _enforce_mass_mismatch_policy(chosen, mass_GeV, MAX_DECAY_FILE_DELTA_GEV)
        return chosen

    if mass_GeV <= low_mass_threshold:
        analytical = [e for e in entries if e.category == "analytical2and3bodydecays"]
        chosen = _nearest_entry(analytical, mass_GeV)
        if chosen is None:
            chosen = _nearest_entry(entries, mass_GeV)
        if chosen is None:
            raise FileNotFoundError(f"No decay files available for flavour '{flavour}'.")
        _enforce_mass_mismatch_policy(chosen, mass_GeV, MAX_DECAY_FILE_DELTA_GEV)
        return chosen

    if mass_GeV < OVERLAY_SWITCH_MASS_GEV:
        external = [e for e in entries if e.source == "external"]
        allowed = [
            e
            for e in external
            if e.category in priorities and e.category != "analytical2and3bodydecays"
        ]
        chosen = _nearest_entry(allowed, mass_GeV)
        if chosen is None:
            chosen = _nearest_entry(external, mass_GeV)
        if chosen is None:
            overlay = [e for e in entries if e.source == "generated"]
            chosen = _nearest_entry(overlay, mass_GeV)
            if chosen is not None:
                warnings.warn(
                    "No external decay files available below overlay switch mass; "
                    f"falling back to generated source for flavour={flavour} mass={mass_GeV:.3f} GeV.",
                    UserWarning,
                )
        if chosen is None:
            raise FileNotFoundError(f"No decay files available for flavour '{flavour}'.")
        _enforce_mass_mismatch_policy(chosen, mass_GeV, MAX_DECAY_FILE_DELTA_GEV)
        return chosen

    # For masses at/above the switch, require generated overlay files and
    # treat them as all-inclusive decay samples.
    overlay = [e for e in entries if e.source == "generated"]
    if not overlay:
        raise FileNotFoundError(
            "No generated decay overlay files found for "
            f"flavour='{flavour}' at mass {mass_GeV:.3f} GeV (required for m >= {OVERLAY_SWITCH_MASS_GEV:.1f} GeV). "
            "Run tools/decay/precompute_decay_library_overlay.py to populate output/decay/generated/."
        )
    chosen = _nearest_entry(overlay, mass_GeV)
    if chosen is None:
        raise FileNotFoundError(f"No decay files available for flavour '{flavour}'.")
    _enforce_mass_mismatch_policy(chosen, mass_GeV, MAX_DECAY_FILE_DELTA_GEV)
    return chosen


def _parse_decay_event_block(lines: List[str]) -> Tuple[List[Tuple[float, float, float, float, float, int]], int]:
    daughters: List[Tuple[float, float, float, float, float, int]] = []
    malformed_rows = 0
    if not lines:
        return daughters, malformed_rows
    for line in lines[1:]:
        parts = [p.strip() for p in line.split(",") if p.strip()]
        if len(parts) < 6:
            continue
        try:
            E, px, py, pz, mass, pid = parts[:6]
            daughters.append(
                (
                    float(E),
                    float(px),
                    float(py),
                    float(pz),
                    float(mass),
                    _parse_pid_token(pid),
                )
            )
        except ValueError:
            malformed_rows += 1
            continue
    return daughters, malformed_rows


@lru_cache(maxsize=32)
def load_decay_events(path: Path) -> List[List[Tuple[float, float, float, float, float, int]]]:
    if not path.exists():
        raise FileNotFoundError(f"Decay file not found: {path}")
    text = path.read_text()
    lines = [line.strip() for line in text.splitlines()]

    events: List[List[Tuple[float, float, float, float, float, int]]] = []
    current: List[str] = []
    malformed_rows = 0
    for line in lines:
        if not line:
            if current:
                event, malformed = _parse_decay_event_block(current)
                malformed_rows += malformed
                if event:
                    events.append(event)
                current = []
            continue
        if line.lower().startswith("format is"):
            continue
        current.append(line)
    if current:
        event, malformed = _parse_decay_event_block(current)
        malformed_rows += malformed
        if event:
            events.append(event)
    if malformed_rows > 0:
        warnings.warn(
            f"Skipped {malformed_rows} malformed daughter row(s) while parsing decay file: {path}",
            UserWarning,
        )
    if not events:
        raise ValueError(f"No decay events parsed from {path}")
    return events


def pick_decay_events(
    rng: np.random.Generator,
    events: List[List[Tuple[float, float, float, float, float, int]]],
    n: int,
) -> List[List[Tuple[float, float, float, float, float, int]]]:
    indices = rng.integers(0, len(events), size=n)
    return [events[i] for i in indices]
