from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, List, Tuple

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
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
        "lightfstau",
        "lightfstauK",
    ),
}

DECAY_CATEGORY_ORDER = (
    "lightfstauK",
    "lightfstau",
    "lightfonly",
    "inclDs",
    "inclDD",
    "inclD",
    "nocharmnoss",
    "nocharm",
    "analytical2and3bodydecays",
)

MAX_DECAY_FILE_DELTA_GEV = 0.5


@dataclass(frozen=True)
class DecayFileEntry:
    path: Path
    mass_GeV: float
    category: str


def _decay_dir(flavour: str) -> Path:
    config = FLAVOUR_CONFIG.get(flavour)
    if not config:
        raise ValueError(f"Unknown flavour '{flavour}' for decay library.")
    return EXTERNAL_ROOT / config["repo"] / config["decay_dir"]


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


@lru_cache(maxsize=6)
def list_decay_files(flavour: str) -> List[DecayFileEntry]:
    decay_dir = _decay_dir(flavour)
    if not decay_dir.exists():
        raise FileNotFoundError(
            f"Decay directory not found: {decay_dir}. "
            "Make sure the MATHUSLA_LLPfiles RHN repositories are present."
        )
    entries: List[DecayFileEntry] = []
    for path in decay_dir.glob("*.txt"):
        mass = _parse_mass_from_name(path)
        if mass is None:
            continue
        entries.append(DecayFileEntry(path=path, mass_GeV=mass, category=_categorize_name(path)))
    if not entries:
        raise FileNotFoundError(f"No decay files found in {decay_dir}.")
    return entries


def _nearest_entry(entries: Iterable[DecayFileEntry], mass_GeV: float) -> DecayFileEntry | None:
    entries = list(entries)
    if not entries:
        return None
    return min(entries, key=lambda e: abs(e.mass_GeV - mass_GeV))


def _warn_if_large_mismatch(chosen: DecayFileEntry, mass_GeV: float, max_delta: float) -> None:
    delta = abs(chosen.mass_GeV - mass_GeV)
    if delta > max_delta:
        warnings.warn(
            f"Large decay-file mass mismatch for {mass_GeV:.3f} GeV: "
            f"selected {chosen.mass_GeV:.3f} GeV (Δ={delta:.3f} GeV) "
            f"category={chosen.category} file={chosen.path.name}",
            UserWarning,
        )


def select_decay_file(flavour: str, mass_GeV: float) -> DecayFileEntry:
    entries = list_decay_files(flavour)
    low_mass_threshold = FLAVOUR_CONFIG[flavour]["low_mass_threshold"]
    priorities = DECAY_PRIORITIES.get(flavour)
    if not priorities:
        chosen = _nearest_entry(entries, mass_GeV)
        if chosen is None:
            raise FileNotFoundError(f"No decay files available for flavour '{flavour}'.")
        _warn_if_large_mismatch(chosen, mass_GeV, MAX_DECAY_FILE_DELTA_GEV)
        return chosen

    if mass_GeV <= low_mass_threshold:
        analytical = [e for e in entries if e.category == "analytical2and3bodydecays"]
        chosen = _nearest_entry(analytical, mass_GeV)
        if chosen is None:
            chosen = _nearest_entry(entries, mass_GeV)
        if chosen is None:
            raise FileNotFoundError(f"No decay files available for flavour '{flavour}'.")
        _warn_if_large_mismatch(chosen, mass_GeV, MAX_DECAY_FILE_DELTA_GEV)
        return chosen

    allowed = [e for e in entries if e.category in priorities and e.category != "analytical2and3bodydecays"]
    chosen = _nearest_entry(allowed, mass_GeV)
    if chosen is None:
        chosen = _nearest_entry(entries, mass_GeV)
    if chosen is None:
        raise FileNotFoundError(f"No decay files available for flavour '{flavour}'.")
    _warn_if_large_mismatch(chosen, mass_GeV, MAX_DECAY_FILE_DELTA_GEV)
    return chosen


def _parse_decay_event_block(lines: List[str]) -> List[Tuple[float, float, float, float, float, int]]:
    daughters: List[Tuple[float, float, float, float, float, int]] = []
    if not lines:
        return daughters
    for line in lines[1:]:
        parts = [p.strip() for p in line.split(",") if p.strip()]
        if len(parts) < 6:
            continue
        try:
            E, px, py, pz, mass, pid = parts[:6]
            daughters.append((float(E), float(px), float(py), float(pz), float(mass), int(pid)))
        except ValueError:
            continue
    return daughters


@lru_cache(maxsize=32)
def load_decay_events(path: Path) -> List[List[Tuple[float, float, float, float, float, int]]]:
    if not path.exists():
        raise FileNotFoundError(f"Decay file not found: {path}")
    text = path.read_text()
    lines = [line.strip() for line in text.splitlines()]

    events: List[List[Tuple[float, float, float, float, float, int]]] = []
    current: List[str] = []
    for line in lines:
        if not line:
            if current:
                event = _parse_decay_event_block(current)
                if event:
                    events.append(event)
                current = []
            continue
        if line.lower().startswith("format is"):
            continue
        current.append(line)
    if current:
        event = _parse_decay_event_block(current)
        if event:
            events.append(event)
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
