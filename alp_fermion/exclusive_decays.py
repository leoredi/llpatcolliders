"""Machine-readable exclusive arXiv:2501.04525 ALP decay definitions."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

DATA_DIR = Path(__file__).resolve().parent / "data" / "senscalc_2501"
HBAR_C_GEV_M = 1.973269804e-16
INV_F_REF = 1.0e-3

# Stable or Pythia-known primary PDG products corresponding positionally to
# decay_channels.json. Pythia subsequently decays unstable mesons and taus and
# showers/hadronizes quark or gluon pairs.
DECAY_PRODUCTS_PDG = {
    "channel_001": (-11, 11),
    "channel_002": (-13, 13),
    "channel_003": (-15, 15),
    "channel_004": (22, 22),
    "channel_005": (211, -211, 111),
    "channel_006": (211, -211, 22),
    "channel_007": (211, 111, -211, 111),
    "channel_008": (211, -211, 211, -211),
    "channel_009": (221, 111, 111),
    "channel_010": (221, 211, -211),
    "channel_011": (323, -323),
    "channel_012": (313, -313),
    "channel_013": (223, 211, -211),
    "channel_014": (111, 111, 111),
    "channel_015": (331, 111, 111),
    "channel_016": (331, 211, -211),
    "channel_017": (223, 223),
    "channel_018": (21, 21),
    "channel_019": (4, -4),
    "channel_020": (3, -3),
    "channel_021": (130, 130, 111),
    "channel_022": (310, 310, 111),
    "channel_023": (130, 310, 111),
    "channel_024": (313, -313),
    "channel_025": (-321, 130, 211),
    "channel_026": (321, 130, -211),
    "channel_027": (-321, 310, 211),
    "channel_028": (321, 310, -211),
    "channel_029": (323, -323),
    "channel_030": (321, -321, 111),
    "channel_031": (113, 113),
    "channel_032": (213, -213),
}

DUPLICATE_CHANNEL_IDS = {"channel_024", "channel_029"}
UNCLASSIFIED_CHANNEL_ID = "unclassified_neutral_remainder"

_BRANCHING_TABLE = None
_TOTAL_WIDTH_TABLE = None


def _branching_table():
    global _BRANCHING_TABLE
    if _BRANCHING_TABLE is None:
        _BRANCHING_TABLE = np.genfromtxt(
            DATA_DIR / "branching_ratios.csv", delimiter=",", names=True
        )
    return _BRANCHING_TABLE


def _total_width_table():
    global _TOTAL_WIDTH_TABLE
    if _TOTAL_WIDTH_TABLE is None:
        metadata = json.loads((DATA_DIR / "widths_metadata.json").read_text())
        total = next(
            entry for entry in metadata["columns"]
            if entry["canonical_name"] == "total"
        )
        table = np.loadtxt(DATA_DIR / "widths_bnt.csv", delimiter=",", skiprows=1)
        _TOTAL_WIDTH_TABLE = (table[:, 0], table[:, total["source_index"] - 1])
    return _TOTAL_WIDTH_TABLE


def exclusive_branching_weights(mass_gev: float) -> dict[str, float]:
    """Unique exclusive BRs plus a conservative neutral missing-width mode."""
    table = _branching_table()
    masses = table["mass_GeV"]
    if not masses[0] <= mass_gev <= masses[-1]:
        raise ValueError(f"mass {mass_gev:g} GeV is outside the decay table")
    weights = {
        channel_id: max(float(np.interp(mass_gev, masses, table[channel_id])), 0.0)
        for channel_id in DECAY_PRODUCTS_PDG
        if channel_id not in DUPLICATE_CHANNEL_IDS
    }
    known = sum(weights.values())
    if known > 1.0:
        weights = {key: value / known for key, value in weights.items()}
        known = 1.0
    weights[UNCLASSIFIED_CHANNEL_ID] = max(1.0 - known, 0.0)
    return {key: value for key, value in weights.items() if value > 0.0}


def ctau_at_reference_coupling(mass_gev: float) -> float:
    masses, coefficients = _total_width_table()
    if not masses[0] <= mass_gev <= masses[-1]:
        raise ValueError(f"mass {mass_gev:g} GeV is outside the width table")
    width = float(np.interp(mass_gev, masses, coefficients)) * INV_F_REF**2
    return HBAR_C_GEV_M / width if width > 0.0 else np.inf

