from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import csv
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

# A normalization key matches how expected_signal.py groups contributions:
# - ("pid", abs(parent_pdg)) for direct parents.
# - ("tau_parent", abs(tau_parent_id)) for fromTau chain rows (parent_pdg==15).
NormKey = tuple[str, int]

ALL_REGIMES = {"all", "combined"}
KAON_PARENT_PDGS = {130, 321}
CHARM_PARENT_PDGS = {411, 421, 431, 4122}
BEAUTY_PARENT_PDGS = {511, 521, 531, 541, 5122, 5232, 5332}
EW_PARENT_PDGS = {23, 24}


@dataclass(frozen=True)
class OverlapSample:
    base_regime: str
    mode: str | None
    is_ff: bool
    qcd_mode: str
    pthat_min: float | None
    path: Path


@dataclass(frozen=True)
class ResolvedOverlapSample:
    sample: OverlapSample
    owned_keys: tuple[NormKey, ...]
    owned_events: int
    total_events: int


@dataclass(frozen=True)
class OverlapResolution:
    samples: tuple[ResolvedOverlapSample, ...]
    warnings: tuple[str, ...]
    total_owned_events: int


def _parse_abs_int(raw: object) -> int:
    if raw is None:
        return 0
    text = str(raw).strip()
    if not text or text.lower() in {"nan", "none"}:
        return 0
    try:
        return abs(int(float(text)))
    except Exception:
        return 0


def _norm_key(parent_abs: int, tau_parent_abs: int) -> NormKey:
    if parent_abs == 15 and tau_parent_abs > 0:
        return ("tau_parent", tau_parent_abs)
    return ("pid", parent_abs)


def format_norm_key(key: NormKey) -> str:
    kind, pdg = key
    if kind == "pid":
        return f"pid:{pdg}"
    if kind == "tau_parent":
        return f"tau_parent:{pdg}"
    return f"{kind}:{pdg}"


def _variant_priority(sample: OverlapSample) -> tuple[int, int, float]:
    if sample.base_regime == "charm" and sample.qcd_mode == "hardccbar":
        qcd_priority = 3
    elif sample.base_regime in {"beauty", "Bc"} and sample.qcd_mode in {"hardbbbar", "hardBc"}:
        qcd_priority = 3
    elif sample.qcd_mode != "auto":
        qcd_priority = 2
    else:
        qcd_priority = 1
    ff_priority = 1 if sample.is_ff else 0
    pthat_priority = float(sample.pthat_min) if sample.pthat_min is not None else -1.0
    return (qcd_priority, ff_priority, pthat_priority)


def _key_sector(key: NormKey) -> str:
    _, pdg = key
    if pdg in EW_PARENT_PDGS:
        return "ew"
    if pdg in KAON_PARENT_PDGS:
        return "kaon"
    if pdg in CHARM_PARENT_PDGS:
        return "charm"
    if pdg in BEAUTY_PARENT_PDGS:
        return "beauty"
    return "other"


def _regime_match_score(base_regime: str, sector: str) -> int:
    if base_regime in ALL_REGIMES:
        return 6

    if sector == "ew":
        return 6 if base_regime == "ew" else 0
    if sector == "kaon":
        return {"kaon": 6, "charm": 2}.get(base_regime, 0)
    if sector == "charm":
        return {"charm": 6, "beauty": 1, "kaon": 1}.get(base_regime, 0)
    if sector == "beauty":
        return {"beauty": 6, "Bc": 6, "charm": 1}.get(base_regime, 0)
    return 1


@lru_cache(maxsize=4096)
def _count_norm_keys_cached(path_text: str, mtime_ns: int, size_bytes: int) -> tuple[tuple[NormKey, int], ...]:
    del mtime_ns, size_bytes  # cache key invalidation only
    counts: dict[NormKey, int] = defaultdict(int)
    path = Path(path_text)
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            parent_abs = _parse_abs_int(row.get("parent_pdg", row.get("parent_id")))
            if parent_abs <= 0:
                continue
            tau_parent_abs = _parse_abs_int(row.get("tau_parent_id"))
            key = _norm_key(parent_abs, tau_parent_abs)
            counts[key] += 1
    return tuple(sorted(counts.items(), key=lambda kv: (kv[0][0], kv[0][1])))


def count_normalization_keys(path: Path) -> dict[NormKey, int]:
    stat = path.stat()
    items = _count_norm_keys_cached(str(path), stat.st_mtime_ns, stat.st_size)
    return dict(items)


def resolve_parent_overlap(
    samples: Iterable[OverlapSample],
    *,
    context: str = "",
    min_events_per_mass: int = 0,
    strict_min_events: bool = False,
) -> OverlapResolution:
    sample_list = list(samples)
    if not sample_list:
        return OverlapResolution(samples=tuple(), warnings=tuple(), total_owned_events=0)

    key_counts: dict[OverlapSample, dict[NormKey, int]] = {}
    total_counts: dict[OverlapSample, int] = {}
    owners_for_key: dict[NormKey, OverlapSample] = {}
    warnings: list[str] = []

    for sample in sample_list:
        counts = count_normalization_keys(sample.path)
        key_counts[sample] = counts
        total_counts[sample] = int(sum(counts.values()))

    all_keys: set[NormKey] = set()
    for counts in key_counts.values():
        all_keys.update(counts.keys())

    for key in sorted(all_keys, key=lambda k: (k[0], k[1])):
        contenders = [s for s in sample_list if key_counts[s].get(key, 0) > 0]
        if not contenders:
            continue
        if len(contenders) == 1:
            owners_for_key[key] = contenders[0]
            continue

        sector = _key_sector(key)

        def _score(sample: OverlapSample) -> tuple[int, int, tuple[int, int, float], int, str]:
            return (
                _regime_match_score(sample.base_regime, sector),
                key_counts[sample][key],
                _variant_priority(sample),
                total_counts[sample],
                str(sample.path),
            )

        winner = max(contenders, key=_score)
        owners_for_key[key] = winner

        contenders_text = ", ".join(
            f"{s.path.name}[{s.base_regime}{'_' + s.mode if s.mode else ''}]"
            for s in sorted(contenders, key=lambda s: str(s.path))
        )
        warnings.append(
            f"[OVERLAP]{' ' + context if context else ''} {format_norm_key(key)} "
            f"appears in {len(contenders)} files; keeping {winner.path.name} "
            f"(owner regime={winner.base_regime}), dropping duplicates from: {contenders_text}"
        )

    owned_by_sample: dict[OverlapSample, set[NormKey]] = {s: set() for s in sample_list}
    for key, owner in owners_for_key.items():
        owned_by_sample[owner].add(key)

    resolved: list[ResolvedOverlapSample] = []
    total_owned_events = 0
    for sample in sample_list:
        owned_keys = tuple(sorted(owned_by_sample[sample], key=lambda k: (k[0], k[1])))
        if not owned_keys:
            continue
        owned_events = int(sum(key_counts[sample].get(k, 0) for k in owned_keys))
        total_events = int(total_counts[sample])
        total_owned_events += owned_events
        resolved.append(
            ResolvedOverlapSample(
                sample=sample,
                owned_keys=owned_keys,
                owned_events=owned_events,
                total_events=total_events,
            )
        )

    if min_events_per_mass > 0 and total_owned_events < int(min_events_per_mass):
        msg = (
            f"[LOW-STAT]{' ' + context if context else ''} total owned events={total_owned_events} "
            f"< min_events_per_mass={int(min_events_per_mass)} after overlap resolution."
        )
        if strict_min_events:
            raise ValueError(msg)
        warnings.append(msg)

    return OverlapResolution(
        samples=tuple(resolved),
        warnings=tuple(warnings),
        total_owned_events=total_owned_events,
    )


def filter_dataframe_by_norm_keys(df: pd.DataFrame, owned_keys: Iterable[NormKey]) -> pd.DataFrame:
    keys = tuple(owned_keys)
    if len(keys) == 0:
        return df.iloc[0:0].copy()

    if "parent_pdg" in df.columns:
        parent_abs = pd.to_numeric(df["parent_pdg"], errors="coerce").fillna(0.0).abs().astype(int).to_numpy()
    elif "parent_id" in df.columns:
        parent_abs = pd.to_numeric(df["parent_id"], errors="coerce").fillna(0.0).abs().astype(int).to_numpy()
    else:
        return df.iloc[0:0].copy()

    if "tau_parent_id" in df.columns:
        tau_parent_abs = (
            pd.to_numeric(df["tau_parent_id"], errors="coerce")
            .fillna(0.0)
            .abs()
            .astype(int)
            .to_numpy()
        )
    else:
        tau_parent_abs = np.zeros(len(df), dtype=int)

    mask_tau_chain = (parent_abs == 15) & (tau_parent_abs > 0)
    keep = np.zeros(len(df), dtype=bool)
    for kind, pdg in keys:
        if kind == "tau_parent":
            keep |= mask_tau_chain & (tau_parent_abs == int(pdg))
            continue
        if kind != "pid":
            raise ValueError(f"Unsupported normalization key kind: {kind}")
        if int(pdg) == 15:
            keep |= (parent_abs == 15) & (~mask_tau_chain)
        else:
            keep |= (~mask_tau_chain) & (parent_abs == int(pdg))

    return df.loc[keep].copy()
