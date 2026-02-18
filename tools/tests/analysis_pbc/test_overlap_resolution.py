#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
import pytest

from analysis_pbc.limits.overlap_resolution import (
    OverlapSample,
    filter_dataframe_by_norm_keys,
    format_norm_key,
    resolve_parent_overlap,
)


def _write_minimal_sim_csv(path: Path, rows: list[tuple[int, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["event", "weight", "parent_pdg", "tau_parent_id"])
        for i, (parent_pdg, tau_parent_id) in enumerate(rows):
            writer.writerow([i, 1.0, parent_pdg, tau_parent_id])


def _sample(path: Path, regime: str, mode: str | None = None) -> OverlapSample:
    return OverlapSample(
        base_regime=regime,
        mode=mode,
        is_ff=False,
        qcd_mode="auto",
        pthat_min=None,
        path=path,
    )


def test_resolve_parent_overlap_prefers_parent_owner(tmp_path: Path) -> None:
    charm_path = tmp_path / "HNL_1p98GeV_electron_charm.csv"
    beauty_path = tmp_path / "HNL_1p98GeV_electron_beauty.csv"
    _write_minimal_sim_csv(charm_path, [(411, 0)] * 5 + [(511, 0)] * 2)
    _write_minimal_sim_csv(beauty_path, [(511, 0)] * 10)

    overlap = resolve_parent_overlap(
        [_sample(charm_path, "charm"), _sample(beauty_path, "beauty")],
        context="m=1.98 (electron)",
    )

    owner_by_key = {}
    for resolved in overlap.samples:
        for key in resolved.owned_keys:
            owner_by_key[key] = resolved.sample.base_regime

    assert owner_by_key[("pid", 411)] == "charm"
    assert owner_by_key[("pid", 511)] == "beauty"
    assert overlap.total_owned_events == 15
    assert any("pid:511" in msg for msg in overlap.warnings)


def test_resolve_parent_overlap_handles_tau_chain_keys(tmp_path: Path) -> None:
    charm_tau = tmp_path / "HNL_1p90GeV_tau_charm_fromTau.csv"
    beauty_tau = tmp_path / "HNL_1p90GeV_tau_beauty_fromTau.csv"
    _write_minimal_sim_csv(charm_tau, [(15, 431)] * 4)
    _write_minimal_sim_csv(beauty_tau, [(15, 531)] * 6 + [(15, 431)] * 1)

    overlap = resolve_parent_overlap(
        [_sample(charm_tau, "charm", "fromTau"), _sample(beauty_tau, "beauty", "fromTau")],
        context="m=1.90 (tau)",
    )

    owner_by_key = {}
    for resolved in overlap.samples:
        for key in resolved.owned_keys:
            owner_by_key[key] = resolved.sample.base_regime

    assert owner_by_key[("tau_parent", 431)] == "charm"
    assert owner_by_key[("tau_parent", 531)] == "beauty"


def test_filter_dataframe_by_norm_keys_respects_tau_chain() -> None:
    df = pd.DataFrame(
        [
            {"parent_pdg": 15, "tau_parent_id": 431},
            {"parent_pdg": 15, "tau_parent_id": 531},
            {"parent_pdg": 511, "tau_parent_id": 0},
            {"parent_pdg": 521, "tau_parent_id": 0},
        ]
    )
    kept = filter_dataframe_by_norm_keys(df, [("tau_parent", 531), ("pid", 511)])
    assert len(kept) == 2
    keys = []
    for _, row in kept.iterrows():
        parent = abs(int(row["parent_pdg"]))
        tau_parent = abs(int(row["tau_parent_id"]))
        key = ("tau_parent", tau_parent) if (parent == 15 and tau_parent > 0) else ("pid", parent)
        keys.append(format_norm_key(key))
    assert set(keys) == {"tau_parent:531", "pid:511"}


def test_resolve_parent_overlap_min_events_strict(tmp_path: Path) -> None:
    low_stat = tmp_path / "HNL_4p40GeV_tau_beauty_direct.csv"
    _write_minimal_sim_csv(low_stat, [(541, 0)] * 3)
    with pytest.raises(ValueError, match="LOW-STAT"):
        resolve_parent_overlap(
            [_sample(low_stat, "beauty", "direct")],
            context="m=4.40 (tau)",
            min_events_per_mass=10,
            strict_min_events=True,
        )
