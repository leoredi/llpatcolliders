#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve()
REPO_ROOT = HERE.parents[3]
ANALYSIS_ROOT = REPO_ROOT / "analysis_pbc"
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from decay import brvis_kappa
from decay import rhn_decay_library as lib
from limits import expected_signal as es
from tools.decay.decay_library_io import write_decay_library_txt


def _geom_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "parent_id": [511, 511],
            "weight": [1.0, 1.0],
            "beta_gamma": [2.0, 2.0],
            "hits_tube": [True, False],
            "entry_distance": [10.0, 0.0],
            "path_length": [2.0, 0.0],
            "eta": [0.1, 0.2],
            "phi": [0.0, 0.1],
        }
    )


def _decay_dir(root: Path, flavour: str) -> Path:
    cfg = lib.FLAVOUR_CONFIG[flavour]
    return root / cfg["repo"] / cfg["decay_dir"]


def _write_kappa_table(path: Path, rows: list[dict]) -> None:
    df = pd.DataFrame(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


@pytest.fixture(autouse=True)
def _clear_caches():
    lib.list_decay_files.cache_clear()
    lib.load_decay_events.cache_clear()
    brvis_kappa._load_kappa_table_cached.cache_clear()
    yield
    lib.list_decay_files.cache_clear()
    lib.load_decay_events.cache_clear()
    brvis_kappa._load_kappa_table_cached.cache_clear()


@pytest.fixture
def roots(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    external = tmp_path / "external"
    generated = tmp_path / "generated"
    monkeypatch.setattr(lib, "EXTERNAL_ROOT", external)
    monkeypatch.setattr(lib, "GENERATED_ROOT", generated)
    return external, generated


def test_brvis_kappa_seed_invariance():
    kwargs = dict(
        geom_df=_geom_df(),
        mass_GeV=4.0,
        eps2=1e-6,
        benchmark="100",
        lumi_fb=1.0,
        separation_m=1e-3,
        decay_mode="brvis_kappa",
        br_vis=0.8,
        kappa_eff=0.5,
        ctau0_m=1.0,
        br_per_parent={511: 0.1},
        br_scale=1.0,
    )
    n_seed_a = es.expected_signal_events(**kwargs, decay_seed=12345)
    n_seed_b = es.expected_signal_events(**kwargs, decay_seed=67890)
    assert n_seed_a == pytest.approx(n_seed_b, rel=0.0, abs=0.0)


def test_library_seed_sensitivity_sanity(monkeypatch: pytest.MonkeyPatch):
    def fake_acceptance(*, geom_df, selection, **kwargs):
        out = np.zeros(len(geom_df), dtype=bool)
        out[0] = (selection.seed % 2) == 0
        return out

    monkeypatch.setattr(es, "compute_decay_acceptance", fake_acceptance)
    monkeypatch.setattr(es, "build_mesh_once", lambda *args, **kwargs: object())

    kwargs = dict(
        geom_df=_geom_df(),
        mass_GeV=4.0,
        eps2=1e-6,
        benchmark="100",
        lumi_fb=1.0,
        separation_m=1e-3,
        decay_mode="library",
        ctau0_m=1.0,
        br_per_parent={511: 0.1},
        br_scale=1.0,
    )

    n_odd = es.expected_signal_events(**kwargs, decay_seed=12345)
    n_even = es.expected_signal_events(**kwargs, decay_seed=12346)
    assert n_odd != n_even
    assert n_even > n_odd


def test_boundary_routing_contract_at_five_gev(roots):
    external, generated = roots
    flavour = "electron"

    write_decay_library_txt(
        _decay_dir(external, flavour) / "vN_Ntoall_inclDs_4.999.txt",
        mass_GeV=4.999,
        event_daughters=[[(1.0, 0.0, 0.0, 0.0, 0.0, 11)]],
    )
    write_decay_library_txt(
        _decay_dir(generated, flavour) / "vN_Ntoall_generated_5.0.txt",
        mass_GeV=5.0,
        event_daughters=[[(1.0, 0.0, 0.0, 0.0, 0.0, 11)]],
    )

    below = lib.select_decay_file(flavour, 4.999)
    at_switch = lib.select_decay_file(flavour, 5.0)
    assert below.source == "external"
    assert at_switch.source == "generated"


def test_strict_mismatch_contract_threshold_and_override(
    roots, monkeypatch: pytest.MonkeyPatch
):
    external, _ = roots
    flavour = "electron"
    write_decay_library_txt(
        _decay_dir(external, flavour) / "vN_Ntoall_inclDs_1.0.txt",
        mass_GeV=1.0,
        event_daughters=[[(1.0, 0.0, 0.0, 0.0, 0.0, 11)]],
    )

    chosen = lib.select_decay_file(flavour, 1.5)
    assert chosen.mass_GeV == pytest.approx(1.0)

    with pytest.raises(ValueError, match="Refusing to extrapolate"):
        lib.select_decay_file(flavour, 1.5001)

    monkeypatch.setenv("HNL_ALLOW_DECAY_MASS_MISMATCH", "1")
    with pytest.warns(UserWarning, match="allowed by HNL_ALLOW_DECAY_MASS_MISMATCH"):
        chosen_override = lib.select_decay_file(flavour, 1.5001)
    assert chosen_override.mass_GeV == pytest.approx(1.0)


def test_parser_neutrino_preservation_for_decimal_pid_tokens(tmp_path: Path):
    path = tmp_path / "decimal_pid.txt"
    path.write_text(
        "\n".join(
            [
                "Format is groups of ...",
                "10.0,0.0,0.0,0.0,5.0,9900012,N",
                "1.0,0.0,0.0,1.0,0.0,16.0,vt",
                "1.0,0.0,0.0,-1.0,0.0,-16.0,vt~",
                "2.0,0.0,0.0,0.0,0.0,11.0,e-",
                "",
                "",
            ]
        )
    )

    events = lib.load_decay_events(path)
    pids = [int(d[5]) for d in events[0]]
    assert pids.count(16) == 1
    assert pids.count(-16) == 1
    assert pids.count(11) == 1

    total_e = sum(float(d[0]) for d in events[0])
    visible_e = sum(float(d[0]) for d in events[0] if abs(int(d[5])) not in {12, 14, 16})
    assert total_e == pytest.approx(4.0)
    assert visible_e / total_e == pytest.approx(0.5)


def test_no_double_counting_contract_library_vs_brvis_kappa():
    geom_df = _geom_df()

    n_lib_a = es.expected_signal_events(
        geom_df=geom_df,
        mass_GeV=4.0,
        eps2=1e-6,
        benchmark="100",
        lumi_fb=1.0,
        separation_m=1e-3,
        decay_mode="library",
        separation_pass=np.array([True, False], dtype=bool),
        br_vis=0.2,
        kappa_eff=0.3,
        ctau0_m=1.0,
        br_per_parent={511: 0.1},
        br_scale=1.0,
    )
    n_lib_b = es.expected_signal_events(
        geom_df=geom_df,
        mass_GeV=4.0,
        eps2=1e-6,
        benchmark="100",
        lumi_fb=1.0,
        separation_m=1e-3,
        decay_mode="library",
        separation_pass=np.array([True, False], dtype=bool),
        br_vis=0.95,
        kappa_eff=1.7,
        ctau0_m=1.0,
        br_per_parent={511: 0.1},
        br_scale=1.0,
    )
    assert n_lib_a == pytest.approx(n_lib_b, rel=0.0, abs=0.0)

    n_brvis_1 = es.expected_signal_events(
        geom_df=geom_df,
        mass_GeV=4.0,
        eps2=1e-6,
        benchmark="100",
        lumi_fb=1.0,
        separation_m=1e-3,
        decay_mode="brvis_kappa",
        br_vis=0.6,
        kappa_eff=0.5,
        ctau0_m=1.0,
        br_per_parent={511: 0.1},
        br_scale=1.0,
    )
    n_brvis_2 = es.expected_signal_events(
        geom_df=geom_df,
        mass_GeV=4.0,
        eps2=1e-6,
        benchmark="100",
        lumi_fb=1.0,
        separation_m=1e-3,
        decay_mode="brvis_kappa",
        br_vis=0.6,
        kappa_eff=1.0,
        ctau0_m=1.0,
        br_per_parent={511: 0.1},
        br_scale=1.0,
    )
    assert n_brvis_2 == pytest.approx(2.0 * n_brvis_1, rel=1e-12, abs=0.0)


def test_kappa_metadata_enforcement_rejects_separation_mismatch(tmp_path: Path):
    path = tmp_path / "kappa.csv"
    _write_kappa_table(
        path,
        [
            {
                "flavour": "electron",
                "mass_GeV": 4.0,
                "kappa": 1.0,
                "p_min_GeV": 0.6,
                "separation_mm": 1.0,
                "source_policy": "hybrid",
                "status": "ok",
            }
        ],
    )

    with pytest.raises(brvis_kappa.KappaTableError, match="separation mismatch"):
        brvis_kappa.lookup_kappa("electron", 4.0, p_min_GeV=0.6, separation_mm=2.0, table_path=path)


def test_kappa_interpolation_monotonic_sanity(tmp_path: Path):
    path = tmp_path / "kappa.csv"
    _write_kappa_table(
        path,
        [
            {
                "flavour": "muon",
                "mass_GeV": 2.0,
                "kappa": 0.8,
                "p_min_GeV": 0.6,
                "separation_mm": 1.0,
                "source_policy": "hybrid",
                "status": "ok",
            },
            {
                "flavour": "muon",
                "mass_GeV": 4.0,
                "kappa": 1.0,
                "p_min_GeV": 0.6,
                "separation_mm": 1.0,
                "source_policy": "hybrid",
                "status": "ok",
            },
            {
                "flavour": "muon",
                "mass_GeV": 6.0,
                "kappa": 1.2,
                "p_min_GeV": 0.6,
                "separation_mm": 1.0,
                "source_policy": "hybrid",
                "status": "ok",
            },
        ],
    )

    kappa_3 = brvis_kappa.lookup_kappa("muon", 3.0, p_min_GeV=0.6, separation_mm=1.0, table_path=path)
    kappa_5 = brvis_kappa.lookup_kappa("muon", 5.0, p_min_GeV=0.6, separation_mm=1.0, table_path=path)

    assert 0.8 < kappa_3 < 1.0
    assert 1.0 < kappa_5 < 1.2
    assert kappa_5 > kappa_3
