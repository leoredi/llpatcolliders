#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path
import sys

import pytest

HERE = Path(__file__).resolve()
REPO_ROOT = HERE.parents[3]
ANALYSIS_ROOT = REPO_ROOT / "analysis_pbc"
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from decay import rhn_decay_library as lib
from tools.decay.decay_library_io import write_decay_library_txt
from tools.decay import validate_decay_overlap as overlap


def _decay_dir(root: Path, flavour: str) -> Path:
    cfg = lib.FLAVOUR_CONFIG[flavour]
    return root / cfg["repo"] / cfg["decay_dir"]


@pytest.fixture(autouse=True)
def _clear_caches():
    lib.list_decay_files.cache_clear()
    lib.load_decay_events.cache_clear()
    yield
    lib.list_decay_files.cache_clear()
    lib.load_decay_events.cache_clear()


@pytest.fixture
def roots(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    external = tmp_path / "external"
    generated = tmp_path / "generated"
    monkeypatch.setattr(lib, "EXTERNAL_ROOT", external)
    monkeypatch.setattr(lib, "GENERATED_ROOT", generated)
    return external, generated


def test_nearest_for_source_uses_requested_source(roots):
    external, generated = roots
    flavour = "electron"
    ext_path = _decay_dir(external, flavour) / "vN_Ntoall_inclDs_4.8.txt"
    gen_path = _decay_dir(generated, flavour) / "vN_Ntoall_generated_4.1.txt"

    write_decay_library_txt(ext_path, mass_GeV=4.8, event_daughters=[[(1.0, 0.0, 0.0, 0.0, 0.0, 11)]])
    write_decay_library_txt(gen_path, mass_GeV=4.1, event_daughters=[[(1.0, 0.0, 0.0, 0.0, 0.0, 11)]])

    ext = overlap._nearest_for_source(flavour, "external", 4.5)
    gen = overlap._nearest_for_source(flavour, "generated", 4.5)
    assert ext is not None
    assert gen is not None
    assert ext.source == "external"
    assert gen.source == "generated"


def test_mean_metrics_computes_visible_fraction_and_daughter_mean():
    events = [
        [(2.0, 0.0, 0.0, 0.0, 0.0, 11), (1.0, 0.0, 0.0, 0.0, 0.0, 12)],
        [(3.0, 0.0, 0.0, 0.0, 0.0, 211)],
    ]
    n_events, mean_daughters, mean_visible_frac = overlap._mean_metrics(events)

    assert n_events == 2
    assert mean_daughters == pytest.approx(1.5)
    # Event 1 visible fraction: 2/3, event 2: 1.0, mean: 5/6.
    assert mean_visible_frac == pytest.approx(5.0 / 6.0)


def test_check_delta_fails_without_override(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.delenv("HNL_ALLOW_DECAY_MASS_MISMATCH", raising=False)
    entry = lib.DecayFileEntry(
        path=Path("/tmp/fake.txt"),
        mass_GeV=3.0,
        category="unknown",
        source="generated",
        source_priority=0,
    )
    ok, delta = overlap._check_delta("generated", "electron", 4.0, entry)
    assert not ok
    assert delta == pytest.approx(1.0)


def test_check_delta_warns_with_override(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setenv("HNL_ALLOW_DECAY_MASS_MISMATCH", "1")
    entry = lib.DecayFileEntry(
        path=Path("/tmp/fake.txt"),
        mass_GeV=3.0,
        category="unknown",
        source="generated",
        source_priority=0,
    )
    with pytest.warns(UserWarning, match="allowed by HNL_ALLOW_DECAY_MASS_MISMATCH"):
        ok, delta = overlap._check_delta("generated", "electron", 4.0, entry)
    assert ok
    assert delta == pytest.approx(1.0)


def test_relative_diff_uses_reference_denominator():
    assert overlap._relative_diff(10.0, 8.0) == pytest.approx(0.25)
    assert overlap._relative_diff(0.0, 0.0) == pytest.approx(0.0)


def test_overlap_main_identical_external_and_generated_passes(
    roots, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    external, generated = roots
    flavour = "electron"
    mass = 4.2
    event_daughters = [[(2.0, 0.0, 0.0, 0.0, 0.0, 11), (1.0, 0.0, 0.0, 0.0, 0.0, 12)]]

    write_decay_library_txt(
        _decay_dir(external, flavour) / "vN_Ntoall_inclDs_4.2.txt",
        mass_GeV=mass,
        event_daughters=event_daughters,
    )
    write_decay_library_txt(
        _decay_dir(generated, flavour) / "vN_Ntoall_generated_4.2.txt",
        mass_GeV=mass,
        event_daughters=event_daughters,
    )

    out_csv = tmp_path / "overlap_ok.csv"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "validate_decay_overlap.py",
            "--flavours",
            flavour,
            "--masses",
            "4.2",
            "--min-mass",
            "4.0",
            "--max-mass",
            "5.0",
            "--out",
            str(out_csv),
        ],
    )
    overlap.main()

    rows = list(csv.DictReader(out_csv.open()))
    assert len(rows) == 1
    assert rows[0]["status"] == "ok"


def test_overlap_main_visible_fraction_shift_fails(
    roots, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    external, generated = roots
    flavour = "electron"
    mass = 4.2

    # External has one invisible daughter -> visible fraction = 0.5.
    write_decay_library_txt(
        _decay_dir(external, flavour) / "vN_Ntoall_inclDs_4.2.txt",
        mass_GeV=mass,
        event_daughters=[[(1.0, 0.0, 0.0, 0.0, 0.0, 11), (1.0, 0.0, 0.0, 0.0, 0.0, 12)]],
    )
    # Generated is fully visible -> visible fraction = 1.0.
    write_decay_library_txt(
        _decay_dir(generated, flavour) / "vN_Ntoall_generated_4.2.txt",
        mass_GeV=mass,
        event_daughters=[[(2.0, 0.0, 0.0, 0.0, 0.0, 11)]],
    )

    out_csv = tmp_path / "overlap_fail.csv"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "validate_decay_overlap.py",
            "--flavours",
            flavour,
            "--masses",
            "4.2",
            "--min-mass",
            "4.0",
            "--max-mass",
            "5.0",
            "--out",
            str(out_csv),
        ],
    )
    with pytest.raises(SystemExit, match="1"):
        overlap.main()

    rows = list(csv.DictReader(out_csv.open()))
    assert len(rows) == 1
    assert "visible_frac_fail" in rows[0]["status"]
