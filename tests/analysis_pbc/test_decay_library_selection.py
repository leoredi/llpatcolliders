#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys

import pytest

HERE = Path(__file__).resolve()
REPO_ROOT = HERE.parents[2]
ANALYSIS_ROOT = REPO_ROOT / "analysis_pbc"
if str(ANALYSIS_ROOT) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_ROOT))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from decay import rhn_decay_library as lib
from tools.decay.decay_library_io import write_decay_library_txt


def _decay_dir(root: Path, flavour: str) -> Path:
    cfg = lib.FLAVOUR_CONFIG[flavour]
    return root / cfg["repo"] / cfg["decay_dir"]


def _write_simple_decay_file(path: Path, mass: float) -> None:
    events = [[(0.1, 0.0, 0.0, 0.0, 0.000511, 11)]]
    write_decay_library_txt(path, mass_GeV=mass, event_daughters=events)


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


def test_overlay_precedence_on_same_basename(roots):
    external, generated = roots
    flavour = "electron"
    base = "vN_Ntoall_generated_6.0.txt"

    ext_path = _decay_dir(external, flavour) / base
    gen_path = _decay_dir(generated, flavour) / base
    _write_simple_decay_file(ext_path, 6.0)
    _write_simple_decay_file(gen_path, 6.0)

    chosen = lib.select_decay_file(flavour, 6.0)
    assert chosen.source == "generated"
    assert chosen.path == gen_path


def test_below_switch_prefers_external_with_legacy_categories(roots):
    external, generated = roots
    flavour = "electron"
    _write_simple_decay_file(_decay_dir(generated, flavour) / "vN_Ntoall_generated_3.0.txt", 3.0)
    _write_simple_decay_file(_decay_dir(external, flavour) / "vN_Ntoall_inclDs_3.01.txt", 3.01)

    chosen = lib.select_decay_file(flavour, 3.0)
    assert chosen.source == "external"
    assert chosen.path.name == "vN_Ntoall_inclDs_3.01.txt"


def test_above_switch_uses_generated_only(roots):
    external, generated = roots
    flavour = "electron"
    _write_simple_decay_file(_decay_dir(external, flavour) / "vN_Ntoall_inclDs_5.0.txt", 5.0)
    _write_simple_decay_file(_decay_dir(generated, flavour) / "vN_Ntoall_generated_5.0.txt", 5.0)

    chosen = lib.select_decay_file(flavour, 5.0)
    assert chosen.source == "generated"
    assert chosen.path.name == "vN_Ntoall_generated_5.0.txt"


def test_above_switch_missing_overlay_raises(roots):
    external, _ = roots
    flavour = "electron"
    _write_simple_decay_file(_decay_dir(external, flavour) / "vN_Ntoall_inclDs_6.0.txt", 6.0)

    with pytest.raises(FileNotFoundError, match="No generated decay overlay files found"):
        lib.select_decay_file(flavour, 6.0)


def test_below_switch_missing_external_warns_and_falls_back_to_generated(roots):
    _, generated = roots
    flavour = "electron"
    _write_simple_decay_file(_decay_dir(generated, flavour) / "vN_Ntoall_generated_3.0.txt", 3.0)

    with pytest.warns(UserWarning, match="falling back to generated source"):
        chosen = lib.select_decay_file(flavour, 3.0)
    assert chosen.source == "generated"


def test_deterministic_tiebreak_for_nearest_overlay(roots):
    _, generated = roots
    flavour = "electron"
    _write_simple_decay_file(_decay_dir(generated, flavour) / "vN_Ntoall_generated_6.0.txt", 6.0)
    _write_simple_decay_file(_decay_dir(generated, flavour) / "vN_Ntoall_generated_6.2.txt", 6.2)

    chosen = lib.select_decay_file(flavour, 6.1)
    assert chosen.mass_GeV == pytest.approx(6.0)


def test_low_mass_prefers_analytical_files(roots):
    external, generated = roots
    flavour = "electron"

    _write_simple_decay_file(_decay_dir(generated, flavour) / "vN_Ntoall_generated_0.30.txt", 0.30)
    _write_simple_decay_file(
        _decay_dir(external, flavour) / "vN_Ntoall_analytical2and3bodydecays_0.31.txt", 0.31
    )
    _write_simple_decay_file(_decay_dir(external, flavour) / "vN_Ntoall_lightfonly_0.30.txt", 0.30)

    chosen = lib.select_decay_file(flavour, 0.30)
    assert chosen.category == "analytical2and3bodydecays"
    assert chosen.source == "external"


def test_strict_mismatch_raises_by_default(roots):
    external, _ = roots
    flavour = "electron"
    _write_simple_decay_file(_decay_dir(external, flavour) / "vN_Ntoall_inclDs_1.0.txt", 1.0)

    with pytest.raises(ValueError, match="Refusing to extrapolate"):
        lib.select_decay_file(flavour, 2.0)


def test_mismatch_override_env_downgrades_to_warning(roots, monkeypatch: pytest.MonkeyPatch):
    external, _ = roots
    flavour = "electron"
    _write_simple_decay_file(_decay_dir(external, flavour) / "vN_Ntoall_inclDs_1.0.txt", 1.0)
    monkeypatch.setenv("HNL_ALLOW_DECAY_MASS_MISMATCH", "1")

    with pytest.warns(UserWarning, match="allowed by HNL_ALLOW_DECAY_MASS_MISMATCH"):
        chosen = lib.select_decay_file(flavour, 2.0)
    assert chosen.mass_GeV == pytest.approx(1.0)


def test_generated_txt_parser_compatibility(tmp_path: Path):
    path = tmp_path / "vN_Ntoall_generated_5.0.txt"
    event_daughters = [
        [(1.0, 0.1, 0.2, 0.3, 0.13957, 211), (2.0, -0.1, -0.2, -0.3, 0.13957, -211)],
        [(3.0, 0.0, 0.0, 1.0, 0.0, 14)],
    ]
    write_decay_library_txt(path, mass_GeV=5.0, event_daughters=event_daughters)
    parsed = lib.load_decay_events(path)
    assert len(parsed) == 2
    assert len(parsed[0]) == 2
    assert len(parsed[1]) == 1
