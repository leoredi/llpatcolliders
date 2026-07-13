"""Tests for the staged SensCalc 2501 export boundary."""

import json

import pytest

from alp_fermion.tools import export_senscalc_2501 as export


def _write_complete_export(path):
    path.mkdir()
    hashes = {}
    for name in export.EXPECTED_OUTPUT_FILES:
        if name == "EXPORT_MANIFEST.json":
            continue
        product = path / name
        product.write_text(f"test product: {name}\n")
        hashes[name] = export.sha256(product)
    (path / "EXPORT_MANIFEST.json").write_text(json.dumps({
        "senscalc_commit": export.SENSCALC_COMMIT,
        "output_sha256": hashes,
    }))


def test_validate_export_accepts_complete_hashed_bundle(tmp_path):
    staged = tmp_path / "staged"
    _write_complete_export(staged)
    export.validate_export(staged)


def test_validate_export_rejects_modified_product(tmp_path):
    staged = tmp_path / "staged"
    _write_complete_export(staged)
    (staged / "widths_bnt.csv").write_text("modified\n")
    with pytest.raises(export.InputError, match="hash mismatch"):
        export.validate_export(staged)


def test_install_export_replaces_destination_bundle(tmp_path):
    staged = tmp_path / "staged"
    destination = tmp_path / "published"
    _write_complete_export(staged)
    destination.mkdir()
    (destination / "widths_bnt.csv").write_text("stale\n")

    export.install_export(staged, destination)

    assert (destination / "widths_bnt.csv").read_text().startswith("test product")
    assert not any(staged.iterdir())
