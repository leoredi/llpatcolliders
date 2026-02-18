#!/usr/bin/env python3
"""Check that key code constants/modes are reflected in project docs.

This is a lightweight drift guard between:
- README.md (minimal operator entrypoint)
- CODING.md (implementation source)
- PHYSICS.md (physics source)
- AGENTS.md (agent execution contract)
"""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

THIS_FILE = Path(__file__).resolve()
REPO_ROOT = THIS_FILE.parents[2]

README = REPO_ROOT / "README.md"
CODING = REPO_ROOT / "CODING.md"
PHYSICS = REPO_ROOT / "PHYSICS.md"
AGENTS = REPO_ROOT / "AGENTS.md"

MASS_GRID_FILE = REPO_ROOT / "config_mass_grid.py"
XSECS_FILE = REPO_ROOT / "analysis_pbc" / "config" / "production_xsecs.py"
RUN_FILE = REPO_ROOT / "analysis_pbc" / "limits" / "run.py"
EXPECTED_SIGNAL_FILE = REPO_ROOT / "analysis_pbc" / "limits" / "expected_signal.py"
PARALLEL_FILE = REPO_ROOT / "production" / "pythia_production" / "run_parallel_production.sh"
MAIN_PYTHIA_FILE = REPO_ROOT / "production" / "pythia_production" / "main_hnl_production.cc"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _safe_eval_number(expr: str) -> float:
    node = ast.parse(expr, mode="eval")

    def _eval(n: ast.AST) -> float:
        if isinstance(n, ast.Expression):
            return _eval(n.body)
        if isinstance(n, ast.Constant) and isinstance(n.value, (int, float)):
            return float(n.value)
        if isinstance(n, ast.BinOp):
            left = _eval(n.left)
            right = _eval(n.right)
            if isinstance(n.op, ast.Add):
                return left + right
            if isinstance(n.op, ast.Sub):
                return left - right
            if isinstance(n.op, ast.Mult):
                return left * right
            if isinstance(n.op, ast.Div):
                return left / right
            if isinstance(n.op, ast.Pow):
                return left ** right
        if isinstance(n, ast.UnaryOp):
            val = _eval(n.operand)
            if isinstance(n.op, ast.UAdd):
                return val
            if isinstance(n.op, ast.USub):
                return -val
        raise RuntimeError(f"Unsupported numeric expression: {expr}")

    return float(_eval(node))


def _must_parse_assignment(text: str, var_name: str) -> float:
    match = re.search(rf"^{re.escape(var_name)}\s*=\s*(.+)$", text, re.M)
    if not match:
        raise RuntimeError(f"Could not parse assignment for {var_name}")
    rhs = match.group(1).split("#", 1)[0].strip()
    return _safe_eval_number(rhs)


def _must_match_float(pattern: str, text: str, label: str) -> float:
    match = re.search(pattern, text)
    if not match:
        raise RuntimeError(f"Could not parse {label} with pattern: {pattern}")
    return float(match.group(1))


def _parse_mass_grid(cfg_text: str) -> list[float]:
    block = re.search(r"MASS_GRID\s*=\s*sorted\(\[(.*?)\]\)", cfg_text, re.S)
    if not block:
        raise RuntimeError("Could not parse MASS_GRID block from config_mass_grid.py")
    vals = [float(tok) for tok in re.findall(r"\b\d+\.\d+\b", block.group(1))]
    if not vals:
        raise RuntimeError("Parsed empty MASS_GRID")
    return vals


def _contains(pattern: str, text: str) -> bool:
    return re.search(pattern, text, re.I | re.S) is not None


def main() -> int:
    failures: list[str] = []

    readme_text = _read(README)
    coding_text = _read(CODING)
    physics_text = _read(PHYSICS)
    agents_text = _read(AGENTS)

    mass_vals = _parse_mass_grid(_read(MASS_GRID_FILE))
    mass_count = len(mass_vals)
    mass_min = min(mass_vals)
    mass_max = max(mass_vals)

    xsec_text = _read(XSECS_FILE)
    sigma_cc_pb = _must_parse_assignment(xsec_text, "SIGMA_CCBAR_PB")
    sigma_bb_pb = _must_parse_assignment(xsec_text, "SIGMA_BBBAR_PB")
    sigma_bc_pb = _must_parse_assignment(xsec_text, "SIGMA_BC_PB")

    sigma_cc_mb = sigma_cc_pb / 1e9
    sigma_bb_ub = sigma_bb_pb / 1e6
    sigma_bc_ub = sigma_bc_pb / 1e6

    run_text = _read(RUN_FILE)
    lumi_fb = _must_match_float(r"L_HL_LHC_FB\s*=\s*([0-9.]+)", run_text, "L_HL_LHC_FB")

    required_run_flags = [
        "--separation-mm",
        "--max-separation-mm",
        "--separation-policy",
        "--geometry-model",
        "--detector-thickness-m",
        "--profile-inset-floor",
    ]
    for flag in required_run_flags:
        if flag not in run_text:
            failures.append(f"Code missing expected run.py flag: {flag}")

    expected_signal_text = _read(EXPECTED_SIGNAL_FILE)
    n_limit = _must_match_float(r"N_limit\s*:\s*float\s*=\s*([0-9.]+)", expected_signal_text, "N_limit")

    parallel_text = _read(PARALLEL_FILE)
    fromtau_threshold = _must_match_float(
        r"FROMTAU_MASS_THRESHOLD=([0-9.]+)",
        parallel_text,
        "FROMTAU_MASS_THRESHOLD",
    )

    main_pythia_text = _read(MAIN_PYTHIA_FILE)
    qcd_modes = ["auto", "hardccbar", "hardbbbar", "hardBc"]
    for mode in qcd_modes:
        if mode not in main_pythia_text:
            failures.append(f"Code does not contain expected QCD mode token: {mode}")

    docs_all = coding_text + "\n" + physics_text + "\n" + agents_text

    required_doc_checks = [
        ("mass-grid count", rf"\b{mass_count}\b"),
        ("mass-grid minimum", rf"\b{mass_min:.2f}\b"),
        ("mass-grid maximum", rf"\b{mass_max:.2f}\b"),
        ("tau fromTau threshold", rf"\b{fromtau_threshold:.2f}\s*GeV\b"),
        ("luminosity", rf"\b{int(lumi_fb)}\s*fb"),
        ("Poisson threshold", rf"\b{n_limit:.3f}\b"),
        ("sigma(ccbar)", rf"\b{sigma_cc_mb:g}\s*mb\b"),
        ("sigma(bbbar)", rf"\b{sigma_bb_ub:g}\s*(?:microbarn|μb|ub)\b"),
        ("sigma(Bc)", rf"\b{sigma_bc_ub:g}\s*(?:microbarn|μb|ub)\b"),
        ("FONLL/LHCb reference label", r"FONLL/LHCb"),
        ("Bc parent filter", r"541"),
        ("HNLCalc scaling script", r"check_hnlcalc_scaling\.py"),
        ("optional utility path (pythia monitor)", r"tools/pythia/monitor_production\.sh"),
        ("optional utility path (EW validator)", r"tools/madgraph/validate_xsec\.py"),
        ("optional utility path (decay generator)", r"tools/decay/generate_hnl_decay_events\.py"),
        ("run flag docs (--max-separation-mm)", r"--max-separation-mm"),
        ("run flag docs (--separation-policy)", r"--separation-policy"),
        ("run flag docs (--geometry-model)", r"--geometry-model"),
        ("separation policy docs (all-pairs-min)", r"all-pairs-min"),
        ("separation policy docs (any-pair-window)", r"any-pair-window"),
        ("geometry model docs (tube)", r"\btube\b"),
        ("geometry model docs (profile)", r"\bprofile\b"),
    ]

    for label, pattern in required_doc_checks:
        if not _contains(pattern, docs_all):
            failures.append(f"Docs missing: {label} (pattern: {pattern})")

    for mode in qcd_modes:
        if mode == "auto":
            continue
        if not _contains(rf"\b{re.escape(mode)}\b", docs_all):
            failures.append(f"Docs missing QCD mode: {mode}")

    readme_checks = [
        ("README points to PHYSICS.md", r"PHYSICS\.md"),
        ("README points to CODING.md", r"CODING\.md"),
    ]
    for label, pattern in readme_checks:
        if not _contains(pattern, readme_text):
            failures.append(f"README missing: {label}")

    agents_checks = [
        ("AGENTS has CANONICAL_TERMS", r"CANONICAL_TERMS"),
        ("AGENTS has DO_NOT_EDIT", r"DO_NOT_EDIT"),
        ("AGENTS mentions meta sidecar", r"meta sidecar"),
        ("AGENTS mentions mg5-hnl", r"\bmg5-hnl\b"),
        ("AGENTS mentions hnlcalc scaling validator", r"check_hnlcalc_scaling\.py"),
        ("AGENTS has tau fromTau threshold", rf"\b{fromtau_threshold:.2f}\s*GeV\b"),
        ("AGENTS has FONLL/LHCb reference", r"FONLL/LHCb"),
        ("AGENTS mentions tools/pythia monitor", r"tools/pythia/monitor_production\.sh"),
        ("AGENTS mentions tools/madgraph validator", r"tools/madgraph/validate_xsec\.py"),
        ("AGENTS mentions tools/decay generator", r"tools/decay/generate_hnl_decay_events\.py"),
    ]
    for label, pattern in agents_checks:
        if not _contains(pattern, agents_text):
            failures.append(f"AGENTS missing: {label}")

    do_not_edit_paths = [
        r"production/pythia_production/pythia8315/",
        r"production/madgraph_production/mg5/",
        r"analysis_pbc/decay/external/",
        r"analysis_pbc/HNLCalc/",
    ]
    for path_pattern in do_not_edit_paths:
        if not _contains(re.escape(path_pattern), agents_text):
            failures.append(f"AGENTS missing do-not-edit path: {path_pattern}")

    print("Docs sync check")
    print(f"  MASS_GRID: {mass_count} points ({mass_min:.2f}-{mass_max:.2f} GeV)")
    print(f"  Cross-sections: ccbar={sigma_cc_mb:g} mb, bbbar={sigma_bb_ub:g} microbarn, Bc={sigma_bc_ub:g} microbarn")
    print(f"  Lumi: {lumi_fb:g} fb^-1")
    print(f"  N_limit: {n_limit:.3f}")
    print(f"  fromTau threshold: {fromtau_threshold:.2f} GeV")

    if failures:
        print("\nFAILED")
        for msg in failures:
            print(f"  - {msg}")
        return 1

    print("\nPASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
