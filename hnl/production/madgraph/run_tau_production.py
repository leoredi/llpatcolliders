#!/usr/bin/env python3
"""
production/madgraph/run_tau_production.py

Prompt-tau HNL production via the chain

    pp -> W -> tau nu_tau  +  pp -> Z/gamma* -> tau+ tau-
    tau -> N + X           (HNLCalc 2-body and 3-body)

Two-stage driver:

  Stage 1  (run once, flavor- and m_N-independent)
      MG5 + SM model + cards/proc_card_tau_production.dat
      LHE -> 6-column tau pool CSV (w, E, px, py, pz, origin)
      where origin = mother PDG (24 = W+, -24 = W-, 23 = Z)
      Output: vendored/tau_pool.csv  (committed: FONLL-style precedent)

  Stage 2  (per flavor x mass, pure numpy)
      Load the pool, split W vs Z taus, decay each via the shared
      production.decay_engine.tau_decay.sample_hnl_from_tau:
        - W±-origin taus: asymmetry = TAU_2BODY_ASYMMETRY (= -1) — same
          chirality as Ds/B -> tau nu, justified by the (V-A) charged
          current and the CP argument in tau_decay.py.
        - Z-origin taus: asymmetry = 0 (unpolarised; documented
          approximation, the residual Z-mediated polarisation from (V-A)
          interference is small and tau_pool.csv preserves origin so a
          finer model can replace this later without re-running MG5).
      Output: output/llp_4vectors/{flavor}/tau/mN_{mass}.csv

Weight chain per HNL row (m_N < m_tau):

    w_i = (sigma_LO_MG5 / N_pool) * K_FACTOR_EW * BR(tau -> N + X | flavor)

The per-event LHE weight already encodes sigma_LO / N_pool (with split-
event-weight handling for the Z -> tau+ tau- branch, see lhe_to_csv.py).
K_FACTOR_EW is the same flat NLO/LO multiplier the W/Z driver applies;
documented in hnl/README.md.

Usage:
    python run_tau_production.py                          # full
    python run_tau_production.py --skip-mg5               # reuse cached pool
    python run_tau_production.py --test                   # smoke: 1k pool, 1 mass
    python run_tau_production.py --flavor Utau --masses 0.5 1.0
"""

import argparse
import random
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "vendored" / "HNLCalc"))

from config_mass_grid import MASS_GRID, N_EVENTS_DEFAULT, format_mass_for_filename
from production.constants import K_FACTOR_EW, M_TAU
from production.decay_engine.tau_decay import (
    init_hnlcalc, compute_tau_production_br_components, sample_hnl_from_tau,
    TAU_2BODY_ASYMMETRY,
)
from production.madgraph._mg5_common import (
    MG5_EXE, LHAPDF_CONFIG, PYTHON_EXE,
    mg5_subprocess_env, patch_me5_configuration,
    patch_rpath_for_lhapdf, force_compile_subprocesses,
    write_process_block, has_five_flavor_proton,
)
from production.madgraph.lhe_to_csv import LHEParser

CARDS_DIR = Path(__file__).parent / "cards"
WORK_DIR = Path(__file__).parent / "work"
OUTPUT_BASE = PROJECT_ROOT / "output" / "llp_4vectors"

# Vendored tau pool location. FONLL-style precedent: a one-time-generated
# upstream input that downstream loops consume. Committed for the same
# reasons the FONLL .dat tables are committed (see vendored/PROVENANCE.md).
POOL_CSV = PROJECT_ROOT / "vendored" / "tau_pool.csv"


# ---------------------------------------------------------------------------
# Stage 1: MG5 → tau pool
# ---------------------------------------------------------------------------


def _build_process_dir():
    """Build the SM tau-pool MG5 process directory (cached at work/tau_pool/)."""
    work_subdir = WORK_DIR / "tau_pool"

    if (
        (work_subdir / "bin" / "generate_events").exists()
        and has_five_flavor_proton(work_subdir)
    ):
        return work_subdir

    if work_subdir.exists():
        shutil.rmtree(work_subdir, ignore_errors=True)

    proc_card = CARDS_DIR / "proc_card_tau_production.dat"
    if not proc_card.exists():
        raise FileNotFoundError(f"Process card not found: {proc_card}")

    WORK_DIR.mkdir(parents=True, exist_ok=True)

    cmd_file = WORK_DIR / "mg5_gen_tau_pool.txt"
    with open(proc_card) as f:
        proc_lines = f.readlines()

    with open(cmd_file, "w") as f:
        # Tau pool uses the built-in SM model (no HeavyN): no BSM needed at
        # the MG5 stage; HNLCalc handles tau -> N analytically downstream.
        f.write("import model sm\n\n")
        f.write("set automatic_html_opening False\n")
        write_process_block(f, proc_lines)
        f.write(f"\noutput {work_subdir} -nojpeg\n")
        f.write("quit\n")

    log_file = WORK_DIR / "mg5_gen_tau_pool.log"
    with open(log_file, "w") as log:
        result = subprocess.run(
            [str(PYTHON_EXE), str(MG5_EXE), str(cmd_file)],
            stdout=log, stderr=subprocess.STDOUT, timeout=600,
            env=mg5_subprocess_env(),
        )

    if result.returncode != 0 or not (work_subdir / "bin" / "generate_events").exists():
        print(f"    FAILED: process generation (see {log_file})")
        return None

    if not force_compile_subprocesses(work_subdir, WORK_DIR / "mg5_compile_tau_pool.log"):
        print(f"    FAILED: SubProcess pre-compile (see "
              f"{WORK_DIR / 'mg5_compile_tau_pool.log'})")
        return None
    patch_rpath_for_lhapdf(work_subdir)
    cmd_file.unlink(missing_ok=True)
    return work_subdir


def _write_pool_cards(work_subdir, n_events):
    """Write run_card.dat for the tau-pool run.

    The pool is generated at U^2 = 1 with the SM model (no HeavyN), so no
    param_card placeholders apply. We reuse run_card_template.dat (which
    carries the NNPDF40 LHAPDF wiring), with only N_EVENTS substituted.
    """
    cards_dir = work_subdir / "Cards"
    cards_dir.mkdir(exist_ok=True)

    run_content = (CARDS_DIR / "run_card_template.dat").read_text()
    run_content = run_content.replace("N_EVENTS_PLACEHOLDER", str(n_events))
    (cards_dir / "run_card.dat").write_text(run_content)
    patch_me5_configuration(cards_dir)


def _run_pool_events(work_subdir, nb_core):
    """Generate events into Events/run_pool/ ; return path to the LHE."""
    run_name = "run_pool"
    log_file = work_subdir / f"generate_events_{run_name}.log"
    cmd = [
        str(PYTHON_EXE), "bin/generate_events",
        run_name,
        "-f", "--laststep=parton",
    ]
    if nb_core > 1:
        cmd += ["--multicore", f"--nb_core={nb_core}"]
    else:
        cmd += ["--nb_core=1"]

    with open(log_file, "w") as log:
        result = subprocess.run(
            cmd, stdout=log, stderr=subprocess.STDOUT,
            cwd=work_subdir, timeout=7200,
            env=mg5_subprocess_env(),
        )
    if result.returncode != 0:
        print(f"    FAILED: event generation (see {log_file})")
        return None

    run_dir = work_subdir / "Events" / run_name
    for lhe in (run_dir / "unweighted_events.lhe.gz",
                run_dir / "unweighted_events.lhe"):
        if lhe.exists():
            return lhe
    return None


def generate_tau_pool(n_events, nb_core=1, pool_path=None):
    """Stage 1: produce a tau pool CSV via one MG5 run.

    pool_path defaults to the vendored production POOL_CSV. Callers can
    pass a tmp path (e.g. --test mode) to avoid overwriting the vendored
    production pool with smoke-sized statistics.

    Returns True on success, False otherwise.
    """
    if pool_path is None:
        pool_path = POOL_CSV
    print(f"\nStage 1: generating tau pool ({n_events} events, {nb_core} core(s)) → {pool_path}")

    work_subdir = _build_process_dir()
    if work_subdir is None:
        return False

    _write_pool_cards(work_subdir, n_events)

    lhe_path = _run_pool_events(work_subdir, nb_core)
    if lhe_path is None:
        return False

    pool_path.parent.mkdir(parents=True, exist_ok=True)
    n_tau = LHEParser(lhe_path).write_tau_csv(pool_path)
    if n_tau == 0:
        print("    FAILED: no taus extracted from LHE")
        return False
    print(f"    OK: {n_tau} tau rows → {pool_path}")

    # K-factor: scale the per-row LHE weight to sigma_NLO. Matches the W/Z
    # driver convention (see hnl/README.md).
    data = np.loadtxt(pool_path, delimiter=",")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    data[:, 0] *= K_FACTOR_EW
    np.savetxt(pool_path, data, delimiter=",", fmt="%.8e")

    # Per-mass cleanup. The process dir itself is kept for re-runs.
    shutil.rmtree(work_subdir / "Events" / "run_pool", ignore_errors=True)
    return True


# ---------------------------------------------------------------------------
# Stage 2: pool → HNL CSV per (flavor, mass)
# ---------------------------------------------------------------------------


def process_flavor(flavor, pool_w, pool_E, pool_px, pool_py, pool_pz, pool_origin,
                   masses, rng):
    """Decay the pool taus into HNLs for one flavor, all masses."""
    hnl = init_hnlcalc(flavor)
    out_dir = OUTPUT_BASE / flavor / "tau"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Polarization split:
    #   W± mother → fully polarised (helicity ±1 by V-A; CP makes the N
    #              energy spectrum the same for both signs, so a single
    #              scalar asymmetry = TAU_2BODY_ASYMMETRY covers both).
    #   anything else (Z, γ*, direct-DY initial-state partons) → unpolarised.
    # We tag by "is W or not"; this correctly covers γ*/Z Drell-Yan, where
    # the LHE mother of the tau is the initial-state quark rather than a
    # Z propagator, and a strict `origin == 23` check would silently drop
    # the bulk of the Drell-Yan contribution.
    is_W = np.abs(pool_origin) == 24
    is_DY = ~is_W

    for m_N in masses:
        mass_label = format_mass_for_filename(m_N)
        csv_path = out_dir / f"mN_{mass_label}.csv"

        if m_N >= M_TAU:
            csv_path.write_text("")
            continue

        br_2body, br_3body, br_total = compute_tau_production_br_components(hnl, m_N)
        if br_total <= 0:
            csv_path.write_text("")
            continue

        out_rows = []
        for sel, asym, tag in (
            (is_W,  TAU_2BODY_ASYMMETRY, "W"),
            (is_DY, 0.0,                  "DY"),
        ):
            if not sel.any():
                continue
            hnl_4v, _, _ = sample_hnl_from_tau(
                pool_E[sel], pool_px[sel], pool_py[sel], pool_pz[sel],
                m_N, br_2body, br_3body, br_total, rng, asymmetry=asym,
            )
            weights = pool_w[sel] * br_total
            block = np.column_stack([
                weights, hnl_4v[:, 0], hnl_4v[:, 1], hnl_4v[:, 2], hnl_4v[:, 3],
            ])
            out_rows.append(block)

        if out_rows:
            data = np.vstack(out_rows)
            np.savetxt(csv_path, data, delimiter=",", fmt="%.8e")
            print(f"    {csv_path.name}: {len(data)} events, BR_total={br_total:.3e}, "
                  f"w_sum={data[:, 0].sum():.3e}")
        else:
            csv_path.write_text("")


def _load_pool(pool_path=None):
    """Load tau pool CSV → (w, E, px, py, pz, origin) arrays.

    pool_path defaults to the vendored production POOL_CSV; --test mode
    passes the per-test tmp path.
    """
    if pool_path is None:
        pool_path = POOL_CSV
    if not pool_path.exists() or pool_path.stat().st_size == 0:
        raise FileNotFoundError(
            f"Tau pool not found at {pool_path}. Run without --skip-mg5 (Stage 1) "
            f"to generate it."
        )
    data = np.loadtxt(pool_path, delimiter=",")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != 6:
        raise ValueError(
            f"Tau pool {pool_path} has {data.shape[1]} columns; expected 6 "
            f"(w, E, px, py, pz, origin). Regenerate with the current driver."
        )
    return data.T


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(description="Prompt-tau HNL production (W → τν, Z → ττ)")
    ap.add_argument("--flavor", choices=["Ue", "Umu", "Utau"], nargs="+",
                    default=["Ue", "Umu", "Utau"])
    ap.add_argument("--masses", type=float, nargs="+", default=None,
                    help="Custom mass list (default: full grid, m < m_tau)")
    ap.add_argument("--nevents", type=int, default=N_EVENTS_DEFAULT,
                    help="Stage 1 tau pool size (default: %(default)d)")
    ap.add_argument("--nb-core", type=int, default=1,
                    help="MG5 cores for Stage 1")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--skip-mg5", action="store_true",
                    help="Skip Stage 1; reuse cached vendored/tau_pool.csv")
    ap.add_argument("--test", action="store_true",
                    help="Smoke test: Umu only, single mass (1 GeV), 1k pool")
    args = ap.parse_args()

    # --test routes Stage 1 + Stage 2 through a separate tmp pool path so a
    # smoke run cannot clobber the vendored production POOL_CSV.
    if args.test:
        flavors = ["Umu"]
        masses = [1.0]
        n_events = 1000
        pool_path = WORK_DIR / "tau_pool_test.csv"
    else:
        flavors = args.flavor
        masses = args.masses if args.masses else MASS_GRID
        masses = [m for m in masses if m < M_TAU]
        n_events = args.nevents
        pool_path = POOL_CSV

    rng = np.random.default_rng(args.seed)
    random.seed(args.seed)  # HNLCalc's 3-body integrator uses stdlib random

    if not MG5_EXE.exists() and not args.skip_mg5:
        print(f"ERROR: MadGraph not found at {MG5_EXE}")
        print("Set $HNL_MG5_EXE or vendor MG5 under hnl/vendored/ "
              "(see hnl/vendored/PROVENANCE.md).")
        return 1
    if not LHAPDF_CONFIG.exists() and not args.skip_mg5:
        print(f"ERROR: lhapdf-config not found at {LHAPDF_CONFIG}")
        return 1

    print(f"Prompt-tau production")
    print(f"  MG5: {MG5_EXE}")
    print(f"  LHAPDF: {LHAPDF_CONFIG}")
    print(f"  Pool CSV: {pool_path}")
    print(f"  Flavors: {flavors}")
    print(f"  Masses: {len(masses)} points")
    print(f"  Stage 1: {'skipped (using cached pool)' if args.skip_mg5 else f'{n_events} events'}")

    if not args.skip_mg5:
        if not generate_tau_pool(n_events, nb_core=args.nb_core, pool_path=pool_path):
            print("\nStage 1 failed. Aborting.")
            return 1

    pool_w, pool_E, pool_px, pool_py, pool_pz, pool_origin = _load_pool(pool_path)
    n_W  = int((np.abs(pool_origin) == 24).sum())
    n_DY = len(pool_origin) - n_W
    print(f"\nStage 2: decaying pool ({len(pool_origin)} taus: W±={n_W} polarised, "
          f"γ*/Z DY={n_DY} unpolarised)")

    for flavor in flavors:
        print(f"\nFlavor: {flavor}")
        process_flavor(flavor, pool_w, pool_E, pool_px, pool_py, pool_pz, pool_origin,
                       masses, rng)

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
