"""Drell-Yan dark-photon spectra via MadGraph5 + the HAHM UFO model.

Validation-only generator for ``p p > zp`` (q qbar -> A') in the Hidden
Abelian Higgs Model. The production spectra are made by ``dy_pythia.py``;
this independent implementation is retained for high-mass cross-checks.

It uses the Py3-clean ATLAS UFO ``HAHM_darkphoton_LJmod_UFO_GFfix`` (vendored
under ``shared/vendored/madgraph_models/``; gitlab.cern.ch/osalin/
MadGraphModels DarkX). Parses the LHE for the Zp four-vectors and the
integrated cross section, and writes a validation-only (pT, y) spectrum +
``sigma_dy_pb`` at ``eps^2 = 1`` (sigma scales exactly as eps^2 at LO, so
sigma(eps^2=1) = sigma_MG / eps_ref^2). The Zp PDG id is read from the UFO's
particles.py at runtime (this model uses 3000001), never hard-coded.

This mirrors the HNL electroweak production, reusing the MG5 process-runner
plumbing in ``hnl/production/madgraph``. BC1 sets a pure kinetic-mixing dark
photon: HIDDEN block mZDinput = m_A', epsilon = eps_ref, kap = 1e-9 (Higgs
portal off), MHSinput = 200 GeV (dark scalar decoupled).

Usage (hnl conda env, from the worktree root):
    python -m dark_photon.generator.dy_madgraph --masses 0.5 --n-events 20000
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[1]           # dark_photon/
_REPO_ROOT = _ROOT.parent                              # worktree root
for _p in (str(_REPO_ROOT), str(_REPO_ROOT / "hnl")):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def _find_vendored(rel):
    for parent in Path(__file__).resolve().parents:
        cand = parent / "shared" / "vendored" / rel
        if cand.exists():
            return cand
    return None


def _find_workspace_hnl(rel):
    """Find a tool/data path in the workspace's pinned native HNL environment."""
    for parent in Path(__file__).resolve().parents:
        cand = parent / "environments" / "conda" / "hnl" / rel
        if cand.exists():
            return cand
    return None


# HNL's MG5 resolver looks inside its own worktree; point it at shared/vendored.
_mg5 = _find_vendored("MG5_aMC_v3_6_6/bin/mg5_aMC")
if _mg5 is not None:
    os.environ.setdefault("HNL_MG5_EXE", str(_mg5))

# Use the same native LHAPDF installation and NNPDF set as HNL even when the
# caller has not activated that conda environment. Docker is not required.
_lhapdf_config = _find_workspace_hnl("bin/lhapdf-config")
_lhapdf_data = _find_workspace_hnl("share/LHAPDF")
if _lhapdf_config is not None:
    os.environ.setdefault("HNL_LHAPDF_CONFIG", str(_lhapdf_config))
if _lhapdf_data is not None:
    os.environ.setdefault("HNL_LHAPDF_DATA", str(_lhapdf_data))

from production.madgraph.runner import ensure_process_dir, run_events  # noqa: E402
from production.madgraph._mg5_common import patch_me5_configuration  # noqa: E402

# Vendored ATLAS HAHM dark-photon UFO, Py3-clean, GFfix variant (aS at
# SMINPUTS code 3, as MG5 3.x requires): gitlab.cern.ch/osalin/MadGraphModels
# DarkX/HAHM_darkphoton_LJmod_UFO_GFfix.
def _hahm_ufo():
    for parent in Path(__file__).resolve().parents:
        cand = (parent / "shared" / "vendored" / "madgraph_models"
                / "HAHM_darkphoton_LJmod_UFO_GFfix")
        if cand.exists():
            return cand
    raise SystemExit("HAHM UFO not found under shared/vendored/madgraph_models")


HAHM_UFO = _hahm_ufo()
WORK_DIR = _ROOT / "tmp" / "madgraph"
PROCESS_DIR = WORK_DIR / "pp_zp"
SPECTRUM_DIR = _ROOT / "data" / "spectra" / "dy_madgraph_validation"
CARDS_DIR = _ROOT / "generator" / "cards"
PROC_CARD = CARDS_DIR / "proc_pp_zp.dat"
RUN_CARD_TEMPLATE = CARDS_DIR / "run_card_template.dat"
EPS_REF = 1.0e-2                # HAHM epsilon used at generation; sigma ~ eps^2
PT_MAX, PT_BINS = 30.0, 150
Y_ABS, Y_BINS = 8.0, 160


def _zp_pdg(ufo=HAHM_UFO):
    """Read the dark-photon (Zp) PDG id from the UFO's particles.py.

    Read rather than hard-coded so it tracks the vendored model: the ATLAS
    GFfix UFO uses 3000001, unlike Curtin's original v3 (1023).
    """
    txt = (ufo / "particles.py").read_text()
    m = re.search(r"Zp\s*=\s*Particle\(\s*pdg_code\s*=\s*(\d+)", txt)
    if not m:
        raise SystemExit(f"could not find Zp pdg_code in {ufo}/particles.py")
    return int(m.group(1))


ZP_PDG = _zp_pdg()


def _mass_label(m_a):
    return f"{m_a:.3f}".replace(".", "p")


def _build_process():
    """One-time MG5 output of the p p > zp process directory (HAHM).

    The explicit 5-flavor proton (b in p) both includes the b qbar -> A'
    contribution and lets the runner's cache check recognize the built dir, so
    the process is generated once and reused across masses.
    """
    return ensure_process_dir(
        label="pp_zp",
        model_import=f"import model {HAHM_UFO}",
        proc_card=PROC_CARD,
        work_dir=WORK_DIR,
        process_dir=PROCESS_DIR,
        generation_timeout=1200,
    )


# NNPDF4.0 NLO validity floor: Q must stay above the grid's Qmin (~1.65 GeV),
# so DY production is physical only for m_A' at/above this scale. Sub-GeV dark
# photons are produced by the meson channel, not Drell-Yan.
PDF_QMIN = 1.65


def _write_run_card(process_dir, n_events, seed, m_a):
    """Write a complete, deterministic per-mass run card.

    Replacing the generated MG5 default is intentional.  Regex-editing that
    default left some settings (notably ``iseed``) at their automatic values,
    making a failed low-mass run impossible to reproduce exactly.
    """
    cards = process_dir / "Cards"
    run = cards / "run_card.dat"
    # Fixed factorization/renormalization scale at the dark-photon mass,
    # clamped into the PDF-valid range (dynamic m_A'/2 falls below Qmin and
    # the integration diverges for a light resonance).
    scale = max(m_a, PDF_QMIN)
    text = RUN_CARD_TEMPLATE.read_text()
    text = text.replace("N_EVENTS_PLACEHOLDER", str(n_events))
    text = text.replace("SEED_PLACEHOLDER", str(seed))
    text = text.replace("SCALE_PLACEHOLDER", f"{scale:.4f}")
    if "PLACEHOLDER" in text:
        raise RuntimeError("unfilled placeholder in BC1 MadGraph run card")
    run.write_text(text)
    patch_me5_configuration(cards)


def _write_param_card(process_dir, m_a):
    """Configure a pure kinetic-mixing A' of mass m_a in the param card.

    Both the HIDDEN inputs (mZDinput/epsilon/kap/MHSinput) AND the derived
    MASS-block entry for the Zp must be set: MadEvent reads the resonance mass
    from the MASS block, and editing only the HIDDEN input leaves a stale
    default there (which mismatches the factorization scale and breaks the
    integration). The Zp width is left Auto for MG5 to recompute at m_a.
    """
    card = process_dir / "Cards" / "param_card.dat"
    lines = card.read_text().splitlines()
    out, block = [], None
    for line in lines:
        s = line.strip().lower()
        if s.startswith("block"):
            block = s.split()[1] if len(s.split()) > 1 else None
        if block == "hidden" and re.match(r"\s*\d+\s+", line):
            code = int(line.split()[0])
            repl = {1: m_a, 2: 200.0, 3: EPS_REF, 4: 1.0e-9}
            if code in repl:
                out.append(f"    {code} {repl[code]:.6e} # HIDDEN {code}")
                continue
        if block == "mass" and re.match(r"\s*\d+\s+", line):
            code = int(line.split()[0])
            if code == ZP_PDG:                       # dark photon mass
                out.append(f"    {ZP_PDG} {m_a:.6e} # zp : mzdinput")
                continue
            if code == 3000005:                      # dark scalar: decouple
                out.append(f"    3000005 2.000000e+02 # hs : mhsinput")
                continue
        if re.match(rf"\s*DECAY\s+{ZP_PDG}\b", line, flags=re.I):
            # Fixed integration width ~3% of the mass. On-shell production
            # sigma(pp->zp) is width-independent (narrow-width approximation),
            # so this does not change the cross section, but the physical
            # eps^2=1e-2 width (~few MeV) is too narrow for MadEvent's s-channel
            # grid to resolve and the integration NaNs.
            out.append(f"DECAY {ZP_PDG} {max(0.5, 0.05 * m_a):.6e}")
            continue
        out.append(line)
    card.write_text("\n".join(out) + "\n")


def _lhe_sigma(lhe_path):
    """Integrated cross section [pb] from an LHE file (<init> block)."""
    opener = gzip.open if str(lhe_path).endswith(".gz") else open
    with opener(lhe_path, "rt") as fh:
        text = fh.read()
    m = re.search(r"<init>(.*?)</init>", text, flags=re.S)
    if m:
        init_lines = [ln for ln in m.group(1).strip().splitlines() if ln.strip()]
        xsec = 0.0
        for ln in init_lines[1:]:               # skip the beam/PDF header line
            try:
                xsec += float(ln.split()[0])
            except (ValueError, IndexError):
                pass
        if xsec > 0:
            return xsec
    mb = re.search(r"Integrated weight \(pb\)\s*:\s*([0-9.eE+-]+)", text)
    return float(mb.group(1)) if mb else None


def _pythia8_config():
    for parent in Path(__file__).resolve().parents:
        cand = parent / "shared" / "vendored" / "pythia8315" / "bin" / "pythia8-config"
        if cand.exists():
            return str(cand)
    raise SystemExit("vendored pythia8315 not found for LHE showering")


def _build_shower_binary(workdir):
    """Compile dp_dy_shower.cc against the vendored Pythia 8.315."""
    import shutil
    cfg = _pythia8_config()
    binary = workdir / "dp_dy_shower"
    src = Path(__file__).resolve().parent / "dp_dy_shower.cc"
    cxxflags = subprocess.check_output([cfg, "--cxxflags"], text=True).split()
    libs = subprocess.check_output([cfg, "--libs"], text=True).split()
    subprocess.run(["c++", str(src), "-o", str(binary), *cxxflags, *libs], check=True)
    return binary, cfg


def _shower_lhe(lhe_path, m_a, seed):
    """Shower the parton LHE with Pythia ISR; return the A' (pT, y) arrays."""
    import shutil
    # Pythia's LHEF reader wants a plain (uncompressed) file.
    lhe_path = Path(lhe_path)
    if lhe_path.suffix == ".gz":
        plain = lhe_path.with_suffix("")
        with gzip.open(lhe_path, "rb") as fi, open(plain, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        lhe_path = plain
    binary, cfg = _build_shower_binary(WORK_DIR)
    env = dict(os.environ)
    env["PYTHIA8DATA"] = subprocess.check_output([cfg, "--xmldoc"], text=True).strip()
    proc = subprocess.run([str(binary), str(lhe_path), f"{m_a:.6f}", str(seed)],
                          capture_output=True, text=True, env=env, check=True)
    meta = {}
    for line in proc.stderr.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0].isupper():
            meta[parts[0]] = parts[1]
    pt, y = [], []
    zp = str(ZP_PDG)
    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) == 3 and parts[0] == zp:
            pt.append(float(parts[1])); y.append(float(parts[2]))
    return np.asarray(pt), np.asarray(y), meta


def make_dy_spectrum(m_a, n_events=20000, seed=42):
    process_dir = _build_process()
    if process_dir is None:
        raise SystemExit("MG5 process generation failed (see tmp/madgraph logs)")
    _write_run_card(process_dir, n_events, seed, m_a)
    _write_param_card(process_dir, m_a)
    lhe = run_events(process_dir, f"dy_{_mass_label(m_a)}", nb_core=1)
    if lhe is None:
        raise SystemExit(f"MG5 event generation failed for m_A={m_a}")
    sigma_mg = _lhe_sigma(lhe)
    if sigma_mg is None:
        raise SystemExit(f"no cross section in LHE for m_A={m_a}")
    # Shower the parton A' (pT=0) with Pythia ISR to get its transverse recoil.
    pt, y, meta = _shower_lhe(lhe, m_a, seed)
    if len(pt) == 0:
        raise SystemExit(f"shower produced no A' for m_A={m_a} (check dp_dy_shower)")
    sigma_dy_pb = sigma_mg / EPS_REF ** 2          # eps^2 = 1
    pt_edges = np.linspace(0.0, PT_MAX, PT_BINS + 1)
    y_edges = np.linspace(-Y_ABS, Y_ABS, Y_BINS + 1)
    hist, _, _ = np.histogram2d(pt, y, bins=[pt_edges, y_edges])
    SPECTRUM_DIR.mkdir(parents=True, exist_ok=True)
    out = SPECTRUM_DIR / f"dy_{_mass_label(m_a)}.npz"
    np.savez_compressed(
        out, hist=hist, pt_edges=pt_edges, y_edges=y_edges,
        sigma_dy_pb=np.float64(sigma_dy_pb), n_sample=np.int64(len(pt)),
        sigma_mg_pb=np.float64(sigma_mg), eps_ref=np.float64(EPS_REF),
        mass_GeV=np.float64(m_a))
    print(f"  m_A={m_a:.3f}  n_zp={len(pt)}  <pT>={pt.mean():.3f} GeV  "
          f"sigma_MG(eps={EPS_REF})={sigma_mg:.4e} pb  "
          f"sigma(eps^2=1)={sigma_dy_pb:.4e} pb  -> {out.name}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--masses", type=float, nargs="+", required=True)
    ap.add_argument("--n-events", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args(argv)
    for m_a in args.masses:
        make_dy_spectrum(m_a, n_events=args.n_events, seed=args.seed)


if __name__ == "__main__":
    main()
