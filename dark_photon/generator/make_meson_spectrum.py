"""Build the committed light-meson (pi0/eta/omega) SoftQCD spectrum for BC1.

Compiles ``dp_meson_softqcd.cc`` against the vendored Pythia 8.315, runs it,
and writes per-species (pT, y) histograms plus the inelastic cross section and
<n_meson>/event to ``data/spectra/meson_softqcd.npz``. The dark-photon
production module samples this committed spectrum (inverse-CDF) and performs
the meson -> A' + X two-body decay analytically, so the heavy Pythia step runs
once and the per-mass generation is fast pure-Python -- the same pattern as the
HNL charged-kaon channel (``make_kaon_spectrum.py``).

Usage:  python generator/make_meson_spectrum.py --n-events 400000
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

SRC = Path(__file__).resolve().parent / "dp_meson_softqcd.cc"
DATA = Path(__file__).resolve().parent.parent / "data" / "spectra" / "meson_softqcd.npz"
PT_MAX, PT_BINS = 20.0, 200
Y_ABS, Y_BINS = 8.0, 160
ECM_GEV = 14000.0
SPECIES = {111: "pi0", 221: "eta", 223: "omega"}


def pythia8_config():
    for parent in Path(__file__).resolve().parents:
        cand = parent / "shared" / "vendored" / "pythia8315" / "bin" / "pythia8-config"
        if cand.exists():
            return str(cand)
    found = shutil.which("pythia8-config")
    if found:
        print(f"WARNING: vendored pythia8315 not found; using {found}", file=sys.stderr)
        return found
    raise SystemExit("pythia8-config not found: need shared/vendored/pythia8315")


def _cfg(cfg, flag):
    return subprocess.check_output([cfg, flag], text=True).split()


def build(cfg, workdir):
    binary = workdir / "dp_meson_softqcd"
    cmd = ["c++", str(SRC), "-o", str(binary),
           *_cfg(cfg, "--cxxflags"), *_cfg(cfg, "--libs")]
    subprocess.run(cmd, check=True)
    return binary


def run(binary, cfg, n_events, seed):
    env = dict(os.environ)
    env["PYTHIA8DATA"] = subprocess.check_output([cfg, "--xmldoc"], text=True).strip()
    proc = subprocess.run([str(binary), str(n_events), str(seed)],
                          capture_output=True, text=True, env=env, check=True)
    meta = {}
    for line in proc.stderr.splitlines():
        parts = line.split()
        if len(parts) == 2 and parts[0].isupper():
            meta[parts[0]] = parts[1]
    valid = {"111", "221", "223"}   # ignore Pythia banner lines on stdout
    pid, pt, y = [], [], []
    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) != 3 or parts[0] not in valid:
            continue
        pid.append(int(parts[0])); pt.append(float(parts[1])); y.append(float(parts[2]))
    return np.asarray(pid), np.asarray(pt), np.asarray(y), meta


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n-events", type=int, default=400_000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("-o", "--output", type=Path, default=DATA)
    args = ap.parse_args(argv)

    cfg = pythia8_config()
    version = subprocess.check_output([cfg, "--version"], text=True).strip()
    pt_edges = np.linspace(0.0, PT_MAX, PT_BINS + 1)
    y_edges = np.linspace(-Y_ABS, Y_ABS, Y_BINS + 1)

    with tempfile.TemporaryDirectory() as tmp:
        binary = build(cfg, Path(tmp))
        pid, pt, y, meta = run(binary, cfg, args.n_events, args.seed)

    n_events = int(meta["N_EVENTS"])
    sigma_inel_mb = float(meta["SIGMA_INEL_MB"])
    out = {
        "pt_edges": pt_edges, "y_edges": y_edges,
        "sigma_inel_mb": np.float64(sigma_inel_mb),
        "n_events": np.int64(n_events),
        "pythia_version": np.asarray(version),
        "ecm_gev": np.float64(ECM_GEV),
    }
    for code, name in SPECIES.items():
        sel = pid == code
        hist, _, _ = np.histogram2d(pt[sel], y[sel], bins=[pt_edges, y_edges])
        n_meson = int(sel.sum())
        out[f"hist_{name}"] = hist
        out[f"n_{name}"] = np.int64(n_meson)
        out[f"n_{name}_per_evt"] = np.float64(n_meson / n_events)
        print(f"  {name:6s}: {n_meson:9d}  <n>/evt = {n_meson / n_events:.4e}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **out)
    print(f"sigma_inel = {sigma_inel_mb:.3f} mb  (Pythia {version}, "
          f"{n_events} events)\nwritten: {args.output}")


if __name__ == "__main__":
    main()
