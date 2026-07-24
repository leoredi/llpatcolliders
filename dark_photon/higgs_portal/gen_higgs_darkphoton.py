"""Generate h -> A' A' dark-photon four-vectors for the GRENDEL higgs-portal reach.

Reproduces the collaboration's dark-scalar sample setup (``pythiaStuff/higgsLL.cmnd``)
but for a dark photon of arbitrary mass: inclusive SM Higgs production
(``HiggsSM:all = on``, ggF-dominated) at 13.6 TeV, with the Higgs forced to decay
100% to a pair of long-lived states (pdg 6000113, as in the scalar sample). The
A' is made stable in Pythia so its production four-vector is recorded verbatim;
the downstream GRENDEL machinery (higgs/decayProbPerEvent_2body.py) imposes its
own proper-lifetime scan and 2-body decay, so the A' width/decay here is
irrelevant to the reach.

Writes a CSV with the columns the machinery expects:
    event, id, pt, eta, phi, momentum, mass

Usage:
    python gen_higgs_darkphoton.py --mass 2.0 --n-events 50000 --out tmp/llp_hAA_2GeV.csv
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pythia8

LLP_ID = 6000113  # same pdg code the scalar sample used


def build_pythia(mass_gev: float, ecm_gev: float, seed: int) -> "pythia8.Pythia":
    pythia = pythia8.Pythia()
    cfg = [
        f"Beams:idA = 2212",
        f"Beams:idB = 2212",
        f"Beams:eCM = {ecm_gev}",
        # Inclusive SM Higgs production (ggF + VBF + VH + ttH), as in higgsLL.cmnd.
        "HiggsSM:all = on",
        "25:m0 = 125",
        # Force h -> A' A' at 100%: turn every SM channel off, then add an ON
        # channel (onMode=1, bRatio=1, meMode=100 phase space). Adding it *after*
        # 25:onMode=off keeps it the only open mode, and avoids hard-coding a
        # channel index (robust across Pythia versions).
        "25:onMode = off",
        f"25:addChannel = 1 1.0 100 {LLP_ID} -{LLP_ID}",
        # Keep the Higgs resonance window wide enough that the 2-body decay opens.
        "25:mMin = 100",
        # Define the long-lived state; make it STABLE so we record its production
        # four-vector directly (the reach machinery supplies its own lifetime).
        f"{LLP_ID}:new = Ap Apbar 1 0 0",
        f"{LLP_ID}:m0 = {mass_gev}",
        f"{LLP_ID}:isResonance = off",
        f"{LLP_ID}:mayDecay = off",
        # Run controls
        "Next:numberCount = 0",
        "Print:quiet = on",
        "Init:showChangedSettings = off",
        "Init:showChangedParticleData = off",
        f"Random:setSeed = on",
        f"Random:seed = {seed}",
    ]
    for line in cfg:
        pythia.readString(line)
    if not pythia.init():
        raise SystemExit("Pythia init failed")
    return pythia


def generate(mass_gev, n_events, ecm_gev, seed, out_path):
    pythia = build_pythia(mass_gev, ecm_gev, seed)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_written = 0
    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["event", "id", "pt", "eta", "phi", "momentum", "mass"])
        for iev in range(n_events):
            if not pythia.next():
                continue
            ev = pythia.event
            found = 0
            for i in range(ev.size()):
                p = ev[i]
                if abs(p.id()) != LLP_ID:
                    continue
                # Record only the final (stable) copies -> exactly the two A'.
                if not p.isFinal():
                    continue
                w.writerow([iev, p.id(),
                            f"{p.pT():.6g}", f"{p.eta():.6g}", f"{p.phi():.6g}",
                            f"{p.pAbs():.6g}", f"{p.m():.6g}"])
                found += 1
                n_written += 1
            if found != 2 and iev < 5:
                print(f"  warn: event {iev} had {found} A' (expected 2)")
    xsec_mb = pythia.infoPython().sigmaGen()  # mb
    print(f"Wrote {n_written} A' rows ({n_written//2} h->A'A' events) to {out_path}")
    print(f"Pythia sigmaGen = {xsec_mb:.4g} mb  (inclusive SM Higgs, informational)")
    return out_path


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--mass", type=float, default=2.0, help="A' mass [GeV]")
    ap.add_argument("--n-events", type=int, default=50000)
    ap.add_argument("--ecm", type=float, default=13600.0)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--out", type=str, default="tmp/llp_hAA_2GeV.csv")
    args = ap.parse_args(argv)
    generate(args.mass, args.n_events, args.ecm, args.seed, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
