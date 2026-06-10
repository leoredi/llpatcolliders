"""
Full signal selection cutflow for all masses, using decayProbPerEvent_2body's
build_cutflow (so it reflects the active cuts incl. SEP_OUT_GATE). No global
pointing cut (the built-in conditional tight pointing for sep_in < 10 cm stays).

Usage:  python run_cutflows.py
"""
import numpy as np
import decayProbPerEvent_2body as sig
from grendel_geometry import mesh_fiducial

LIFETIME_S = 1e-7
N_PER = 200
MASSES = [('0.5 GeV', 'LLP0p5GeVSmall.csv'),
          ('1 GeV',   'LLP1GeVSmall.csv'),
          ('15 GeV',  'LLPSmall.csv'),
          ('40 GeV',  'LLP40GeVSmall.csv')]


def cutflow(csv):
    geo = sig.cache_geometry(csv, mesh_fiducial, [0, 0, 0])
    mc = sig.sample_separations(geo, LIFETIME_S, n_samples_per_particle=N_PER)
    return sig.build_cutflow(
        mc['sep'], mc['pointing'], mc['weights'], mc['momenta'],
        mc['p_soft'], mc['dca'], mc['vtx_in'], mc['open_angle'],
        mc['sep_outer'], mc['collin'], on_tracker=mc['on_tracker'],
        timing_chi2=mc['timing_chi2'])


def main():
    print(f"SEP_OUT_GATE = {sig.SEP_OUT_GATE*100:.0f} cm, "
          f"COLLIN_FRAC = {sig.COLLIN_FRAC} "
          f"(collin > {sig.COLLIN_FRAC*sig.DETECTOR_THICKNESS*1000:.0f} mm), "
          f"no global pointing cut\n")
    for label, csv in MASSES:
        rows = cutflow(csv)
        print("=" * 72)
        print(f"{label}   ({csv})")
        print(f"  {'cut':<46}{'cum.eff':>9}{'marg.eff':>10}")
        print("  " + "-" * 65)
        for r in rows:
            print(f"  {r['cut']:<46}{r['efficiency']:>9.4f}"
                  f"{r['marginal_efficiency']:>10.4f}")


if __name__ == '__main__':
    main()
