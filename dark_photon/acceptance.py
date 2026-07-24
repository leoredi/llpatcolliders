"""BC1 dark-photon signal acceptance: A' decays -> shared GRENDEL reco.

Mirrors ``scalar/acceptance.py`` (the BC4 template): the shared GRENDEL mesh
(``higgs/grendel_geometry``) and the bounded 4-hit two-track reconstruction +
selection (``hnl/analysis/decay_reco_acceptance``) are imported, never copied.
Only the decay physics differs -- the A' decay-mode composition comes from
``model.branching_ratios`` (DeLiVeR VMD table below 1.7 GeV) and each mode is
mapped to a two-charged-track proxy exactly as BC4 maps its scalar modes.

Per A' four-vector that hits the fiducial air volume:
  1. sample a decay vertex along the flight path (uniform; the decay-density
     weight is folded in by ``scan_u2``'s lifetime reweighting);
  2. sample a decay MODE by branching ratio, generate the two charged daughters
     as a two-body A' -> x xbar in the rest frame and boost to the lab (neutral
     sub-modes yield no tracks and fail, folding the visible fraction into the
     reconstruction outcome -- no separate BR_vis multiplier);
  3. ray-cast both daughters to the wall, reconstruct, apply ``selection_mask``.

The decay MC is independent of eps^2 (every partial width scales as eps^2, so
the mode composition is coupling-independent), so it is built ONCE per mass and
``scan_u2`` reweights it to every eps^2 by the production (linear) and lifetime
(exponential) factors.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parent
_REPO_ROOT = _ROOT.parent
# Put the repo root (so ``dark_photon`` is importable as a package) and the
# shared hnl/ + higgs/ trees on the path -- but NOT dark_photon/ itself, which
# would shadow the hnl ``production`` package with our ``production.py``.
for _p in (str(_REPO_ROOT), str(_REPO_ROOT / "hnl"), str(_REPO_ROOT / "higgs")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from analysis.decay_reco_acceptance import (          # noqa: E402
    reconstruct_decays, selection_mask, scan_u2,
    signal_contribution_diagnostics,
    HIT_RESOLUTION, P_CUT, CMS_ORIGIN)
from analysis._engine import (                          # noqa: E402
    compute_geometry, _eta_phi_to_directions_batch)
from analysis.exclusion import find_exclusion_band_refined  # noqa: E402
from analysis.format_bridge import load_combined_csv     # noqa: E402
from analysis.constants import L_INT_PB, N_THRESHOLD     # noqa: E402
from production.decay_engine.kinematics import _boost_to_lab  # noqa: E402
from reco_common import SIGMA_T_DEFAULT                  # noqa: E402

from dark_photon import model  # noqa: E402

# Decay mode -> (daughter mass for the 2-body proxy, charged-track fraction).
# The charged fraction encodes neutral sub-modes that give no tracks; multi-body
# hadronic modes are proxied by their two leading charged hadrons (pion mass),
# exactly as BC4 handles its quark/gluon continuum. K0 K0bar gives no prompt
# tracks (K_S displaced, K_L invisible); pi0/eta+gamma modes are fully neutral.
_MODE_DAUGHTER = {
    "ee":       (model.M_ELECTRON, 1.0),
    "mumu":     (model.M_MUON,     1.0),
    "tautau":   (model.M_TAU,      1.0),
    "tau":      (model.M_TAU,      1.0),   # DeLiVeR column name
    "pipi":     (model.M_PIPLUS,   1.0),
    "pi3":      (model.M_PIPLUS,   1.0),   # pi+pi-pi0: leading two charged
    "pi4c":     (model.M_PIPLUS,   1.0),   # 2pi+2pi-
    "pi4n":     (model.M_PIPLUS,   1.0),   # pi+pi-2pi0
    "KKc":      (model.M_KPLUS,    1.0),   # K+K-
    "KKn":      (model.M_KPLUS,    0.0),   # K0 K0bar: no prompt charged tracks
    "PiGamma":  (model.M_PI0,      0.0),   # pi0 gamma: neutral
    "EtaGamma": (model.M_PI0,      0.0),   # eta gamma: mostly neutral
    "ppbar":    (model.M_PROTON,   1.0),
    "had_other":(model.M_PIPLUS,   1.0),
    "hadrons":  (model.M_PIPLUS,   1.0),   # inclusive mode outside the VMD range
}

EVENT_RECO_CHUNK = 20000


def _mode_table(m_a):
    """(modes, probs, daughter_mass, charged_frac) for one A' mass."""
    br = model.branching_ratios(m_a)
    modes = [k for k in br if br[k] > 0.0 and k in _MODE_DAUGHTER]
    probs = np.array([br[k] for k in modes], float)
    probs = probs / probs.sum()
    m_d = np.array([_MODE_DAUGHTER[k][0] for k in modes], float)
    chf = np.array([_MODE_DAUGHTER[k][1] for k in modes], float)
    return modes, probs, m_d, chf


def build_event_mc(p4, direction, entry_d, exit_d, m_a, n_samples, rng,
                   sigma_hit=HIT_RESOLUTION, sigma_t=SIGMA_T_DEFAULT,
                   origin=CMS_ORIGIN, reco_chunk=EVENT_RECO_CHUNK,
                   return_mc=False):
    """Sample decay vertices + A' decays and reconstruct, for a batch of A'
    four-vectors. Returns ``(d, passed)`` each ``(n_events, n_samples)``; both
    are eps^2-independent (``scan_u2`` applies the coupling)."""
    origin = np.asarray(origin, float)
    p4 = np.asarray(p4, float)
    direction = np.asarray(direction, float)
    n_ev = len(entry_d)

    d = rng.uniform(np.asarray(entry_d)[:, None], np.asarray(exit_d)[:, None],
                    size=(n_ev, n_samples))
    M = n_ev * n_samples
    vtx = origin[None, :] + d.reshape(M, 1) * np.repeat(direction, n_samples, axis=0)
    parent = np.repeat(p4, n_samples, axis=0)

    modes, probs, m_d_tab, chf_tab = _mode_table(m_a)
    mi = rng.choice(len(modes), size=M, p=probs)
    m_d = m_d_tab[mi]
    charged = rng.random(M) < chf_tab[mi]
    pstar = np.sqrt(np.clip(0.25 * m_a ** 2 - m_d ** 2, 0.0, None))
    valid = charged & (pstar > 0)

    passed = np.zeros(M, dtype=bool)
    mc_full = None
    if valid.any():
        idx = np.where(valid)[0]
        cos_t = rng.uniform(-1.0, 1.0, len(idx))
        sin_t = np.sqrt(np.clip(1.0 - cos_t ** 2, 0.0, None))
        phi = rng.uniform(0.0, 2 * np.pi, len(idx))
        nhat = np.column_stack([sin_t * np.cos(phi), sin_t * np.sin(phi), cos_t])
        ps = pstar[idx][:, None]
        E_d = np.full(len(idx), m_a / 2.0)
        d1_rest = np.column_stack([E_d, ps * nhat])
        d2_rest = np.column_stack([E_d, -ps * nhat])
        pe, ppx, ppy, ppz = (parent[idx, 0], parent[idx, 1],
                             parent[idx, 2], parent[idx, 3])
        d1 = _boost_to_lab(d1_rest, pe, ppx, ppy, ppz)
        d2 = _boost_to_lab(d2_rest, pe, ppx, ppy, ppz)
        p1 = np.linalg.norm(d1[:, 1:], axis=1)
        p2 = np.linalg.norm(d2[:, 1:], axis=1)
        dir1 = d1[:, 1:] / p1[:, None]
        dir2 = d2[:, 1:] / p2[:, None]
        p_soft = np.minimum(p1, p2)
        b1 = p1 / np.maximum(d1[:, 0], p1)
        b2 = p2 / np.maximum(d2[:, 0], p2)
        ok = (p_soft > P_CUT)
        if ok.any():
            sub = idx[ok]
            ok_idx = np.where(ok)[0]
            for start in range(0, len(sub), reco_chunk):
                stop = min(start + reco_chunk, len(sub))
                take = ok_idx[start:stop]
                mc = reconstruct_decays(
                    vtx[sub[start:stop]], dir1[take], dir2[take], p_soft[take],
                    sigma_hit, sigma_t, rng, beta1=b1[take], beta2=b2[take])
                passed[sub[start:stop]] = selection_mask(mc)

    return d, passed.reshape(n_ev, n_samples)


def _geometry(csv_path, m_a, eta, phi, mesh):
    """Ray-cast with an on-disk cache per mass (mirrors BC4/BC10/HNL)."""
    cache_dir = csv_path.parent.parent / "geometry_cache"
    cache = cache_dir / f"geom_{csv_path.stem}.npz"
    if cache.exists() and cache.stat().st_mtime >= csv_path.stat().st_mtime:
        d = np.load(cache)
        return d["hits"].astype(bool), d["entry_d"], d["exit_d"]
    hits, entry_d, exit_d = compute_geometry(
        eta, phi, mesh, CMS_ORIGIN, batch_label=f"[mA={m_a:.3f}]")
    cache_dir.mkdir(parents=True, exist_ok=True)
    tmp = cache.with_suffix(cache.suffix + ".tmp")
    with open(tmp, "wb") as fh:
        np.savez_compressed(fh, hits=hits, entry_d=entry_d, exit_d=exit_d)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, cache)
    return hits, entry_d, exit_d


def process_mass_point(m_a, mesh, csv_path, eps2_grid, n_samples=100, rng=None):
    """Acceptance + eps^2 scan for one A' mass. Returns its exclusion band dict
    augmented with mass / bookkeeping, or None if production is closed / no hits."""
    csv_path = Path(csv_path)
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return None
    data = load_combined_csv(csv_path, m_a)
    if len(data["weight"]) == 0:
        return None
    if rng is None:
        rng = np.random.default_rng(12345)

    hits, entry_d, exit_d = _geometry(csv_path, m_a, data["eta"], data["phi"], mesh)
    idx = np.where(hits & np.isfinite(entry_d) & np.isfinite(exit_d))[0]
    n_hits = len(idx)
    base = {"mass_GeV": m_a, "n_events": int(len(data["weight"])), "n_hits": n_hits}
    if n_hits == 0:
        return {**base, "eps2_min": np.nan, "eps2_max": np.nan,
                "eps2_min_open": False, "eps2_max_open": False,
                "peak_N": 0.0, "peak_eps2": np.nan, "has_sensitivity": False}

    direction = _eta_phi_to_directions_batch(data["eta"][idx], data["phi"][idx])
    p_mag = data["beta_gamma"][idx] * m_a
    energy = data["gamma"][idx] * m_a
    p4 = np.column_stack([energy, p_mag[:, None] * direction])

    d, passed = build_event_mc(p4, direction, entry_d[idx], exit_d[idx],
                               m_a, n_samples, rng)
    ctau1 = model.ctau_eps2_1(m_a)
    _, N_grid = scan_u2(
        d, passed, exit_d[idx] - entry_d[idx], data["weight"][idx],
        data["beta_gamma"][idx], ctau1, L_INT_PB, eps2_grid)

    def evaluate(eps2):
        _, signal = scan_u2(
            d, passed, exit_d[idx] - entry_d[idx], data["weight"][idx],
            data["beta_gamma"][idx], ctau1, L_INT_PB, np.asarray([eps2]))
        return signal[0]

    result = find_exclusion_band_refined(eps2_grid, N_grid, evaluate, N_THRESHOLD)
    # find_exclusion_band_refined returns u2_min/u2_max keys; relabel to eps2.
    result = dict(result)
    for src, dst in (("u2_min", "eps2_min"), ("u2_max", "eps2_max"),
                     ("u2_min_open", "eps2_min_open"),
                     ("u2_max_open", "eps2_max_open"),
                     ("peak_u2", "peak_eps2")):
        if src in result:
            result[dst] = result.pop(src)
    result.update(base)
    result["ctau_eps2_1_m"] = ctau1
    return result
