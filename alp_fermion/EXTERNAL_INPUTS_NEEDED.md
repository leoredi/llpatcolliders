# External inputs — BC10 fermiophilic ALP

_Updated 2026-07-08._ The two physics inputs originally flagged here as
analytic placeholders have been **replaced by external data products**; what
remains below is the residual systematics they carry and the (deliberately
dropped) overlay curves.

## 1. Absolute B → K a production normalization — SATISFIED

The leading-log placeholder (`C_TOP = 1`) was replaced by the one-loop RG
coefficient of GKOZ (arXiv:2310.03524), evaluated with ALPINIST's
implementation (`tools/compute_cbs_alpinist.py`, BSD-3 port; pinned commit in
`data/alpinist/COMMIT_SHA.txt`):

    g_bs = CBS_EFF * (1/f),  CBS_EFF = 3.518383e-4  at Lambda_UV = 1 TeV.

The finite one-loop + RG terms partially cancel the small log at 1 TeV: the
proper coefficient is ~2.5x *smaller* in amplitude (6.3x in rate) than the
leading-log placeholder. Production now also spans the **full kaon tower**
(K, K0*(700/1430), K*(892/1410/1680), K1(1270/1400), K2*(1430)) for B+ and B0
with the Boiarska et al. (arXiv:1904.10447) form factors — GKOZ's "almost 4x"
uplift over K + K*(892) (we get 3.9-4.1x depending on m_a).

**Residual:** cross-checking CBS_EFF against GKOZ Table 1 (|C_bs| = 1.8e-3 at
f_GKOZ = 1 GeV) agrees to 6-20% in amplitude depending on the m_b scheme used
to unfold their C_bs = c_bs m_b/2f definition. The island's lower
(production-limited) edge inherits <~20% in 1/f. Bs → phi a is not included
(matching ALPINIST); it would add O(10%) production.

## 2. Data-driven hadronic width — SATISFIED

`model.alp_partial_widths` now reads the digitized GKOZ per-channel width
tables (via ALPINIST, `data/alpinist/`, provenance + normalisation
cross-checks in `data/alpinist/PROVENANCE.md`): total hadronic + gamma gamma
over 0.01-3.01 GeV, with the physical eta/eta' mixing poles and the 2 m_c
onset. Above 3.01 GeV the perturbative quark-level sum continues the width,
normalised to the table at the seam (raw mismatch ~1.5x, charm mass-scheme).
All-neutral final states (3pi0, K0 K0bar pi0, pi0 pi0 eta(') with neutral
eta(') decays) are excluded from the visible channels
(`model.visible_fraction`, 0.77-1.0 across the island).

**Residual:** the sharp artificial 2 m_c cliff of the old perturbative
placeholder is gone (the island's upper edge is now data-driven up to
3.01 GeV); above 3.01 GeV the upper edge still rests on the matched
perturbative continuation.

**Convention note (for future overlays / axis labels):** our 1/f axis is the
BNT (arXiv:1708.00443) convention g_aff = c_f m_f / f with c_f = 1. GKOZ
normalise with 1/(2 f_GKOZ), so 1/f_here = 2/f_GKOZ. Any comparison curve
digitized from GKOZ/PBC BC10 figures must be mapped accordingly.

## 3. Competitor + existing-bound curves — DROPPED (by decision)

The GRENDEL island is the deliverable; overlay curves (unified-calc
competitors, CHARM/E137, LHCb/Belle II B→K(*) mumu / B→K+inv) are not part of
it. `plot.py` still accepts an optional long-format overlay CSV
(`curve, m_a_GeV, invf_GeV_inv`) if that decision is ever revisited.
