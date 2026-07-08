# External inputs needed — BC10 fermiophilic ALP

The model layer (`model.py`) uses standard analytic formulae from the
fermiophilic-ALP literature with every constant made explicit. Two pieces are
genuine **external data products** that should be supplied to replace the
documented analytic placeholders before the island is quoted as final. They are
NOT fabricated here; the placeholders are clearly bounded and their effect on
the reach is stated.

## 1. Absolute B → K a production normalization
`model.h_sb` / `model.width_B_to_K_a` implement the top-W penguin at
**leading log** with an O(1) coefficient `C_TOP = 1` and a fixed UV scale
`LAMBDA_UV = 1 TeV`, plus a single-pole B→K scalar form factor
`f_0(q²)` (`F0_B_K_AT_0 = 0.33`, pole² = 37.5 GeV²).

Needed: the BC10 value of BR(B→K a)(1/f) as tabulated in the unified FIP
calculation **arXiv:2311.00507** (the source the 2025 PBC report
arXiv:2505.00947 uses), or equivalently the matched Wilson coefficient and the
lattice f_0(q²) (e.g. ALPINIST, arXiv:2105.10806; HPQCD B→K form factors).

Impact: BR(B→K a) ∝ (C_TOP · f_0)², so the **production yield ∝ (C_TOP f_0)²**
and the island's coupling edges move as `1/f ∝ (C_TOP f_0)^{-1}`. A factor-2
error in C_TOP·f_0 shifts the whole island by ~2× in 1/f. The reported reach
should be read with this scaling until the normalization is matched.

Extension (uplift, not included): the vector channels B→K\* a and B_s→φ a add an
O(1) factor to the production rate; only the pseudoscalar B→K a channels (B⁺→K⁺,
B⁰→K⁰) are in the baseline.

## 2. Data-driven hadronic width / spectral function (0.3–~2 GeV)
`model.alp_partial_widths["hadronic"]` is the **perturbative quark-level** sum
(N_c=3, current masses, 1+α_s/π), gated at m_a > 2 m_π. In the resonance region
the true a→hadrons rate is data-driven (ALP–π⁰/η/η′ mixing, R-ratio), which is
what arXiv:2311.00507 provides.

Needed: the data-driven a→hadrons partial width (or the hadronic R-ratio /
exclusive ChPT modes π⁺π⁻π⁰, ηππ, …) vs m_a.

Impact: this sets cτ and the visible-track mix above 2 m_π. The perturbative
model has a **sharp, artificial cc̄ onset at m_a ≈ 2 m_c ≈ 2.55 GeV** that
over-shortens cτ and pinches the island closed near 2.5–3 GeV. The
**clean, robust part of the result is the low-mass μμ-dominated window**
(2 m_μ ≲ m_a ≲ 1 GeV, and especially the pure-dimuon 2 m_μ–2 m_π band); above
~2.5 GeV the upper island edge is placeholder-limited.

## 3. Competitor + existing-bound curves for the overlay plot
`plot.py` draws the GRENDEL island and overlays curves from an optional
long-format CSV (`curve, m_a_GeV, invf_GeV_inv`). The unified-calc competitor
(arXiv:2311.00507) and the existing-bound contours (CHARM/E137 beam dumps,
LHCb/Belle II B→K(\*)μμ & B→K+inv, etc.) are external digitizations to drop in
there. Without them the plot shows only the GRENDEL island.
