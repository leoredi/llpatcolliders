# External inputs needed — BC4 dark scalar (`scalar/`)

The package runs end-to-end today, but two physics inputs are external data
products (not things to invent). They are requested here rather than fabricated.

## 1. Winkler dispersive hadronic width table (improves the model, 0.5–2 GeV)

`scalar/model.py` computes the `pi pi` / `K K` widths in `2 m_pi < m_S < 2 GeV`
with **leading-order ChPT form factors** (Winkler arXiv:1809.01876 eq. 21 → 15;
the "Voloshin" curve of his Fig. 1/3). This is a real ChPT calculation but
misses the `f0(980)` enhancement near 1 GeV and overestimates `K K` toward
2 GeV, and it leaves a ~factor-few seam where it hands over to the perturbative
spectator model at 2 GeV.

**Requested:** Winkler's *dispersive* result — the table behind his Fig. 4:
`m_S` (0.2–2.0 GeV) vs the per-channel widths `Gamma_xx / sin^2 theta`
(`pi pi`, `K K`, total hadronic), or equivalently the dispersive form factors of
his Fig. 2 (`Gamma_pi, Delta_pi, Theta_pi, Gamma_K, Delta_K, Theta_K`). Source:
arXiv:1809.01876 supplementary / the unified FIP calculation arXiv:2311.00507
(which ships the BC4 width tables). Drop in as
`scalar/data/winkler_widths.csv` and switch `model.partial_widths` to read it
below 2 GeV.

## 2. BC4 competitor / existing-bound curves (for the money plot)

`scalar/plot_exclusion.py` overlays competitor curves if present in
`scalar/data/competitors/<name>.csv` (columns `m_S_GeV,sin2theta`). They are
digitized published curves — external data products. Needed files:

| file               | curve                         | type      |
|--------------------|-------------------------------|-----------|
| `CHARM.csv`        | CHARM beam dump               | existing  |
| `LHCb_BKmumu.csv`  | LHCb B→K(*)μμ (displaced)      | existing  |
| `MATHUSLA.csv`     | MATHUSLA projection           | projection|
| `CODEXb.csv`       | CODEX-b projection            | projection|
| `ANUBIS.csv`       | ANUBIS projection             | projection|
| `SHiP.csv`         | SHiP projection               | projection|

**Source:** the 2025 PBC report BC4 figure (arXiv:2505.00947) and the unified
FIP calculation (arXiv:2311.00507); CHARM/LHCb also in Winkler Fig. 8. Digitize
each curve to `(m_S [GeV], sin^2 theta)`. The plot is produced with or without
them (missing ones are skipped with a note).
