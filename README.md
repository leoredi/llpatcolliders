# llpatcolliders

HNL sensitivity projections for a transverse LLP detector near CMS at the HL-LHC.

## Documentation

- **[CODING.md](CODING.md)**: **Start here.** Detailed installation, production, and analysis instructions.
- **[PHYSICS.md](PHYSICS.md)**: Physics assumptions, cross-sections, and formulas.
- **[AGENTS.md](AGENTS.md)**: Machine-readable context and file contracts.

## Analysis Approach

This repository implements a complete simulation-to-limit pipeline:

1.  **Production**: Generates HNL events from meson decays (using Pythia 8 for $D, B, B_c$ mesons) and electroweak boson decays (using MadGraph 5 for $W, Z$).
2.  **Geometry & Decay**: Applies detector geometric acceptance and calculates decay probabilities using either detailed precomputed libraries or a calibrated `brvis-kappa` fast surrogate.
3.  **Limit Setting**: Scans the coupling strength $|U|^2$ across a dense mass grid to determine exclusion intervals. The final output is a "money plot" showing the sensitivity "island" in the $(m_N, |U|^2)$ plane where the expected signal exceeds the exclusion threshold.
