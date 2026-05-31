"""
production/constants.py

Physical constants for HNL production at LHC 14 TeV.

Masses via the `particle` package (PDG 2024).
Cross-sections from FONLL NLO+NLL.
Fragmentation fractions from PDG/HFLAV/ALICE.

References:
  - FONLL: Cacciari, Greco, Nason (NLO+NLL heavy-quark production)
  - Fragmentation fractions: PDG 2024, ALICE D-meson measurements, HFLAV
  - Bc: σ(pp → Bc) from BCVEGPY/FONLL (CMS/LHCb)
"""

from particle import Particle

# ==========================================================================
# Particle masses (GeV) from PDG via `particle` package
# ==========================================================================

M_B0 = Particle.from_pdgid(511).mass * 1e-3      # B0
M_BPLUS = Particle.from_pdgid(521).mass * 1e-3    # B+
M_BS = Particle.from_pdgid(531).mass * 1e-3       # Bs
M_BC = Particle.from_pdgid(541).mass * 1e-3       # Bc+
M_D0 = Particle.from_pdgid(421).mass * 1e-3       # D0
M_DPLUS = Particle.from_pdgid(411).mass * 1e-3    # D+
M_DS = Particle.from_pdgid(431).mass * 1e-3       # Ds+
M_TAU = Particle.from_pdgid(15).mass * 1e-3       # tau
M_ELECTRON = Particle.from_pdgid(11).mass * 1e-3  # electron
M_MUON = Particle.from_pdgid(13).mass * 1e-3      # muon
M_PION = Particle.from_pdgid(211).mass * 1e-3     # pi+

# ==========================================================================
# FONLL inclusive cross-sections at 14 TeV (pb)
# ==========================================================================

# σ(pp → Bc) ~ 0.9 μb at 14 TeV (BCVEGPY/FONLL, CMS/LHCb)
SIGMA_BC_PB = 0.9e6

# ==========================================================================
# Fragmentation fractions
# ==========================================================================

# Central pp fragmentation policy for the HL-LHC off-axis study.
#
# The FONLL tables are generated with fragmentation fraction 1 and provide the
# kinematic meson-shape input. Physical species fractions are applied here in
# the event-weight layer. Species sampling normalizes over the simulated meson
# subset, but the weights below keep the absolute fragmentation fractions.
#
# Charm: ALICE pp, sqrt(s)=13 TeV, |y|<0.5 fragmentation fractions.
# Bottom: LHCb pp, sqrt(s)=13 TeV, 2<eta<5, 4<pT<25 GeV averages with
# fs/(fu+fd)=0.122, f_Lambdab/(fu+fd)=0.259, and fu=fd.
#
# Baryons are tracked as omitted fractions because the current HNLCalc meson
# production layer does not simulate Lambda_c/Xi_c/Lambda_b parents.

# Beauty fragmentation, probabilities per b quark.
FRAG_B = {
    521: 0.36205648081100655,  # B+/B-
    511: 0.36205648081100655,  # B0/B0bar
    531: 0.08834178131788560,  # Bs
}

# Charm fragmentation, probabilities per c quark.
FRAG_C = {
    421: 0.382,   # D0/D0bar, includes strong D* feeddown convention
    411: 0.191,   # D+/D-
    431: 0.061,   # Ds+/Ds-
}

OMITTED_FRAG_B = {
    # Aggregate omitted bottom-baryon remainder inferred from the Lambda_b ratio.
    "b_baryons": 0.18754525706010140,
}

OMITTED_FRAG_C = {
    "Lambda_c+": 0.168,
    "Xi_c0": 0.099,
    "Xi_c+": 0.096,
    "J/psi": 0.0037,
}

FRAGMENTATION_POLICY = {
    "bottom": {
        "source": "LHCb Phys. Rev. D100 (2019) 031102, arXiv:1902.06794",
        "scope": "pp 13 TeV, 2<eta<5, 4<pT<25 GeV; fu=fd; Lambda_b kept as omitted baryon fraction",
        "meson_fraction_sum": sum(FRAG_B.values()),
        "omitted_fraction_sum": sum(OMITTED_FRAG_B.values()),
    },
    "charm": {
        "source": "ALICE charm fragmentation fractions in pp 13 TeV, arXiv:2308.04877",
        "scope": "pp 13 TeV, |y|<0.5; D* counted as feeddown to D0/D+ rather than as an independent weak parent",
        "meson_fraction_sum": sum(FRAG_C.values()),
        "omitted_fraction_sum": sum(OMITTED_FRAG_C.values()),
    },
}

# ==========================================================================
# Meson species lookup tables
# ==========================================================================

# PDG ID → mass (GeV)
MESON_MASSES = {
    521: M_BPLUS, 511: M_B0, 531: M_BS, 541: M_BC,
    421: M_D0, 411: M_DPLUS, 431: M_DS,
}

# Lepton PDG ID → mass (GeV), keyed by flavor label
LEPTON_MASSES = {
    'Ue': M_ELECTRON,
    'Umu': M_MUON,
    'Utau': M_TAU,
}

# Flavor label → lepton PDG ID
FLAVOR_TO_LEPTON_PDG = {
    'Ue': 11,
    'Umu': 13,
    'Utau': 15,
}

# Flavor label → MadGraph flavor name
FLAVOR_TO_MG5 = {
    'Ue': 'electron',
    'Umu': 'muon',
    'Utau': 'tau',
}

# Quark type → list of (pdg_id, frag_fraction) for meson species
QUARK_MESON_MAP = {
    'bottom': [(521, FRAG_B[521]), (511, FRAG_B[511]), (531, FRAG_B[531])],
    'charm':  [(421, FRAG_C[421]), (411, FRAG_C[411]), (431, FRAG_C[431])],
}
