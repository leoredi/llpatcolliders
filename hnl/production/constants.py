"""Physical constants for HNL production at LHC 14 TeV."""

from particle import Particle

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
M_KAON = Particle.from_pdgid(321).mass * 1e-3     # K+

# Inclusive pp -> Bc+ + Bc-; already both charges, so no downstream factor 2.
SIGMA_BC_PB = 0.9e6

K_FACTOR_EW = 1.3

# Approximate K+ + K- flux: replace with measured K+- (pT,y) when available.
SIGMA_KAON_PB = 3.0e11  # ~80 mb × ~4 K± per inelastic event (approximate)

KAON_TSALLIS_T = 0.17   # GeV
KAON_TSALLIS_N = 7.0
KAON_PT_MAX = 5.0       # GeV, sampling ceiling (kaons of interest are soft)
KAON_RAPIDITY_SIGMA = 2.5
KAON_E_MAX = 7000.0  # GeV, half of sqrt(s) at LHC 14 TeV

# Meson fractions are absolute weights; baryon fractions are tracked as omitted.
FRAG_B = {
    521: 0.36205648081100655,  # B+/B-
    511: 0.36205648081100655,  # B0/B0bar
    531: 0.08834178131788560,  # Bs
}

FRAG_C = {
    421: 0.382,   # D0/D0bar, includes strong D* feeddown convention
    411: 0.191,   # D+/D-
    431: 0.061,   # Ds+/Ds-
}

OMITTED_FRAG_B = {
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

MESON_MASSES = {
    521: M_BPLUS, 511: M_B0, 531: M_BS, 541: M_BC,
    421: M_D0, 411: M_DPLUS, 431: M_DS,
}

LEPTON_MASSES = {
    'Ue': M_ELECTRON,
    'Umu': M_MUON,
    'Utau': M_TAU,
}

FLAVOR_TO_LEPTON_PDG = {
    'Ue': 11,
    'Umu': 13,
    'Utau': 15,
}

FLAVOR_TO_MG5 = {
    'Ue': 'electron',
    'Umu': 'muon',
    'Utau': 'tau',
}

QUARK_MESON_MAP = {
    'bottom': [(521, FRAG_B[521]), (511, FRAG_B[511]), (531, FRAG_B[531])],
    'charm':  [(421, FRAG_C[421]), (411, FRAG_C[411]), (431, FRAG_C[431])],
}
