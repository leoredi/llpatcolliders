"""
hnl/analysis/constants.py

Analysis-level constants for the GRENDEL/GARGOYLE HNL sensitivity scan.

Inherited (with small adjustments) from llpatcolliders_FONLL/analysis/constants.py:
  - DEFAULT_FLAVORS defaults to Umu only (current production focus: muon-mixing
    money plot).
  - FONLL_MASS_MAX raised to 10 GeV so the W/Z-only high-mass region is
    included in the scan (matches the upstream "MATT" hnl_alaship default).
"""

# HL-LHC integrated luminosity
L_INT_FB = 3000.0           # fb-1
L_INT_PB = L_INT_FB * 1e3   # pb-1 = 3e6 pb-1

# Exclusion threshold: N_signal >= N_THRESHOLD (95% CL Poisson, zero bkg)
N_THRESHOLD = 3.0

# U^2 scan range (log10 space)
LOG_U2_MIN = -12.0
LOG_U2_MAX = -1.0
N_U2_POINTS = 200

# CMS IP5 origin for ray-casting
CMS_ORIGIN = (0.0, 0.0, 0.0)

# 2-body acceptance cuts (must match the daughter ray-cast cuts upstream)
SEP_MIN = 0.001    # m -- minimum track separation (1 mm)
SEP_MAX = 1.0      # m -- maximum track separation
P_CUT = 0.600      # GeV/c -- minimum daughter momentum

# Electron mass for 2-body acceptance kinematics
M_ELECTRON = 0.000511  # GeV/c^2

# Full HNL scan range. FONLL meson channels close around 5 GeV but W/Z
# production reaches ~80 GeV, so we keep the full grid.
FONLL_MASS_MAX = 10.0  # GeV

# Flavors known to the analysis (used by reference-curve loader, plot panels).
FLAVORS = ["Ue", "Umu", "Utau"]

# Production focus: muon-only scenario (HNL pattern 010).
DEFAULT_FLAVORS = ["Umu"]
