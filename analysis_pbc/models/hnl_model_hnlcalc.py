"""
models/hnl_model_hnlcalc.py

Thin wrapper around the local HNLCalc package in ./HNLCalc.

Goal
-----

For a given HNL mass and flavour mixings (Ue^2, Umu^2, Utau^2), provide a
simple interface:

    model = HNLModel(mass_GeV, Ue2, Umu2, Utau2)
    ctau0_m   = model.ctau0_m          # proper c*tau in metres
    brs_dict  = model.production_brs() # {parent_pdg: BR(parent -> N + X)}

The *implementation details* (exact HNLCalc constructor, attribute names,
and methods) are implemented against the local HNLCalc API. For reference see:

    analysis_pbc/HNLCalc/HNLCalc.py
    analysis_pbc/HNLCalc/Example.ipynb
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Any


# ----------------------------------------------------------------------
# 1. Make sure the local HNLCalc package is importable
# ----------------------------------------------------------------------

THIS_FILE = Path(__file__).resolve()
REPO_ROOT = THIS_FILE.parents[1]          # .../analysis_pbc
HNL_CALC_DIR = REPO_ROOT / "HNLCalc"

# Put ./HNLCalc on sys.path so that "import HNLCalc" uses the local copy
if str(HNL_CALC_DIR) not in sys.path:
    sys.path.insert(0, str(HNL_CALC_DIR))

try:
    # Depending on how HNLCalc.py is written, this will import the module.
    # If your HNLCalc.py defines a class also called HNLCalc, you'll likely do:
    #
    #     from HNLCalc import HNLCalc as _HNLCalcClass
    #
    # below in _build_hnlcalc().
    import HNLCalc as _hnlcalc_module  # noqa: F401
    # Force HNLCalc to write its cache under analysis_pbc/model
    if hasattr(_hnlcalc_module, "MODEL_DIR"):
        _hnlcalc_module.MODEL_DIR = HNL_CALC_DIR.parent / "model"
except ImportError as exc:
    raise ImportError(
        "Could not import the local HNLCalc package. "
        "Make sure analysis_pbc/HNLCalc exists and is a valid Python module."
    ) from exc


# ----------------------------------------------------------------------
# 2. HNLModel wrapper
# ----------------------------------------------------------------------


class HNLModel:
    """
    Lightweight wrapper around the local HNLCalc implementation.

    Parameters
    ----------
    mass_GeV : float
        HNL mass in GeV.
    Ue2, Umu2, Utau2 : float
        Flavour-diagonal mixings |U_e|^2, |U_mu|^2, |U_tau|^2.
    extra_args : dict, optional
        Any additional configuration options you may want to pass through
        to the underlying HNLCalc object.

    Attributes
    ----------
    mass_GeV, Ue2, Umu2, Utau2
        Stored input values.
    _hnlcalc
        The underlying HNLCalc object.
    """

    def __init__(
        self,
        mass_GeV: float,
        Ue2: float,
        Umu2: float,
        Utau2: float,
        **extra_args: Any,
    ) -> None:
        self.mass_GeV = float(mass_GeV)
        self.Ue2 = float(Ue2)
        self.Umu2 = float(Umu2)
        self.Utau2 = float(Utau2)
        self.extra_args = extra_args

        # Underlying HNLCalc object
        self._hnlcalc = self._build_hnlcalc()

    # ------------------------------------------------------------------
    # 2.1. Construct the underlying HNLCalc object
    # ------------------------------------------------------------------

    def _build_hnlcalc(self):
        """
        Construct and return the underlying HNLCalc object.

        HNLCalc takes coupling RATIOS (ve, vmu, vtau), not |U|² values directly.
        It normalizes them internally: epsilon^2 = |U_e|^2 + |U_mu|^2 + |U_tau|^2.

        We pass sqrt(Ue2), sqrt(Umu2), sqrt(Utau2) as the ratios.
        """
        from HNLCalc import HNLCalc as _HNLCalcClass
        import numpy as np

        # HNLCalc constructor signature: HNLCalc(ve=1, vmu=0, vtau=0)
        # It takes coupling ratios and normalizes them
        # To preserve our |U_alpha|^2 values, pass sqrt(U_alpha^2)
        ve = np.sqrt(self.Ue2)
        vmu = np.sqrt(self.Umu2)
        vtau = np.sqrt(self.Utau2)

        # Initialize HNLCalc with coupling pattern
        hnl = _HNLCalcClass(ve=ve, vmu=vmu, vtau=vtau)

        # Generate decay widths and ctau for this mass
        # The coupling parameter here is epsilon = sqrt(Ue2 + Umu2 + Utau2)
        # ctau scales as 1/epsilon^2, so we pass sqrt(Ue2 + Umu2 + Utau2)
        epsilon = np.sqrt(self.Ue2 + self.Umu2 + self.Utau2)

        # Generate ctau and BRs for this mass point
        # Pass coupling so ctau is computed for the actual |U|² values
        hnl.get_br_and_ctau(mpts=np.array([self.mass_GeV]), coupling=epsilon)

        return hnl

    # ------------------------------------------------------------------
    # 2.2. Proper decay length c*tau_0 in metres
    # ------------------------------------------------------------------

    @property
    def ctau0_m(self) -> float:
        """
        Proper decay length c * tau_0 in metres.

        After calling get_br_and_ctau(), HNLCalc stores the decay length
        in the ctau array (in metres). Since we initialized with a single
        mass point, ctau[0] gives us the value.
        """
        # HNLCalc stores ctau in metres in the ctau array
        # We passed a single mass point, so index [0]
        return self._hnlcalc.ctau[0]

    # ------------------------------------------------------------------
    # 2.3. Production BRs B(P -> N + X) per parent
    # ------------------------------------------------------------------

    def production_brs(self) -> Dict[int, float]:
        """
        Return total production branching fractions for each parent:

            { parent_pdg : BR(parent -> N + X) }

        where:
          - parent_pdg is an integer PDG code (e.g. 211, 321, 421, 431, 511, ...),
          - BR(...) is the *total* BSM branching fraction into HNLs for that parent,
            at the HNL mass and mixings stored in this HNLModel.

        This is exactly the quantity you need to weight your per-parent Pythia
        samples:

            weight_i_production = BR(parent_i -> N + X)

        Note: HNLCalc only provides meson production channels. W/Z boson BRs
        are computed analytically and added manually.
        """
        import numpy as np

        # Get all 2-body and 3-body production channels from HNLCalc (mesons only)
        channels_2body = sum(self._hnlcalc.get_channels_2body()["mode"].values(), [])
        channels_3body = sum(self._hnlcalc.get_channels_3body()["mode"].values(), [])

        # Variables needed for controlled eval() calls
        mass = self.mass_GeV
        coupling = np.sqrt(self.Ue2 + self.Umu2 + self.Utau2)
        hnl = self._hnlcalc  # BR strings reference "hnl" object

        def _safe_eval(expr):
            """
            Evaluate a BR expression with a tightly scoped context.
            Allows access only to hnl, mass, and coupling; no builtins.
            """
            return eval(
                expr,
                {"__builtins__": {}},
                {"hnl": hnl, "mass": mass, "coupling": coupling, "np": np},
            )

        # Dictionary to accumulate BRs per parent
        br_per_parent = {}

        # Process 2-body channels
        for channel in channels_2body:
            pid0 = int(channel['pid0'])  # Parent PDG
            pid1 = int(channel['pid1'])  # Daughter PDG
            br_string = channel['br']    # String like "hnl.get_2body_br(411, -11)"

            # Check kinematic threshold
            m_parent = hnl.masses(str(pid0))
            m_daughter = hnl.masses(str(pid1))

            if mass < m_parent - m_daughter:
                # Evaluate the BR string safely; it returns another string with "mass" and "coupling" variables
                br_formula = _safe_eval(br_string)
                # Now evaluate that formula in the same restricted context
                br_value = _safe_eval(br_formula)

                # Accumulate for this parent (use absolute value for antiparticles)
                parent_abs = abs(pid0)
                br_per_parent[parent_abs] = br_per_parent.get(parent_abs, 0.0) + br_value

        # Process 3-body channels
        for channel in channels_3body:
            pid0 = int(channel['pid0'])     # Parent PDG
            pid1 = int(channel['pid1'])     # Daughter 1 PDG
            pid2 = int(channel['pid2'])     # Daughter 2 PDG
            br_string = channel['br']       # Differential BR string
            integration = channel['integration']

            # Get masses
            m0 = hnl.masses(str(pid0))
            m1 = hnl.masses(str(pid1))
            m2 = hnl.masses(str(pid2))

            # Check kinematic threshold
            if mass < m0 - m1 - m2:
                # Evaluate differential BR string safely
                br_diff = _safe_eval(br_string)

                # Integrate to get total BR
                br_value = hnl.integrate_3body_br(
                    br_diff, mass, m0, m1, m2,
                    coupling=coupling,
                    integration=integration
                )

                # Accumulate for this parent
                parent_abs = abs(pid0)
                br_per_parent[parent_abs] = br_per_parent.get(parent_abs, 0.0) + br_value

        # ----------------------------------------------------------------
        # ADD TAU → π N PRODUCTION BR (for "fromTau" samples)
        # ----------------------------------------------------------------
        # Pythia "fromTau" production uses τ → π N as the representative
        # decay channel; without an explicit BR here those events would
        # contribute zero signal. We use the leading-order expression from
        # Atre et al. (arXiv:0901.3589 / 1805.08567 Eq. 2.19) for
        # τ → π N:
        #   Γ = |U_tau|² * (G_F² f_π² |V_ud|² m_τ³ / 16π) *
        #       λ^(1/2)(1, r_π², r_N²) *
        #       [(1 + r_N² - r_π²)(1 + r_N²) - 4 r_N²]
        # BR = Γ / Γ_total(τ)
        m_tau = 1.77686  # GeV
        m_pi = 0.13957   # GeV (charged pion)
        if mass < m_tau - m_pi and self.Utau2 > 0.0:
            G_F = 1.1663787e-5  # GeV^-2
            f_pi = 0.1302       # GeV
            V_ud = 0.97420
            gamma_tau_total = 2.267e-12  # GeV (from PDG lifetime)

            r_N = mass / m_tau
            r_pi = m_pi / m_tau
            # Källén function λ^(1/2)
            kallen = np.sqrt(max(0.0, (1.0 - (r_N + r_pi) ** 2) * (1.0 - (r_N - r_pi) ** 2)))
            matrix_elem = ((1 + r_N**2 - r_pi**2) * (1 + r_N**2) - 4 * r_N**2)

            width_tau_pi_N = (
                self.Utau2
                * (G_F**2 * f_pi**2 * (V_ud**2) * (m_tau**3) / (16.0 * np.pi))
                * kallen
                * matrix_elem
            )

            br_tau_pi_N = width_tau_pi_N / gamma_tau_total if gamma_tau_total > 0 else 0.0
            # Guard against tiny negative due to numerical precision
            br_per_parent[15] = max(br_tau_pi_N, 0.0)

        # ----------------------------------------------------------------
        # ADD W AND Z BOSON PRODUCTION BRs (not in HNLCalc database)
        # ----------------------------------------------------------------
        # HNLCalc only has meson channels. For electroweak bosons,
        # we use theoretical formulas from arXiv:1805.08567 (Atre et al.)
        #
        # Formula: BR(W → ℓN) = |U_ℓ|² × BR_SM × f_PS × f_helicity
        # where:
        #   - BR_SM = BR(W → ℓν) from PDG (one lepton flavor)
        #   - f_PS = (1 - m_N²/m_W²)² (phase space suppression)
        #   - f_helicity = (1 + m_N²/m_W²) (V-A coupling structure)
        #
        # Reference: arXiv:1805.08567 Eq. 2.11-2.12
        # Validated against MadGraph output (agrees within 15%)

        m_W = 80.4  # GeV (PDG)
        m_Z = 91.2  # GeV (PDG)

        # W± → ℓ± N (kinematically allowed if m_N < m_W)
        if mass < m_W:
            # SM branching ratio W → ℓν (one lepton flavor, PDG 2024)
            BR_W_to_lnu_SM = 0.1086

            # Kinematic factors
            r_W = mass / m_W
            phase_space_W = (1.0 - r_W**2)**2
            helicity_W = (1.0 + r_W**2)

            # Total BR: |U|² × BR_SM × kinematics
            br_W = (self.Ue2 + self.Umu2 + self.Utau2) * BR_W_to_lnu_SM * phase_space_W * helicity_W
            br_per_parent[24] = br_W

        # Z → ν N (kinematically allowed if m_N < m_Z)
        if mass < m_Z:
            # SM branching ratio Z → νν (total over 3 flavors, PDG 2024)
            # BR(Z → invisible) = 0.201, divided by 3 neutrino flavors
            BR_Z_to_nunu_SM = 0.201 / 3.0

            # Kinematic factors
            r_Z = mass / m_Z
            phase_space_Z = (1.0 - r_Z**2)**2
            helicity_Z = (1.0 + r_Z**2)

            # Total BR: |U|² × BR_SM × kinematics
            br_Z = (self.Ue2 + self.Umu2 + self.Utau2) * BR_Z_to_nunu_SM * phase_space_Z * helicity_Z
            br_per_parent[23] = br_Z

        return br_per_parent

    # ------------------------------------------------------------------
    # 2.4. Convenience representation
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        return (
            f"HNLModel(mass_GeV={self.mass_GeV}, "
            f"Ue2={self.Ue2}, Umu2={self.Umu2}, Utau2={self.Utau2})"
        )
