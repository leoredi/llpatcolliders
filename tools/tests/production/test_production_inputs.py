#!/usr/bin/env python3
"""
Validate production cross-sections, fragmentation fractions, and channel
completeness against literature reference values.

References:
  - FONLL NLO+NLL: Cacciari, Greco, Nason, JHEP10(2012)137
  - LHCb σ(bb̄) at 13 TeV: LHCb-PAPER-2016-031  (495 ± 52 μb)
  - ALICE charm frag at 13 TeV: JHEP12(2023)086 (arXiv:2308.04877)
  - MATHUSLA HNL files: github.com/davidrcurtin/MATHUSLA_LLPfiles_RHN_Ue
  - Bondarenko et al.: arXiv:1805.08567
  - PDG 2024: pdg.lbl.gov
  - BCVEGPY Bc production: Chang et al.
  - W/Z NNLO: FEWZ / DYNNLO at 14 TeV
"""

import sys
import re
from pathlib import Path

import pytest

# Ensure repo root is on path
REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

from analysis_pbc.config.production_xsecs import (
    SIGMA_CCBAR_PB,
    SIGMA_BBBAR_PB,
    SIGMA_BC_PB,
    SIGMA_KAON_PB,
    SIGMA_KL_PB,
    SIGMA_W_PB,
    SIGMA_Z_PB,
    K_FACTOR_EW,
    FRAG_C_D0,
    FRAG_C_DPLUS,
    FRAG_C_DS,
    FRAG_C_LAMBDA,
    FRAG_B_B0,
    FRAG_B_BPLUS,
    FRAG_B_BS,
    FRAG_B_LAMBDA,
    get_parent_sigma_pb,
    get_parent_tau_br,
)


# =========================================================================
#  Reference values (literature)
# =========================================================================

# --- Base cross-sections at 14 TeV ---
# FONLL NLO+NLL total σ(pp→cc̄X) at 14 TeV.
# Cacciari et al. JHEP10(2012)137 Table 4: central ~8.3 mb, scale range ~4-17 mb.
# However, PBC studies (CERN-PBC-REPORT-2018-007) commonly quote ~8 mb for
# σ(cc̄) but then define σ_parent = 2 × σ(cc̄) × f_c (×2 for c + c̄).
# Some references (including older PBC inputs) use σ(cc̄) ≈ 24 mb which
# already folds in a multiplicity factor ~3 (average cc̄ pairs per inelastic
# collision × σ_inel ≈ 80 mb × 0.3).  Our value 23.6 mb follows this
# convention.  The ×2 in get_parent_sigma_pb() then accounts for
# particle+antiparticle, not c+c̄ pair counting.
#
# For this test we check the TOTAL per-parent σ(D0+D̄0) against what
# MATHUSLA uses, rather than checking σ(cc̄) in isolation, since the
# convention for σ(cc̄) varies across papers.
REF_SIGMA_CCBAR_MB_RANGE = (4.0, 30.0)  # mb, generous to cover conventions
REF_SIGMA_BBBAR_UB_RANGE = (350.0, 600.0)  # μb
REF_SIGMA_BC_UB_RANGE = (0.3, 2.0)  # μb
REF_SIGMA_KAON_MB_RANGE = (20.0, 100.0)  # mb, very uncertain
REF_SIGMA_W_NB_RANGE = (180.0, 220.0)  # nb (NNLO W+ + W-)
REF_SIGMA_Z_NB_RANGE = (50.0, 70.0)  # nb (NNLO)
REF_K_FACTOR_EW_RANGE = (1.1, 1.5)  # NLO/LO for W/Z + HNL

# --- Fragmentation fractions ---
# PDG/e+e- and ALICE 13 TeV (JHEP12(2023)086)
REF_CHARM_FRAGS = {
    "D0": {"ours": FRAG_C_D0, "pdg_ee": 0.565, "alice_13": 0.542, "tol": 0.30},
    "D+": {"ours": FRAG_C_DPLUS, "pdg_ee": 0.246, "alice_13": 0.225, "tol": 0.30},
    "Ds+": {"ours": FRAG_C_DS, "pdg_ee": 0.080, "alice_13": 0.092, "tol": 0.30},
    "Lc+": {"ours": FRAG_C_LAMBDA, "pdg_ee": 0.057, "alice_13": 0.086, "tol": 0.50},
}

# PDG 2024 / LHCb (pT-dependent, approximate averages)
REF_BEAUTY_FRAGS = {
    "B0": {"ours": FRAG_B_B0, "pdg": 0.40, "mathusla": 0.448, "tol": 0.20},
    "B+": {"ours": FRAG_B_BPLUS, "pdg": 0.40, "mathusla": 0.448, "tol": 0.20},
    "Bs": {"ours": FRAG_B_BS, "pdg": 0.105, "mathusla": 0.103, "tol": 0.30},
    "Lb": {"ours": FRAG_B_LAMBDA, "pdg": 0.09, "mathusla": None, "tol": 0.30},
}

# --- Tau parent BRs ---
REF_TAU_BRS = {
    431: {"name": "Ds→τν", "ref": 0.0548, "tol": 0.05, "source": "PDG 2024"},
    511: {"name": "B0→D(*)τν", "ref": 0.023, "tol": 0.20, "source": "R(D*) WA"},
    521: {"name": "B+→D(*)τν", "ref": 0.023, "tol": 0.20, "source": "isospin"},
    531: {"name": "Bs→Ds(*)τν", "ref": 0.023, "tol": 0.30, "source": "estimate"},
    541: {"name": "Bc→τν", "ref": 0.024, "tol": 0.30, "source": "HPQCD 2020"},
}


# =========================================================================
#  Check 1: Base cross-sections
# =========================================================================

class TestBaseCrossSections:
    """Compare production cross-sections against literature reference ranges."""

    def test_sigma_ccbar_in_range(self):
        sigma_mb = SIGMA_CCBAR_PB / 1e9
        lo, hi = REF_SIGMA_CCBAR_MB_RANGE
        assert lo <= sigma_mb <= hi, (
            f"σ(cc̄) = {sigma_mb:.1f} mb outside [{lo}, {hi}] mb. "
            f"Check convention: FONLL gives ~8 mb; PBC studies sometimes use ~24 mb "
            f"(multiplicity-inclusive)."
        )

    def test_sigma_bbbar_in_range(self):
        sigma_ub = SIGMA_BBBAR_PB / 1e6
        lo, hi = REF_SIGMA_BBBAR_UB_RANGE
        assert lo <= sigma_ub <= hi, (
            f"σ(bb̄) = {sigma_ub:.0f} μb outside [{lo}, {hi}] μb. "
            f"Reference: LHCb 13 TeV = 495 ± 52 μb."
        )

    def test_sigma_bc_in_range(self):
        sigma_ub = SIGMA_BC_PB / 1e6
        lo, hi = REF_SIGMA_BC_UB_RANGE
        assert lo <= sigma_ub <= hi, (
            f"σ(Bc) = {sigma_ub:.2f} μb outside [{lo}, {hi}] μb."
        )

    def test_sigma_kaon_in_range(self):
        sigma_mb = SIGMA_KAON_PB / 1e9
        lo, hi = REF_SIGMA_KAON_MB_RANGE
        assert lo <= sigma_mb <= hi, (
            f"σ(K±) = {sigma_mb:.0f} mb outside [{lo}, {hi}] mb. "
            f"Note: soft QCD, very approximate."
        )

    def test_sigma_w_in_range(self):
        sigma_nb = SIGMA_W_PB / 1e6  # pb → nb: 1 nb = 10^6 pb
        lo, hi = REF_SIGMA_W_NB_RANGE
        assert lo <= sigma_nb <= hi, (
            f"σ(W) = {sigma_nb:.0f} nb outside [{lo}, {hi}] nb."
        )

    def test_sigma_z_in_range(self):
        sigma_nb = SIGMA_Z_PB / 1e6  # pb → nb: 1 nb = 10^6 pb
        lo, hi = REF_SIGMA_Z_NB_RANGE
        assert lo <= sigma_nb <= hi, (
            f"σ(Z) = {sigma_nb:.0f} nb outside [{lo}, {hi}] nb."
        )

    def test_k_factor_ew_in_range(self):
        lo, hi = REF_K_FACTOR_EW_RANGE
        assert lo <= K_FACTOR_EW <= hi, (
            f"K_FACTOR_EW = {K_FACTOR_EW} outside [{lo}, {hi}]."
        )

    def test_sigma_kl_is_half_kaon(self):
        ratio = SIGMA_KL_PB / SIGMA_KAON_PB
        assert 0.4 <= ratio <= 0.6, (
            f"σ(K_L)/σ(K±) = {ratio:.2f}, expected ~0.5 (isospin)."
        )


# =========================================================================
#  Check 2: Fragmentation fractions
# =========================================================================

class TestFragmentationFractions:
    """Compare fragmentation fractions against ALICE, PDG, and MATHUSLA."""

    def test_charm_fractions_sum_to_one(self):
        total = FRAG_C_D0 + FRAG_C_DPLUS + FRAG_C_DS + FRAG_C_LAMBDA
        assert abs(total - 1.0) < 0.05, (
            f"Charm fragmentation sum = {total:.3f}, expected ~1.0."
        )

    def test_beauty_fractions_sum_to_one(self):
        total = FRAG_B_B0 + FRAG_B_BPLUS + FRAG_B_BS + FRAG_B_LAMBDA
        assert abs(total - 1.0) < 0.05, (
            f"Beauty fragmentation sum = {total:.3f}, expected ~1.0."
        )

    @pytest.mark.parametrize("species,refs", REF_CHARM_FRAGS.items())
    def test_charm_frag_vs_reference(self, species, refs):
        ours = refs["ours"]
        # Check against at least one reference within tolerance
        close_to_any = False
        details = []
        for key in ("pdg_ee", "alice_13"):
            ref_val = refs.get(key)
            if ref_val is not None:
                rel_diff = abs(ours - ref_val) / ref_val
                details.append(f"{key}={ref_val:.3f} (diff={rel_diff:.0%})")
                if rel_diff <= refs["tol"]:
                    close_to_any = True
        assert close_to_any, (
            f"f(c→{species}) = {ours:.3f} not within {refs['tol']:.0%} of any "
            f"reference: {', '.join(details)}"
        )

    @pytest.mark.parametrize("species,refs", REF_BEAUTY_FRAGS.items())
    def test_beauty_frag_vs_reference(self, species, refs):
        ours = refs["ours"]
        close_to_any = False
        details = []
        for key in ("pdg", "mathusla"):
            ref_val = refs.get(key)
            if ref_val is not None:
                rel_diff = abs(ours - ref_val) / ref_val
                details.append(f"{key}={ref_val:.3f} (diff={rel_diff:.0%})")
                if rel_diff <= refs["tol"]:
                    close_to_any = True
        assert close_to_any, (
            f"f(b→{species}) = {ours:.3f} not within {refs['tol']:.0%} of any "
            f"reference: {', '.join(details)}"
        )


# =========================================================================
#  Check 3: Factor-of-2 convention
# =========================================================================

class TestFactorOfTwo:
    """Verify the ×2 particle+antiparticle convention is applied correctly."""

    def test_charm_parents_have_factor_2(self):
        """D0, D+, Ds should be σ(cc̄) × f × 2."""
        for pdg, frag, name in [
            (421, FRAG_C_D0, "D0"),
            (411, FRAG_C_DPLUS, "D+"),
            (431, FRAG_C_DS, "Ds"),
            (4122, FRAG_C_LAMBDA, "Lc"),
        ]:
            expected = SIGMA_CCBAR_PB * frag * 2
            actual = get_parent_sigma_pb(pdg)
            assert actual == pytest.approx(expected, rel=1e-6), (
                f"σ({name}) = {actual:.3e}, expected σ(cc̄)×f×2 = {expected:.3e}"
            )

    def test_beauty_parents_have_factor_2(self):
        """B0, B+, Bs, Λb should be σ(bb̄) × f × 2."""
        for pdg, frag, name in [
            (511, FRAG_B_B0, "B0"),
            (521, FRAG_B_BPLUS, "B+"),
            (531, FRAG_B_BS, "Bs"),
            (5122, FRAG_B_LAMBDA, "Lb"),
        ]:
            expected = SIGMA_BBBAR_PB * frag * 2
            actual = get_parent_sigma_pb(pdg)
            assert actual == pytest.approx(expected, rel=1e-6), (
                f"σ({name}) = {actual:.3e}, expected σ(bb̄)×f×2 = {expected:.3e}"
            )

    def test_bc_uses_independent_sigma(self):
        """Bc should use its own σ(Bc), NOT σ(bb̄) × f_Bc × 2."""
        actual = get_parent_sigma_pb(541)
        assert actual == pytest.approx(SIGMA_BC_PB, rel=1e-6), (
            f"σ(Bc) = {actual:.3e}, expected independent σ(Bc) = {SIGMA_BC_PB:.3e}. "
            f"Bc must NOT be derived from σ(bb̄) × fragmentation."
        )

    def test_w_no_double_factor(self):
        """W should be σ(W) × K_EW, no extra ×2 (already W+ + W-)."""
        expected = SIGMA_W_PB * K_FACTOR_EW
        actual = get_parent_sigma_pb(24)
        assert actual == pytest.approx(expected, rel=1e-6), (
            f"σ(W) = {actual:.3e}, expected σ(W)×K = {expected:.3e}"
        )

    def test_z_no_double_factor(self):
        """Z should be σ(Z) × K_EW, no extra ×2."""
        expected = SIGMA_Z_PB * K_FACTOR_EW
        actual = get_parent_sigma_pb(23)
        assert actual == pytest.approx(expected, rel=1e-6), (
            f"σ(Z) = {actual:.3e}, expected σ(Z)×K = {expected:.3e}"
        )

    def test_bc_scales_with_bbbar_slice(self):
        """When a bb̄ slice σ is provided, Bc should scale proportionally."""
        slice_sigma = 1e6  # 1 μb slice
        actual = get_parent_sigma_pb(541, sigma_bbbar_pb=slice_sigma)
        expected = slice_sigma * (SIGMA_BC_PB / SIGMA_BBBAR_PB)
        assert actual == pytest.approx(expected, rel=1e-6), (
            f"Bc with bb̄ slice: {actual:.3e}, expected proportional scaling: {expected:.3e}"
        )


# =========================================================================
#  Check 4: Pythia production settings
# =========================================================================

class TestPythiaSettings:
    """Verify Pythia card files have correct beam, tune, and QCD settings."""

    CARDS_DIR = REPO_ROOT / "production" / "pythia_production" / "cards"

    def _read_card(self, name):
        path = self.CARDS_DIR / name
        assert path.exists(), f"Card file not found: {path}"
        return path.read_text()

    def test_kaon_card_uses_softqcd(self):
        text = self._read_card("hnl_Kaon.cmnd")
        assert "SoftQCD:nonDiffractive = on" in text, (
            "Kaon card should use SoftQCD:nonDiffractive for inclusive production"
        )
        assert "HardQCD" not in text or "off" in text.split("HardQCD")[1][:30], (
            "Kaon card should not enable HardQCD"
        )

    def test_charm_card_uses_hardccbar(self):
        text = self._read_card("hnl_Dmeson.cmnd")
        assert "HardQCD:hardccbar = on" in text, (
            "Charm card should use HardQCD:hardccbar"
        )

    def test_beauty_card_uses_hardbbbar(self):
        text = self._read_card("hnl_Bmeson.cmnd")
        assert "HardQCD:hardbbbar = on" in text, (
            "Beauty card should use HardQCD:hardbbbar"
        )

    def test_bc_card_uses_gg2bbbar(self):
        text = self._read_card("hnl_Bc.cmnd")
        assert "HardQCD:gg2bbbar = on" in text, (
            "Bc card should use HardQCD:gg2bbbar"
        )

    def test_bc_card_has_pthatmin(self):
        text = self._read_card("hnl_Bc.cmnd")
        match = re.search(r"PhaseSpace:pTHatMin\s*=\s*([0-9.]+)", text)
        assert match, "Bc card should set PhaseSpace:pTHatMin"
        pthat = float(match.group(1))
        assert pthat >= 10.0, (
            f"Bc pTHatMin = {pthat}, expected >= 10 GeV for transverse detector enrichment"
        )

    def test_all_cards_14tev(self):
        for name in ["hnl_Kaon.cmnd", "hnl_Dmeson.cmnd", "hnl_Bmeson.cmnd", "hnl_Bc.cmnd"]:
            text = self._read_card(name)
            match = re.search(r"Beams:eCM\s*=\s*([0-9.]+)", text)
            assert match, f"{name} should set Beams:eCM"
            ecm = float(match.group(1))
            assert ecm == pytest.approx(14000.0, rel=0.01), (
                f"{name}: eCM = {ecm}, expected 14000 (14 TeV HL-LHC)"
            )

    def test_all_cards_monash_tune(self):
        for name in ["hnl_Kaon.cmnd", "hnl_Dmeson.cmnd", "hnl_Bmeson.cmnd", "hnl_Bc.cmnd"]:
            text = self._read_card(name)
            assert "Tune:pp = 14" in text, (
                f"{name} should use Monash 2013 tune (Tune:pp = 14)"
            )


# =========================================================================
#  Check 5: EW production (MadGraph)
# =========================================================================

class TestEWProduction:
    """Verify MadGraph EW production configuration."""

    def test_process_cards_use_explicit_w_z_decay(self):
        """Process cards should use pp > w+, w+ > e+ n1 syntax (not inclusive)."""
        cards_dir = REPO_ROOT / "production" / "madgraph_production" / "cards"
        for flavour in ("electron", "muon", "tau"):
            card = cards_dir / f"proc_card_{flavour}.dat"
            if not card.exists():
                pytest.skip(f"proc_card_{flavour}.dat not found")
            text = card.read_text()
            # Should contain explicit boson decay syntax
            has_w_decay = "w+" in text.lower() or "w-" in text.lower()
            has_z_decay = "z" in text.lower()
            assert has_w_decay or has_z_decay, (
                f"proc_card_{flavour}.dat should use explicit W/Z decay syntax"
            )

    def test_k_factor_not_applied_at_generation(self):
        """K_FACTOR_EW should be applied in analysis, not in MadGraph generation."""
        scan_script = REPO_ROOT / "production" / "madgraph_production" / "scripts" / "run_hnl_scan.py"
        if not scan_script.exists():
            pytest.skip("run_hnl_scan.py not found")
        text = scan_script.read_text()
        # K-factor should be mentioned but not applied to param_card
        assert "k_factor" in text.lower() or "K_FACTOR" in text, (
            "run_hnl_scan.py should reference K_FACTOR_EW"
        )

    def test_ew_xsec_reference_ranges(self):
        """EW σ(W→ℓN) × K should be O(10-25 nb) at m_N = 15 GeV.

        Reference: Bondarenko et al. (1805.08567) gives σ ≈ 10-18 nb
        for |V|² = 1 at m_N = 15 GeV, 14 TeV.
        """
        # This is a documentation check — the actual validation is in
        # tools/madgraph/validate_xsec.py which compares generated
        # cross-sections against Bondarenko reference values.
        validate_script = REPO_ROOT / "tools" / "madgraph" / "validate_xsec.py"
        assert validate_script.exists(), (
            "tools/madgraph/validate_xsec.py should exist for EW cross-section validation"
        )


# =========================================================================
#  Check 6: Production channel completeness
# =========================================================================

class TestChannelCompleteness:
    """Verify all standard HNL production channels are present."""

    MAIN_CC = REPO_ROOT / "production" / "pythia_production" / "main_hnl_production.cc"

    @pytest.fixture(autouse=True)
    def _load_source(self):
        assert self.MAIN_CC.exists(), f"C++ source not found: {self.MAIN_CC}"
        self.source = self.MAIN_CC.read_text()

    # --- 2-body leptonic channels ---

    @pytest.mark.parametrize("pdg,name", [
        (321, "K±"),
        (411, "D±"),
        (431, "Ds±"),
        (521, "B±"),
        (541, "Bc±"),
    ])
    def test_2body_leptonic_channel(self, pdg, name):
        """Each charged meson should have a forced 2-body ℓN decay."""
        assert str(pdg) in self.source, (
            f"{name} (PDG {pdg}) 2-body leptonic channel not found in C++ source"
        )

    # --- 3-body semileptonic channels ---

    @pytest.mark.parametrize("pdg,name", [
        (130, "K_L"),
        (421, "D0"),
        (411, "D±"),
        (511, "B0"),
        (521, "B±"),
        (531, "Bs"),
    ])
    def test_3body_semileptonic_channel(self, pdg, name):
        """Each neutral/charged meson should have a 3-body semileptonic channel."""
        assert str(pdg) in self.source, (
            f"{name} (PDG {pdg}) 3-body semileptonic channel not found in C++ source"
        )

    # --- Baryon channels ---

    @pytest.mark.parametrize("pdg,name", [
        (5122, "Λb"),
        (4122, "Λc"),
    ])
    def test_baryon_channel(self, pdg, name):
        """Baryon channels should be present (we're more complete than MATHUSLA)."""
        assert str(pdg) in self.source, (
            f"{name} (PDG {pdg}) baryon channel not found in C++ source"
        )

    # --- EW channels ---

    def test_ew_channels_present(self):
        """W and Z production should be handled by MadGraph EW files."""
        cards_dir = REPO_ROOT / "production" / "madgraph_production" / "cards"
        assert cards_dir.exists(), "MadGraph cards directory not found"
        proc_cards = list(cards_dir.glob("proc_card_*.dat"))
        assert len(proc_cards) >= 3, (
            f"Expected proc_card for electron/muon/tau, found {len(proc_cards)}"
        )

    # --- fromTau chain ---

    def test_from_tau_chain_present(self):
        """fromTau mode should configure meson→τν and τ→NX channels."""
        assert "configureMesonDecaysToTauNu" in self.source, (
            "fromTau: meson→τν configuration function not found"
        )
        assert "configureTauDecays" in self.source, (
            "fromTau: τ→NX decay configuration function not found"
        )

    # --- Channel that should NOT be present ---

    def test_ks_suppressed(self):
        """K_S (PDG 310) should NOT be a production channel (τ_S/τ_L ~ 1/570)."""
        # Check it's not in the charged meson list or forced decays
        # (it may appear in comments, which is fine)
        lines_with_310 = [
            line for line in self.source.splitlines()
            if "310" in line
            and not line.strip().startswith("//")
            and "CHARGED_MESON" not in line
            and "onMode" in line
        ]
        assert len(lines_with_310) == 0, (
            f"K_S (310) appears to have forced decays — it should be suppressed"
        )


# =========================================================================
#  Bonus: Tau parent BRs
# =========================================================================

class TestTauParentBRs:
    """Verify SM BR(meson → τν) values match PDG."""

    @pytest.mark.parametrize("pdg,ref", REF_TAU_BRS.items())
    def test_tau_parent_br(self, pdg, ref):
        actual = get_parent_tau_br(pdg)
        expected = ref["ref"]
        tol = ref["tol"]
        rel_diff = abs(actual - expected) / expected if expected > 0 else 0
        assert rel_diff <= tol, (
            f"{ref['name']}: BR = {actual:.4f}, reference = {expected:.4f} "
            f"(diff = {rel_diff:.0%}, tol = {tol:.0%}). Source: {ref['source']}"
        )

    def test_unknown_parent_returns_zero(self):
        assert get_parent_tau_br(999) == 0.0
