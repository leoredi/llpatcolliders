# MadGraph HNL Pipeline - Installation Guide

This guide covers the setup required before running the MadGraph HNL production pipeline.

## Prerequisites

✅ **Already installed:**
- MadGraph5_aMC@NLO 3.6.6 (`mg5/`)
- Conda environment `mg5env`
- Python 3.11+

❌ **Required (manual installation):**
- HeavyN UFO model

---

## Installing the HeavyN UFO Model

The pipeline requires the `SM_HeavyN_CKM_AllMasses_LO` UFO model, which is not included in the default MadGraph installation.

### Method 1: Download from FeynRules Database

1. Visit the FeynRules Model Database:
   - https://feynrules.irmp.ucl.ac.be/wiki/ModelDatabaseMainPage

2. Search for "HeavyN" or "Sterile Neutrino" models

3. Download the UFO model (e.g., `HeavyN_UFO.tar.gz`)

4. Extract to MadGraph models directory:
   ```bash
   cd production/madgraph/mg5/models/
   tar -xzf /path/to/HeavyN_UFO.tar.gz
   # Rename if needed to SM_HeavyN_CKM_AllMasses_LO
   ```

### Method 2: Clone from GitHub/GitLab

Some HeavyN models are hosted on repositories like:
- https://github.com/
- https://gitlab.com/

Search for "HeavyN UFO" or "Sterile Neutrino UFO" and clone to `mg5/models/`.

### Method 3: FeynRules Generation (Advanced)

If you have access to the FeynRules source for the HeavyN model:

1. Install FeynRules and Mathematica
2. Load the HeavyN model in Mathematica
3. Export to UFO format
4. Copy to `mg5/models/SM_HeavyN_CKM_AllMasses_LO/`

### Verification

After installation, verify the model is available:

```bash
cd production/madgraph
ls mg5/models/ | grep -i heavy

# Should show:
# SM_HeavyN_CKM_AllMasses_LO
```

Test import in MadGraph:

```bash
conda run -n mg5env mg5/bin/mg5_aMC

# In MadGraph prompt:
import model SM_HeavyN_CKM_AllMasses_LO
display particles
# Should show n1, n2, n3 (HNLs)
exit
```

---

## Model Requirements

The HeavyN model must provide:

### Particles
- **n1** (PDG 9900012): Heavy neutral lepton N1
- **n2** (PDG 9900014): Heavy neutral lepton N2 (decoupled)
- **n3** (PDG 9900016): Heavy neutral lepton N3 (decoupled)

### Parameters (in param_card)
- **Block MASS**: Masses for n1, n2, n3
- **Block HNLMIXING**: Mixing parameters
  - ve1, vmu1, vtau1 (N1 mixing with e, μ, τ)
  - ve2, vmu2, vtau2 (N2 mixing - set to 0)
  - ve3, vmu3, vtau3 (N3 mixing - set to 0)

### Couplings
- Must allow processes: pp → W → ℓ N and pp → Z → ν N
- Should handle W/Z off-shell propagators

---

## Alternative: Simplified Model (If HeavyN Unavailable)

If the full HeavyN model is unavailable, you can use a simplified approach:

1. **Use SM + ADD model** (if available):
   - Some ADD/extra dimension models include heavy neutrinos
   - May require process card adjustments

2. **Contact model authors**:
   - Reach out to LLP physics groups (MATHUSLA, CODEX-b)
   - Request access to their MadGraph models

3. **Create custom UFO** (advanced):
   - Minimal SM extension with HNL
   - Requires FeynRules expertise

---

## Known Compatible Models

Based on LLP literature, these models should work:

| Model Name | Source | Notes |
|------------|--------|-------|
| `SM_HeavyN_CKM_AllMasses_LO` | FeynRules DB | Recommended |
| `HeavyN_3gen` | Various | 3 generations |
| `Minimal_HNL` | Custom | Simplified |

**Note**: Exact model names may vary by source. The pipeline expects PDG code 9900012 for N1.

---

## Troubleshooting

### "ImportError: No module named SM_HeavyN_CKM_AllMasses_LO"

**Cause**: Model not installed or wrong name

**Solution**:
```bash
# Check available models
ls mg5/models/

# If model exists but different name, update process cards:
# Edit cards/proc_card_*.dat and change:
# import model YOUR_MODEL_NAME
```

### "Unknown particle n1"

**Cause**: Model doesn't define n1 particle

**Solution**:
- Check model particle definitions
- Verify PDG code matches (should be 9900012)
- May need to edit process cards to use different particle name

### Model installed but MadGraph doesn't recognize it

**Solution**:
```bash
# Clear MadGraph cache
rm -rf mg5/.mg5_cache/

# Restart MadGraph
conda run -n mg5env mg5/bin/mg5_aMC
```

---

## Quick Start After Installation

Once HeavyN model is installed:

```bash
cd production/madgraph

# Test run (15 GeV muon, 1000 events)
conda run -n mg5env python scripts/run_hnl_scan.py --test

# If successful, you should see:
# - MadGraph generating events
# - LHE → CSV conversion
# - Output in csv/muon/HNL_mass_15GeV_muon_EW.csv
```

---

## Additional Resources

- **MadGraph Documentation**: https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki
- **FeynRules**: https://feynrules.irmp.ucl.ac.be/
- **LLP Physics**:
  - MATHUSLA: arXiv:1811.00927
  - CODEX-b: arXiv:1911.00481
  - ANUBIS: arXiv:1909.13022

---

## Contact

For model-specific questions:
- Check with your LLP physics group
- Contact MadGraph mailing list
- See main project documentation in `../../CLAUDE.md`
