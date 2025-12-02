# MadGraph HNL Production

**Electroweak HNL production (5-80 GeV) using MadGraph5 in Docker**

---

## Status

✅ **Production complete**: 96 CSV files (32 masses × 3 flavors, ~4.8M events)

---

## Quick Start

```bash
# 1. Build Docker image (once, ~10 min)
cd /path/to/llpatcolliders/production/madgraph_production
docker build -t mg5-hnl .

# 2. Test run (15 GeV muon, 1000 events, ~5 min)
cd /path/to/llpatcolliders
docker run --rm -it -v "$(pwd)":/work mg5-hnl bash
cd /work/production/madgraph_production
python3 scripts/run_hnl_scan.py --test

# 3. Full production
python3 scripts/run_hnl_scan.py --flavour muon  # Single flavor (~8 hrs)
python3 scripts/run_hnl_scan.py                 # All flavors (~30 hrs)
```

---

## Output

**Location**: `csv/HNL_{mass}GeV_{flavor}_EW.csv`

**Format**: Event data with HNL 4-vectors, parent PDG codes, weights

**Summary**: `csv/summary_HNL_EW_production.csv` (cross-sections and metadata)

---

## Usage

```bash
python3 scripts/run_hnl_scan.py [OPTIONS]

Options:
  --test                      Test mode (15 GeV muon, 1k events)
  --flavour {e,mu,tau}        Single flavor
  --masses M1 M2 ...          Custom masses in GeV
  --nevents N                 Events per point (default: 50k)
```

**Examples**:
```bash
# Custom masses
python3 scripts/run_hnl_scan.py --masses 10 20 30

# Low statistics
python3 scripts/run_hnl_scan.py --flavour electron --nevents 10000
```

---

## Why MadGraph?

Pythia fails for m ≥ 5 GeV (off-shell W/Z, kinematic blocking). MadGraph handles W*/Z* propagators correctly via matrix element calculation.

---

## Integration

**Mass coverage**: 5-80 GeV (complements Pythia meson production at 0.2-5 GeV)

**Analysis pipeline**: CSV format compatible with `analysis_pbc/` geometry and limit calculators

---

## Troubleshooting

**Image not found**: `docker build -t mg5-hnl .`

**Hangs on generation**: First run downloads PDFs (~10 min wait), then normal

**Empty /work**: Mount from repo root, not from `madgraph_production/`

**Low event count**: Check `work/hnl_*/Cards/run_card.dat`

---

## Files

- `Dockerfile` - Container definition
- `scripts/run_hnl_scan.py` - Main driver
- `cards/*.dat` - MadGraph input templates
- `csv/` - Output event files (403 MB)
- `work/` - Process directories (26 GB, temporary)

---

## Physics

**Processes**: pp → W±/Z → ℓ N at √s = 14 TeV

**Mass grid**: 32 points from 5-80 GeV

**Mixing**: |U_ℓ|² = 1 at generation (scale in analysis)

**Cross-sections**: 600-25,600 pb (decreasing with mass)

**K-factor**: 1.3 for NLO corrections

---

## References

- MadGraph: [arXiv:1405.0301](https://arxiv.org/abs/1405.0301)
- HNL phenomenology: [arXiv:1805.08567](https://arxiv.org/abs/1805.08567)
- MATHUSLA: [arXiv:1811.00927](https://arxiv.org/abs/1811.00927)

**Main documentation**: `../../CLAUDE.md`
