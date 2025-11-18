# Lifetime Independence Validation Results

**Mass Point:** 39GeV

**Validation Date:** 2025-11-18 11:31:37

## Executive Summary

- **KS Test p-values:** min=1.0000, max=1.0000, mean=1.0000
- **Acceptance rate variation:** 0.00%

**VALIDATION STATUS: ✅ PASSED**

Kinematics are statistically independent of lifetime. All KS test p-values > 0.05 and acceptance variation < 10%.

## Dataset Statistics

| Lifetime | Events | Particles | Acceptance Rate |
|----------|--------|-----------|----------------|
| 1ns | 5,000 | 5,000 | 93.12% |
| 30ns | 5,000 | 5,000 | 93.12% |
| 1us | 5,000 | 5,000 | 93.12% |

## Statistical Tests

### Kolmogorov-Smirnov Tests

Tests the null hypothesis that two samples come from the same distribution. p-value > 0.05 indicates distributions are statistically identical.

#### PT

| Comparison | KS Statistic | p-value | Result |
|------------|--------------|---------|--------|
| 1ns_vs_30ns | 0.0000 | 1.0000 | ✅ Same |
| 1ns_vs_1us | 0.0000 | 1.0000 | ✅ Same |
| 30ns_vs_1us | 0.0000 | 1.0000 | ✅ Same |

#### ETA

| Comparison | KS Statistic | p-value | Result |
|------------|--------------|---------|--------|
| 1ns_vs_30ns | 0.0000 | 1.0000 | ✅ Same |
| 1ns_vs_1us | 0.0000 | 1.0000 | ✅ Same |
| 30ns_vs_1us | 0.0000 | 1.0000 | ✅ Same |

#### PHI

| Comparison | KS Statistic | p-value | Result |
|------------|--------------|---------|--------|
| 1ns_vs_30ns | 0.0000 | 1.0000 | ✅ Same |
| 1ns_vs_1us | 0.0000 | 1.0000 | ✅ Same |
| 30ns_vs_1us | 0.0000 | 1.0000 | ✅ Same |

#### MOMENTUM

| Comparison | KS Statistic | p-value | Result |
|------------|--------------|---------|--------|
| 1ns_vs_30ns | 0.0000 | 1.0000 | ✅ Same |
| 1ns_vs_1us | 0.0000 | 1.0000 | ✅ Same |
| 30ns_vs_1us | 0.0000 | 1.0000 | ✅ Same |

## Interpretation

### ✅ Kinematics are Lifetime-Independent

The validation study confirms that HNL kinematics at production do not depend on the particle lifetime (tau0 parameter in PYTHIA). This means:

1. **Option 1 is valid:** We can generate kinematics at one lifetime and re-weight for different lifetimes in post-processing
2. **Option 2 remains valid:** Brute-force approach with separate simulations per lifetime point still works
3. **Option 3 is unnecessary:** We don't need to modify PYTHIA code for custom lifetime handling

### Recommendation

**Proceed with Option 1** (generate once, re-weight in analysis) as it provides:
- Maximum computational efficiency
- Statistical validity (confirmed by this study)
- Clean separation between event generation and physics analysis

## Plots Generated

1. `pt_comparison.png` - Transverse momentum distributions
2. `eta_comparison.png` - Pseudorapidity distributions
3. `phi_comparison.png` - Azimuthal angle distributions
4. `momentum_comparison.png` - Total momentum distributions
5. `acceptance_summary.png` - Detector acceptance rates
6. `correlation_analysis.png` - Angular correlations

## Next Steps

1. Implement Option 1 workflow:
   - Generate large statistics (1M events) at single reference lifetime
   - Modify `decayProbPerEvent.py` to scan over lifetime range
   - Generate exclusion limits with lifetime re-weighting
2. Validate results against Option 2 at a few mass points
3. Document the workflow in repository README
