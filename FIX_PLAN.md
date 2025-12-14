# Remaining Improvements / TODO (Actionable)

This file tracks **only what is left to do** in `llpatcolliders`. Completed items are intentionally removed to keep this list short and current.

---

## Remaining TODOs (short list)

None (as of this update).

---

## Recent completions

- Geometry preprocessing is now batched/vectorized and handles multi-intersection rays: `analysis_pbc/geometry/per_parent_efficiency.py`.
- Cross-section file now documents assumptions and cites baseline sources: `analysis_pbc/config/production_xsecs.py`.
- HNLCalc wrapper no longer uses `eval()` for BR expressions; an AST whitelist evaluator is used instead: `analysis_pbc/models/hnl_model_hnlcalc.py`.
- HNLCalc on-disk cache is now treated as generated data: `analysis_pbc/model/` is gitignored and `analysis_pbc/model/ctau/ctau.txt` is no longer tracked.

## Quick verification

```bash
cd analysis_pbc
conda run -n llpatcolliders python tests/closure_anubis/test_expected_signal_events_kernel.py
conda run -n llpatcolliders python tests/test_26gev_muon.py
```
