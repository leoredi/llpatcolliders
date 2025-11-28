#!/bin/bash
# Production run script for HNL benchmark validation
# Tests dual-mode tau production with production-level statistics

set -e  # Exit on error

echo "============================================"
echo "HNL Production Benchmark Suite"
echo "============================================"
echo ""

# Configuration
NEVENTS=200000

# Benchmark 1: 2.6 GeV muon (B-meson dominated, from validation tests)
echo "[1/4] Running 2.6 GeV muon (BC7 benchmark)..."
./main_hnl_production 2.6 muon $NEVENTS
echo ""

# Benchmark 2: 3.0 GeV tau, direct mode (B-meson, MODE A)
echo "[2/4] Running 3.0 GeV tau, direct mode (BC8 MODE A)..."
./main_hnl_production 3.0 tau $NEVENTS direct
echo ""

# Benchmark 3: 1.0 GeV tau, fromTau mode (tau decay, MODE B)
echo "[3/4] Running 1.0 GeV tau, fromTau mode (BC8 MODE B)..."
./main_hnl_production 1.0 tau $NEVENTS fromTau
echo ""

# Benchmark 4: 1.0 GeV tau, direct mode (for comparison)
echo "[4/4] Running 1.0 GeV tau, direct mode (BC8 MODE A)..."
./main_hnl_production 1.0 tau $NEVENTS direct
echo ""

echo "============================================"
echo "Benchmark Complete!"
echo "============================================"
echo ""
echo "Output files:"
ls -lh HNL_*.csv
echo ""
echo "Next: Run analysis pipeline to verify physics results"
