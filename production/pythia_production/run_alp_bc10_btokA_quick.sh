#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHIA_ROOT="$SCRIPT_DIR/pythia8315"

export DYLD_LIBRARY_PATH="$PYTHIA_ROOT/lib:${DYLD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$PYTHIA_ROOT/lib:${LD_LIBRARY_PATH:-}"
export PYTHIA8DATA="$PYTHIA_ROOT/share/Pythia8/xmldoc"

OUTPUT_DIR="$SCRIPT_DIR/../../output/csv/simulation"
LOG_DIR="$SCRIPT_DIR/../../output/logs/simulation"
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

NEVENTS="${NEVENTS:-20000}"
BENCH="${BENCH:-BC10}"

echo "Building main_alp_production..."
pushd "$SCRIPT_DIR" >/dev/null

if [ ! -x ./main_alp_production ]; then
  CXX=${CXX:-g++}
  FLAGS="$($PYTHIA_ROOT/bin/pythia8-config --cxxflags --libs)"
  $CXX $FLAGS -o main_alp_production main_alp_production.cc
fi

echo "Running ALP BC10 B->Ka samples (NEVENTS=$NEVENTS)..."

masses=(0.50 1.00 2.00 3.00 4.00)
for m in "${masses[@]}"; do
  echo "  m_a=$m GeV"
  ./main_alp_production "$m" "$BENCH" "B_to_Ka" "$NEVENTS" > "${LOG_DIR}/ALP_${m}_${BENCH}_B_to_Ka.log" 2>&1
  mass_label=$(printf "%.2f" "$m" | tr '.' 'p')
  f="ALP_${mass_label}GeV_${BENCH}_B_to_Ka.csv"
  if [ -f "$f" ]; then
    mv "$f" "$OUTPUT_DIR/"
  fi
done

popd >/dev/null

echo "Done. Produced:"
ls -1 "$OUTPUT_DIR"/ALP_*_"$BENCH"_B_to_Ka.csv 2>/dev/null || true
