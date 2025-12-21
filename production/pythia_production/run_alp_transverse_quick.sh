#!/bin/bash
set -euo pipefail

# Quick ALP production suite for transverse detectors:
# - BC10: B_to_Ka (low mass)
# - BC10: h_to_aa (higher mass, Higgs portal style)
# - BC9 : Z_to_gamma_a (photon-coupled style)
#
# Outputs CSVs into output/csv/simulation/ as:
#   ALP_{mass}GeV_{benchmark}_{mode}.csv
#
# Defaults keep runtime modest; tune NEVENTS_* as needed.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHIA_ROOT="$SCRIPT_DIR/pythia8315"

export DYLD_LIBRARY_PATH="$PYTHIA_ROOT/lib:${DYLD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$PYTHIA_ROOT/lib:${LD_LIBRARY_PATH:-}"
export PYTHIA8DATA="$PYTHIA_ROOT/share/Pythia8/xmldoc"

OUTPUT_DIR="$SCRIPT_DIR/../../output/csv/simulation"
LOG_DIR="$SCRIPT_DIR/../../output/logs/simulation"
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

NEVENTS_B="${NEVENTS_B:-20000}"
NEVENTS_H="${NEVENTS_H:-20000}"
NEVENTS_Z="${NEVENTS_Z:-20000}"

build() {
  pushd "$SCRIPT_DIR" >/dev/null
  CXX=${CXX:-g++}
  FLAGS="$($PYTHIA_ROOT/bin/pythia8-config --cxxflags --libs)"
  $CXX $FLAGS -o main_alp_production main_alp_production.cc
  popd >/dev/null
}

run_one() {
  local mass="$1"
  local bench="$2"
  local mode="$3"
  local nevents="$4"
  local mass_label
  mass_label=$(printf "%.2f" "$mass" | tr '.' 'p')

  pushd "$SCRIPT_DIR" >/dev/null
  ./main_alp_production "$mass" "$bench" "$mode" "$nevents" > "${LOG_DIR}/ALP_${mass_label}GeV_${bench}_${mode}.log" 2>&1
  local f="ALP_${mass_label}GeV_${bench}_${mode}.csv"
  if [ -f "$f" ]; then
    mv "$f" "$OUTPUT_DIR/"
  fi
  popd >/dev/null
}

echo "Building..."
build

echo "Running BC10 B_to_Ka (NEVENTS_B=$NEVENTS_B)..."
for m in 0.50 1.00 2.00 3.00 4.00; do
  echo "  m_a=$m"
  run_one "$m" "BC10" "B_to_Ka" "$NEVENTS_B"
done

echo "Running BC10 h_to_aa (NEVENTS_H=$NEVENTS_H)..."
for m in 1.00 5.00 10.00 20.00 40.00 60.00; do
  echo "  m_a=$m"
  run_one "$m" "BC10" "h_to_aa" "$NEVENTS_H"
done

echo "Running BC9 Z_to_gamma_a (NEVENTS_Z=$NEVENTS_Z)..."
for m in 0.10 1.00 5.00 10.00 30.00; do
  echo "  m_a=$m"
  run_one "$m" "BC9" "Z_to_gamma_a" "$NEVENTS_Z"
done

echo "Done. Produced:"
ls -1 "$OUTPUT_DIR"/ALP_*.csv | tail -n 20

