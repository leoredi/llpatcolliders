#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

HNL_ENV=/Volumes/sandbox/conda/envs/hnl
export PATH="$HNL_ENV/bin:$PATH"
export CONDA_PREFIX="$HNL_ENV"
export HNL_RUN_TAG="${HNL_RUN_TAG:-full_20260606_all}"
export MPLCONFIGDIR="$PWD/tmp/.mplconfig"
export XDG_CACHE_HOME="$PWD/tmp/.cache"
export PYTHONUNBUFFERED=1

mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

python -u run_all.py \
  --flavor Ue Umu Utau \
  --no-wz \
  --skip-combine \
  --workers 6 \
  --prompt-tau-nb-core 12 \
  2>&1 | tee "tmp/${HNL_RUN_TAG}_prod_nonwz.log"

for flavor in Ue Umu Utau; do
  python -u production/madgraph/run_wz_sharded.py \
    --flavor "$flavor" \
    --jobs 8 \
    --nb-core 1 \
    --skip-existing \
    2>&1 | tee "tmp/${HNL_RUN_TAG}_prod_wz_${flavor}.log"
done

python -m production.combine_channels --flavor Ue Umu Utau \
  2>&1 | tee "tmp/${HNL_RUN_TAG}_combine.log"

python -u run_analysis.py \
  --flavor Ue Umu Utau \
  --workers 3 \
  2>&1 | tee "tmp/${HNL_RUN_TAG}_analysis.log"
