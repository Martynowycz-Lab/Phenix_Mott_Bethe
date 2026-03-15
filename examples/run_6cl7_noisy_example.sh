#!/usr/bin/env bash
set -euo pipefail

PACKAGE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO_ROOT="$(cd "$(dirname "${PACKAGE_DIR}")/.." && pwd)"

MODEL="$REPO_ROOT/runs/6cl7_stress_test/thickness_runs/neutral_rmerge_20_thickness_100nm/reference_model.pdb"
DATA="$REPO_ROOT/runs/6cl7_stress_test/thickness_runs/neutral_rmerge_20_thickness_100nm/simulated_obs.mtz"
OUT_DIR="${1:-$REPO_ROOT/runs/6cl7_mb_neutral_example/refine_mb}"

"$PACKAGE_DIR/bin/phenix.refine.mb" \
  --mb-xray-table wk1995 \
  --mb-electron-voltage-kv 200 \
  "$MODEL" \
  "$DATA" \
  data_manager.model.type=electron \
  data_manager.miller_array.labels.name=FSIM,SIGFSIM \
  data_manager.miller_array.labels.name=FREE \
  refinement.main.scattering_table=electron \
  refinement.main.number_of_macro_cycles=3 \
  refinement.output.write_maps=False \
  refinement.output.write_map_coefficients=True \
  overwrite=true \
  output.prefix=phenix_mb_run
