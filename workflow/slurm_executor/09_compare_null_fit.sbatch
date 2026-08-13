#!/bin/bash
#SBATCH --job-name=power_09_compare_null_fit
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Step 9 -- compare the null_fit subset (step 8) against the same splits from the as_is run (step 4).
#
# Answers: how much did testing simulated counts against real-data coefficients move the reported
# power? The comparison is exactly paired -- the two runs share their simulated count matrices,
# because seeds derive from (seed, target, rep, effect_size) and not from the null model -- so any
# systematic difference is bias, not sampling noise.
#
# Only the splits present in both runs are compared, so this can be run before the step 8 array has
# finished, on whatever has landed.
#
# Usage:  sbatch workflow/slurm_executor/09_compare_null_fit.sbatch

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

ES="${EFFECT_SIZES[0]}"

log_header "Step 9: null_fit vs as_is at effect size ${ES}"

load_pixi
require_file "${PREPARED_DIR}/discovery_threshold.txt" "discovery_threshold.txt from step 01"
require_dir_nonempty "${NULL_FIT_SIM_DIR}/es${ES}"     "null_fit subset from step 08"
require_dir_nonempty "${SIM_DIR}/es${ES}"              "as_is simulation output from step 04"

mkdir -p "${COMPARISON_DIR}"

# Compare only the splits the null_fit subset covers: step 04 has all 1,000 and step 08 has ~20, and
# feeding the full as_is set in would merge on pairs that have no null_fit counterpart. Restricting
# the as_is side to the same basenames keeps the merge honest.
as_is_files=""
n_matched=0
for f in "${NULL_FIT_SIM_DIR}/es${ES}"/split_*.tsv; do
  counterpart="${SIM_DIR}/es${ES}/$(basename "${f}")"
  if [[ -s "${counterpart}" ]]; then
    as_is_files="${as_is_files:+${as_is_files},}${counterpart}"
    n_matched=$(( n_matched + 1 ))
  else
    echo "  NOTE: $(basename "${f}") has no as_is counterpart yet, skipping"
  fi
done

echo "  splits compared: ${n_matched}"
if [[ "${n_matched}" -eq 0 ]]; then
  echo "ERROR: no split has output from both runs." >&2
  exit 1
fi

null_fit_files=""
for f in "${NULL_FIT_SIM_DIR}/es${ES}"/split_*.tsv; do
  if [[ -s "${SIM_DIR}/es${ES}/$(basename "${f}")" ]]; then
    null_fit_files="${null_fit_files:+${null_fit_files},}${f}"
  fi
done

out="${COMPARISON_DIR}/null_fit_vs_as_is_es${ES}.tsv"

pixi_run Rscript workflow/compare_null_fit.R \
  --as-is "${as_is_files}" \
  --null-fit "${null_fit_files}" \
  --threshold-file "${PREPARED_DIR}/discovery_threshold.txt" \
  --out "${out}"

echo
require_file "${out}" "comparison table"
echo ">>> largest per-pair differences"
head -6 "${out}" | sed 's/^/  /'

report_usage
