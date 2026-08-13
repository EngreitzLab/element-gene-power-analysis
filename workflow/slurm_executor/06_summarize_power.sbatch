#!/bin/bash
#SBATCH --job-name=power_06_summarize
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Step 6 -- merge the per-effect-size power tables into one wide table.
#
# One row per element-gene pair, one power column per effect size, plus the smallest tested
# effect size at which each pair reaches --power-threshold.
#
# With a single effect size configured this is a thin transformation -- it becomes useful once
# 0.2 is added to EFFECT_SIZES in config.sh, which is also what enables the monotonicity check.
#
# Usage:  sbatch workflow/slurm_executor/06_summarize_power.sbatch

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

log_header "Step 6: summarize power across effect sizes"

load_pixi

power_files=()
for es in "${EFFECT_SIZES[@]}"; do
  f="${POWER_DIR}/power_es${es}.tsv"
  require_file "${f}" "power table for effect size ${es} from step 05"
  power_files+=("${f}")
done

# Comma-separated, in the order the effect sizes are configured.
power_list=$(IFS=,; echo "${power_files[*]}")
out="${OUTDIR}/power_summary.tsv"

echo "inputs:  ${power_list}"
echo

pixi_run Rscript bin/summarize_power.R \
  --power "${power_list}" \
  --power-threshold 0.8 \
  --out "${out}"

require_file "${out}" "power summary"

echo
echo ">>> summary"
echo "  pairs:  $(($(wc -l < "${out}") - 1))"
echo "  columns:"
head -1 "${out}" | tr '\t' '\n' | sed 's/^/    /'
echo "  first rows:"
head -4 "${out}" | cut -c1-200 | sed 's/^/    /'

if [[ "${#EFFECT_SIZES[@]}" -lt 2 ]]; then
  echo
  echo "  Note: only ${#EFFECT_SIZES[@]} effect size configured, so min_detectable_effect_size is"
  echo "  either 0.15 or NA and the monotonicity check is not available. Add 0.2 to"
  echo "  EFFECT_SIZES in config.sh and re-run steps 04-06 to get it."
fi

report_usage
