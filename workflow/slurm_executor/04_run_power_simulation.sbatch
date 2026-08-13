#!/bin/bash
#SBATCH --job-name=power_04_sim
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/refactor/%x-%A_%a.out
#SBATCH --error=logs/refactor/%x-%A_%a.err
#
# Step 4 -- the real simulation. One array task per (split, effect size).
#
# SUBMIT THIS WITH 04_submit.sh, NOT WITH sbatch DIRECTLY.
#
# The array range depends on N_SPLITS x the number of effect sizes, and #SBATCH directives are
# read by sbatch before the script runs, so they cannot reference a variable sourced at runtime.
# 04_submit.sh computes the range and passes it as --array, which overrides any directive here.
# There is deliberately no --array directive, so a bare `sbatch` fails loudly rather than
# silently running a stale range.
#
# WHY (split x effect size) AND NOT A LOOP OVER EFFECT SIZES
#
# Looping over effect sizes inside each task would share the per-task setup -- the 18 MB
# sim_input.rds and 50 MB sceptre_template.rds reads, plus container and R startup -- across all
# of them. But that setup is ~1 second against a task measured in tens of minutes, while the loop
# multiplies each task's duration by the number of effect sizes and divides parallelism by the
# same factor. With thousands of idle cores on `owners`, trading a little repeated I/O for full
# parallelism is much the better deal: a 10-effect-size sweep goes from ~11 hours to ~70 minutes.
#
# --cpus-per-task=1 is deliberate: sceptre is called with parallel = FALSE and all the
# parallelism comes from the array fan-out. Asking for more cores would idle them.
#
# --time=4h against a measured ~49 min at 1000 splits. The headroom is not superstition: the
# per-target cost varies 4x with pair count (10.3s to 41.4s for 2 reps, measured in the step 3
# smoke test), so a split that happens to collect several large targets runs well over the
# median. Raise it if you lower N_SPLITS -- at 290 splits a task is ~2.8 hours.
#
# Resubmitting: completed tasks are skipped (see the resume check below), so re-running
# 04_submit.sh after a partial failure is safe and picks up only what is missing. To redo
# specific tasks:
#   sbatch --array=17,42,103 workflow/slurm_executor/04_run_power_simulation.sbatch

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: no SLURM_ARRAY_TASK_ID -- this script must be submitted as an array." >&2
  echo "       Use: workflow/slurm_executor/04_submit.sh" >&2
  exit 1
fi

## DECOMPOSE THE TASK ID INTO (split, effect size) ================================================
#
# Split-major within each effect size: tasks 1..N_SPLITS are effect size 1, N_SPLITS+1..2N are
# effect size 2, and so on. So a single effect size can be resubmitted by array range.

if [[ "${SLURM_ARRAY_TASK_ID}" -gt "${TOTAL_SIM_TASKS}" ]]; then
  echo "ERROR: array task ${SLURM_ARRAY_TASK_ID} exceeds ${TOTAL_SIM_TASKS} = ${N_SPLITS} splits x ${N_EFFECT_SIZES} effect sizes." >&2
  echo "       The array range and config.sh have diverged; resubmit with 04_submit.sh." >&2
  exit 1
fi

idx=$(( SLURM_ARRAY_TASK_ID - 1 ))
split_idx=$(( idx % N_SPLITS + 1 ))
es_idx=$(( idx / N_SPLITS ))
ES="${EFFECT_SIZES[${es_idx}]}"

log_header "Step 4: split ${split_idx}/${N_SPLITS}, effect size ${ES}"

load_pixi
require_file "${PREPARED_DIR}/sim_input.rds"        "sim_input.rds from step 01"
require_file "${PREPARED_DIR}/sceptre_template.rds" "sceptre_template.rds from step 01"
require_file "${PREPARED_DIR}/grna_targets.tsv"     "grna_targets.tsv from step 01"
require_dir_nonempty "${SPLITS_DIR}"                "splits directory from step 02"

n_available=$(find "${SPLITS_DIR}" -name 'split_*.tsv' | wc -l)
if [[ "${n_available}" -ne "${N_SPLITS}" ]]; then
  echo "ERROR: ${n_available} split files present but N_SPLITS=${N_SPLITS} in config.sh." >&2
  echo "       Re-run step 02 after changing N_SPLITS." >&2
  exit 1
fi

# By position, not by formatting the index: split_pairs.R pads filenames to the width of
# --n-splits, so 1000 splits give split_0001.tsv but 280 would give split_001.tsv.
SPLIT=$(ls "${SPLITS_DIR}"/split_*.tsv | sed -n "${split_idx}p")
SPLIT_NAME=$(basename "${SPLIT}" .tsv)
n_pairs=$(($(wc -l < "${SPLIT}") - 1))
echo "split file:  ${SPLIT}"
echo "pairs:       ${n_pairs}"

## SIMULATE =======================================================================================

out_dir="${SIM_DIR}/es${ES}"
mkdir -p "${out_dir}"
out="${out_dir}/${SPLIT_NAME}.tsv"

# Resume support. A task that already produced a complete output is skipped, so the array can be
# resubmitted wholesale after a partial failure without redoing finished work.
#
# "Complete" is judged by row count: reps x pairs data rows plus a header. A task killed mid-write
# -- or preempted on `owners` -- leaves a short file, which is caught here rather than silently
# under-counting replicates in compute_power.R.
if [[ -s "${out}" ]]; then
  expected=$(( NUM_REPLICATES * n_pairs + 1 ))
  actual=$(wc -l < "${out}")
  if [[ "${actual}" -eq "${expected}" ]]; then
    echo
    echo ">>> already complete (${actual} lines), skipping"
    report_usage
    exit 0
  fi
  echo
  echo ">>> existing output is incomplete (${actual} of ${expected} lines), redoing"
  rm -f "${out}"
fi

echo
echo ">>> ${NUM_REPLICATES} reps at effect size ${ES} -> ${out}"
pixi_run Rscript bin/run_power_simulation.R \
  --sim-input "${PREPARED_DIR}/sim_input.rds" \
  --sceptre-template "${PREPARED_DIR}/sceptre_template.rds" \
  --pairs "${SPLIT}" \
  --grna-targets "${PREPARED_DIR}/grna_targets.tsv" \
  --effect-size "${ES}" \
  --reps "${NUM_REPLICATES}" \
  --rep-offset 0 \
  --guide-sd "${GUIDE_SD}" \
  --seed "${SEED}" \
  --out "${out}"

require_file "${out}" "simulation output"
echo "    wrote $(($(wc -l < "${out}") - 1)) rows"

report_usage
