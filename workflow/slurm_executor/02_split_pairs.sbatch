#!/bin/bash
#SBATCH --job-name=power_02_split
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Step 2 -- split the QC-passing pairs into balanced per-task chunks.
#
# Whole targets stay together: the simulation builds perturbation status and guide-level effect
# sizes once per target, so splitting a target across tasks would duplicate that setup.
#
# n_splits affects scheduling only, not results -- the pipeline is seeded per
# (seed, target, rep, effect_size) and so is invariant to the split layout.
#
# Usage:  sbatch workflow/slurm_executor/02_split_pairs.sbatch

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

log_header "Step 2: split pairs into ${N_SPLITS} chunks"

load_pixi
require_file "${PREPARED_DIR}/pairs.tsv" "pairs.tsv from step 01"

# Splits are regenerated from scratch: a stale split file left over from a different --n-splits
# would be silently picked up by the simulation array's `ls | sed -n` selection.
rm -rf "${SPLITS_DIR}"
mkdir -p "${SPLITS_DIR}"

pixi_run Rscript src/split_pairs.R \
  --pairs "${PREPARED_DIR}/pairs.tsv" \
  --n-splits "${N_SPLITS}" \
  --outdir "${SPLITS_DIR}"

## CHECKS =========================================================================================

echo
echo ">>> splits written"
n_written=$(find "${SPLITS_DIR}" -name 'split_*.tsv' | wc -l)
echo "  files:  ${n_written} (expected ${N_SPLITS})"
if [[ "${n_written}" -ne "${N_SPLITS}" ]]; then
  echo "ERROR: wrote ${n_written} splits but expected ${N_SPLITS}." >&2
  exit 1
fi

# Every pair must appear exactly once across the splits, or the power estimates will be computed
# from the wrong number of replicates. Header rows are excluded from the count.
pairs_in=$(($(wc -l < "${PREPARED_DIR}/pairs.tsv") - 1))
pairs_out=$(cat "${SPLITS_DIR}"/split_*.tsv | grep -vc '^grna_target' || true)
echo "  pairs in:   ${pairs_in}"
echo "  pairs out:  ${pairs_out}"
if [[ "${pairs_in}" -ne "${pairs_out}" ]]; then
  echo "ERROR: pair count changed during splitting (${pairs_in} -> ${pairs_out})." >&2
  exit 1
fi

echo
echo ">>> balance (pairs per split)"
for f in "${SPLITS_DIR}"/split_*.tsv; do echo $(($(wc -l < "$f") - 1)); done \
  | sort -n \
  | awk '{ v[NR] = $1; s += $1 }
         END { printf "  min %d  median %d  max %d  mean %.1f  imbalance %.3fx\n",
                      v[1], v[int((NR + 1) / 2)], v[NR], s / NR, v[NR] / (s / NR) }'

echo
echo "  First split:"
head -3 "$(ls "${SPLITS_DIR}"/split_*.tsv | head -1)" | sed 's/^/    /'

report_usage
