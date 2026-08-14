#!/bin/bash
#SBATCH --job-name=power_13_test_data
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Build the synthetic sceptre object, then run the whole pipeline on it.
#
# The generator validates its own output against read_sceptre_object()'s requirements, but that only
# proves the object is well formed. What matters is that the five pipeline steps run on it, because
# the point of the object is to let CI exercise the sceptre path without lab data. So this builds it
# and then runs prepare -> split -> null models -> simulate -> power end to end, at a scale that
# finishes in minutes.
#
# Usage:  sbatch workflow/slurm_executor/13_make_test_data.sbatch

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

log_header "Step 13: synthetic test data, and the pipeline on it"

load_pixi

OBJ="${REPO_ROOT}/tests/data/sceptre_object.rds"
WORK="${SCRATCH}/synthetic_pipeline_check"
rm -rf "${WORK}"
mkdir -p "${WORK}"

echo
echo ">>> generating"
pixi_run Rscript src/make_test_data.R --out "${OBJ}"

echo
echo ">>> step 1: prepare"
pixi_run Rscript src/prepare_sim_input.R --sceptre-object "${OBJ}" --outdir "${WORK}/prepared"

# n_splits must not exceed the number of targets; the synthetic object has far fewer than a real
# one, which is exactly the kind of mismatch a CI run should catch rather than inherit.
n_targets=$(( $(cut -f1 "${WORK}/prepared/pairs.tsv" | tail -n +2 | sort -u | wc -l) ))
n_splits=$(( n_targets < 4 ? n_targets : 4 ))
echo "  targets: ${n_targets}, using ${n_splits} splits"

echo
echo ">>> step 2: split"
pixi_run Rscript src/split_pairs.R --pairs "${WORK}/prepared/pairs.tsv" \
  --n-splits "${n_splits}" --outdir "${WORK}/splits"

echo
echo ">>> step 2b/2c: null models (2 replicates)"
pixi_run Rscript src/fit_null_models.R \
  --sim-input "${WORK}/prepared/sim_input.rds" \
  --sceptre-template "${WORK}/prepared/sceptre_template.rds" \
  --grna-targets "${WORK}/prepared/grna_targets.tsv" \
  --reps 2 --rep-offset 0 --seed 1 --out "${WORK}/null_chunk.rds"
pixi_run Rscript src/merge_null_models.R \
  --inputs "${WORK}/null_chunk.rds" --reps 2 --out "${WORK}/null_precomputations.rds"

echo
echo ">>> step 4: simulate (2 replicates, first split)"
first_split=$(ls "${WORK}/splits"/split_*.tsv | head -1)
pixi_run Rscript src/run_power_simulation.R \
  --sim-input "${WORK}/prepared/sim_input.rds" \
  --sceptre-template "${WORK}/prepared/sceptre_template.rds" \
  --pairs "${first_split}" \
  --grna-targets "${WORK}/prepared/grna_targets.tsv" \
  --null-precomputations "${WORK}/null_precomputations.rds" \
  --effect-size 0.5 --reps 2 --rep-offset 0 --guide-sd 0.13 --seed 1 \
  --out "${WORK}/sim.tsv"

echo
echo ">>> step 5: power"
pixi_run Rscript src/compute_power.R \
  --simulations "${WORK}/sim.tsv" \
  --threshold-file "${WORK}/prepared/discovery_threshold.txt" \
  --conf-level 0.95 --out "${WORK}/power.tsv"

echo
echo ">>> step 6: summarize"
pixi_run Rscript src/summarize_power.R \
  --power "${WORK}/power.tsv" --power-threshold 0.8 --out "${WORK}/summary.tsv"

echo
echo ">>> results"
printf '  %-16s %s rows\n' "power.tsv" "$(( $(wc -l < "${WORK}/power.tsv") - 1 ))"
printf '  %-16s %s rows\n' "summary.tsv" "$(( $(wc -l < "${WORK}/summary.tsv") - 1 ))"
echo "  object: $(du -h "${OBJ}" | cut -f1)"

echo
echo "SYNTHETIC PIPELINE CHECK PASSED -- every step runs without lab data."

report_usage
