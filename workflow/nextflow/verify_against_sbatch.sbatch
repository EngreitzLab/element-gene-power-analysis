#!/bin/bash
#SBATCH --job-name=power_nf_verify
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/nextflow/%x-%A.out
#SBATCH --error=logs/nextflow/%x-%A.err
#
# Check the Nextflow runner against the sbatch runner, output by output.
#
# The two runners exist side by side on purpose: the sbatch scripts produced every measurement in
# docs/status.md, so they are the reference until Nextflow reproduces them. This is the check that
# decides when that has happened. Run it after a Nextflow run and it says, per artefact, whether the
# two agree.
#
# Byte-identity is the bar for everything except the null-model bundle, which the two runners may
# cover to different depths (a test run fits 2 replicates, a production run 100). For that one the
# comparison is over the shared replicates, and it should still be exact -- null models are seeded
# from (seed, rep) alone, so replicate 1 is the same object either way.
#
# Usage:
#   sbatch workflow/nextflow/verify_against_sbatch.sbatch [nextflow_outdir]
#
# Defaults to the test tree; pass results/nextflow for a production run.

set -uo pipefail

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
REPO_ROOT="$(cd "$(dirname "${_cmd:-${BASH_SOURCE[0]}}")/../.." && pwd)"
cd "${REPO_ROOT}"
mkdir -p logs/nextflow

source workflow/slurm_executor/config.sh

NF_OUT="${1:-${REPO_ROOT}/results/nextflow_test}"
NF_DIR="${NF_OUT}/${SAMPLE}"

log_header "Verify: Nextflow (${NF_OUT}) against sbatch (${OUTDIR})"
load_pixi

status=0

echo
echo ">>> prepared/ -- byte comparison"
for f in pairs.tsv grna_targets.tsv discovery_threshold.txt sim_input.rds sceptre_template.rds; do
  a="${PREPARED_DIR}/${f}"
  b="${NF_DIR}/prepared/${f}"
  if [[ ! -s "${b}" ]]; then
    printf '  %-26s SKIP (not produced by this run)\n' "${f}"
  elif cmp -s "${a}" "${b}"; then
    printf '  %-26s IDENTICAL\n' "${f}"
  else
    printf '  %-26s DIFFERS\n' "${f}"
    status=1
  fi
done

echo
echo ">>> splits/ -- byte comparison of every split"
if [[ -d "${NF_DIR}/splits" ]]; then
  n_diff=$(diff -rq "${SPLITS_DIR}" "${NF_DIR}/splits" 2>&1 | wc -l)
  n_files=$(find "${NF_DIR}/splits" -name 'split_*.tsv' | wc -l)
  if [[ "${n_diff}" -eq 0 ]]; then
    echo "  ${n_files} splits, all IDENTICAL"
  else
    echo "  ${n_diff} difference(s) across ${n_files} splits"
    diff -rq "${SPLITS_DIR}" "${NF_DIR}/splits" 2>&1 | head -5 | sed 's/^/    /'
    status=1
  fi
else
  echo "  SKIP (not produced by this run)"
fi

echo
echo ">>> null_precomputations.rds -- shared replicates"
nf_null="${NF_DIR}/prepared/null_precomputations.rds"
if [[ -s "${nf_null}" && -s "${NULL_MODELS}" ]]; then
  pixi_run Rscript workflow/compare_null_models.R \
    --a "${NULL_MODELS}" \
    --b "${nf_null}" | sed 's/^/  /'
  [[ "${PIPESTATUS[0]}" -ne 0 ]] && status=1
else
  echo "  SKIP (one side missing)"
fi

echo
echo ">>> sim/ -- per-replicate rows, over the replicates both runners covered"
#
# The two runners may simulate different numbers of replicates, so compare the overlap. It should be
# exact: seeds derive from (seed, target, rep, effect_size), which involves neither the runner nor
# how many replicates it went on to do, so replicate 1 of split 1 is the same draw either way.
ES="${EFFECT_SIZES[0]}"
nf_sim_dir="${NF_DIR}/sim/es${ES}"
sb_sim_dir="${SIM_DIR}/es${ES}"
if [[ -d "${nf_sim_dir}" && -d "${sb_sim_dir}" ]]; then
  n_checked=0
  n_bad=0
  for f in "${nf_sim_dir}"/split_*.tsv; do
    [[ -s "${f}" ]] || continue
    # split_0001_es0.15_rep0.tsv -> split_0001.tsv
    base="$(basename "${f}")"
    counterpart="${sb_sim_dir}/${base%%_es*}.tsv"
    [[ -s "${counterpart}" ]] || continue

    # Restrict both sides to the replicates the Nextflow run actually produced, then compare.
    max_rep=$(awk -F'\t' 'NR>1 && $12+0 > m { m = $12+0 } END { print m+0 }' "${f}")
    if ! diff -q \
        <(awk -F'\t' -v m="${max_rep}" 'NR==1 || $12+0 <= m' "${counterpart}" | sort) \
        <(sort "${f}") >/dev/null; then
      [[ "${n_bad}" -lt 3 ]] && echo "    DIFFERS: ${base} (vs $(basename "${counterpart}"), reps 1-${max_rep})"
      n_bad=$(( n_bad + 1 ))
    fi
    n_checked=$(( n_checked + 1 ))
  done
  if [[ "${n_checked}" -eq 0 ]]; then
    echo "  SKIP (no overlapping split outputs)"
  elif [[ "${n_bad}" -eq 0 ]]; then
    echo "  ${n_checked} split(s) compared, all IDENTICAL over the shared replicates"
  else
    echo "  ${n_bad} of ${n_checked} split(s) differ"
    status=1
  fi
else
  echo "  SKIP (one side missing)"
fi

echo
if [[ "${status}" -eq 0 ]]; then
  echo "VERIFY PASSED -- the two runners agree on everything compared."
else
  echo "VERIFY FAILED -- see above."
fi

report_usage
exit ${status}
