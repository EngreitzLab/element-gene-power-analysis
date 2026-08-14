#!/bin/bash
#SBATCH --job-name=power_12_tests
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Run the unit tests.
#
# A job rather than a login-node command because it starts R, which does not belong on a shared
# login node -- and because the tests must run inside the container for the same glibc reason
# everything else does.
#
# The suite covers the pure functions only: no sceptre, no lab data, nothing that needs the cluster.
# It is fast and it is the part that can regress silently, which is why it exists -- every function
# tested here has a corresponding entry in docs/status.md under "Correctness fixes".
#
# Usage:  sbatch workflow/slurm_executor/12_unit_tests.sbatch

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

log_header "Step 12: unit tests"

load_pixi

# testthat is in the `dev` feature, not the default environment, so this runs there. Invoking the
# `test` task rather than an inline Rscript keeps this script and `pixi run test` from drifting --
# there is one definition of what "the tests" are, in pixi.toml.
set +e
pixi_run_env dev test
status=$?
set -e

echo
if [[ "${status}" -eq 0 ]]; then
  echo "TESTS PASSED"
else
  echo "TESTS FAILED"
fi

report_usage
exit ${status}
