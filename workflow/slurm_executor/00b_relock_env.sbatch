#!/bin/bash
#SBATCH --job-name=power_00b_relock
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Step 0b -- re-solve pixi.lock after pixi.toml's dependencies change.
#
# Only needed when [dependencies] changes, which in practice means when the sceptre pin moves to a
# version with new DESCRIPTION Imports (R CMD INSTALL enforces Imports, so they must be in the
# environment). Step 0 runs `pixi install --frozen` and deliberately fails if pixi.toml and
# pixi.lock have diverged, so this job has to run first and its result has to be committed.
#
# WHY THIS IS A JOB AND NOT A LOGIN-NODE COMMAND
#
# Two reasons, and the second is the one that actually forces it:
#
#   1. The solve is not free -- it reads thousands of repodata shards.
#   2. Sherlock login nodes cap open file descriptors at 256, hard (`ulimit -Hn`), and pixi's
#      shard reader exhausts that mid-solve:
#
#        × failed to read from file .../repodata/shards-v1/<hash>.msgpack:
#          No file descriptors available (os error 24)
#
#      The limit cannot be raised on a login node. Compute nodes allow far more, so the solve
#      simply works here.
#
# It needs outbound network for the repodata fetch, which Sherlock compute nodes have.
#
# Usage:  sbatch workflow/slurm_executor/00b_relock_env.sbatch
#         then commit the updated pixi.lock, then sbatch 00_setup_env.sbatch

_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

log_header "Step 0b: re-solve pixi.lock"

load_pixi
cd "${REPO_ROOT}"

echo "open file limit: $(ulimit -n) soft / $(ulimit -Hn) hard  (login nodes cap at 256, which is too few)"
echo

# Keep a copy, so a bad solve can be reverted without relying on git alone.
backup="${REPO_ROOT}/pixi.lock.before-relock"
cp "${REPO_ROOT}/pixi.lock" "${backup}"
echo "backup: ${backup}"

echo
echo ">>> pixi lock"
container_run pixi lock

## VERIFY ==========================================================================================

# The point of the re-solve is that `--frozen` works again afterwards. If this fails, the lock and
# the manifest are still out of step and step 0 would fail in exactly the same way.
echo
echo ">>> pixi install --frozen (default)"
container_run pixi install --frozen

echo
echo ">>> pixi install --frozen (build)"
container_run pixi install --frozen --environment build

echo
echo ">>> what changed"
if command -v git >/dev/null 2>&1; then
  git --no-pager diff --stat -- pixi.lock || true
  echo
  echo "packages added or removed:"
  diff <(grep -oE '^\s+- (conda|pypi): .*' "${backup}" | sort -u) \
       <(grep -oE '^\s+- (conda|pypi): .*' "${REPO_ROOT}/pixi.lock" | sort -u) \
    | grep -E '^[<>]' | sed 's/^/    /' || echo "    (none)"
fi

echo
echo "Re-locked. Commit pixi.toml and pixi.lock together, then run:"
echo "  sbatch workflow/slurm_executor/00_setup_env.sbatch"
report_usage
