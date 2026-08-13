#!/bin/bash
#SBATCH --job-name=power_00_setup
#SBATCH --partition=engreitz,owners
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/refactor/%x-%A.out
#SBATCH --error=logs/refactor/%x-%A.err
#
# Step 0 -- build the pixi environment and install the pinned sceptre.
#
# Not a pipeline step, but a hard prerequisite: nothing else in this directory can run until
# .pixi/ exists and sceptre is installed into .pixi/rlibs.
#
# This is a job rather than a login-node command because `pixi run setup` compiles sceptre from
# source, and compiling on a shared login node is not acceptable. It does need outbound network
# access, for both the conda solve and the GitHub checkout -- Sherlock compute nodes have it, but
# if this job fails on a network error, the fallback is to run `pixi install` on a login node
# (download only, no compilation) and then resubmit this job for the compile.
#
# Usage:  sbatch workflow/slurm_executor/00_setup_env.sbatch

# sbatch copies the script into a spool directory, so $0 is not the path that was submitted.
# Ask Slurm where the script came from, and fall back to $0 when run directly with bash.
_cmd="$(scontrol show job "${SLURM_JOB_ID:-0}" 2>/dev/null | sed -n 's/.*Command=\([^[:space:]]*\).*/\1/p' | head -1)"
source "$(cd "$(dirname "${_cmd:-$0}")" && pwd)/config.sh"

log_header "Step 0: build pixi environment and install pinned sceptre"

load_pixi
cd "${REPO_ROOT}"

echo "PIXI_CACHE_DIR:  ${PIXI_CACHE_DIR}"
echo "repo root:       ${REPO_ROOT}"
echo "host glibc:      $(ldd --version | head -1 | awk '{print $NF}')  (too old for the lock; hence the container)"
echo "container glibc: $(container_run ldd --version | head -1 | awk '{print $NF}')"
echo

## SOLVE / DOWNLOAD ===============================================================================

# --frozen: install exactly what pixi.lock specifies. If this errors saying the lockfile is out
# of date, that is a real signal -- pixi.toml and pixi.lock have diverged and someone needs to
# re-lock deliberately, rather than have a compute node do it silently.
echo ">>> pixi install (default environment)"
container_run pixi install --frozen

echo
echo ">>> pixi install (build environment; carries the compiler toolchain and r-remotes)"
container_run pixi install --frozen --environment build

## INSTALL SCEPTRE ================================================================================

# `setup` is scoped to the build feature in pixi.toml, so it resolves to the build environment
# without naming it. It installs sceptre from the SHA pinned in [activation.env] into
# .pixi/rlibs, and is a no-op if that exact SHA is already present.
echo
echo ">>> pixi run setup (compile sceptre from the pinned commit)"
container_run pixi run --frozen setup

## VERIFY =========================================================================================

# The gate. The pipeline reads several *unexported* S4 slots of sceptre_object; this check is
# what proves the pinned version still exposes them. If it fails, nothing downstream is
# trustworthy, so fail the job here.
echo
echo ">>> pixi run check-api"
pixi_run check-api

echo
echo "Environment ready:"
echo "  R:       $(pixi_run Rscript -e 'cat(R.version.string)' 2>/dev/null)"
echo "  sceptre: $(pixi_run Rscript -e 'cat(as.character(utils::packageVersion("sceptre")))' 2>/dev/null)"
echo "  rlibs:   ${REPO_ROOT}/.pixi/rlibs"

report_usage
