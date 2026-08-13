---
title: Status and handoff
nav_order: 7
---

# Status and handoff

State of the refactor, written so the work can be picked up by someone else.

**Short version:** the five pipeline steps are complete, run standalone, and have now run on
Sherlock as well as on a laptop. `workflow/slurm_executor/` holds one plain `sbatch` script per step
— deliberately no orchestration layer, so a failure is attributable to the step rather than the
runner — and steps 0–3 (environment, prepare, split, smoke test) are validated there. The full
simulation array at effect size 0.15 went out on 2026-08-13 as job `38849611`, 1,000 tasks at 100
replicates.

Still missing: Nextflow (Step 7), a test suite (Step 8), and the old-vs-new comparison, which is
scripted in `workflow/compare_old_new.R` but cannot run until the array lands. The comparison is now
straightforward because the reference sample `day0_grna20_no_shuffle` was produced by the old
pipeline with the size-factor shuffle already removed, so only the RNG stream differs between the two
implementations — see the note in `bin/lib/simulate.R`.

---

## Done

### The pipeline

Five standalone executables in `bin/`, each with `--help`, none referencing `snakemake@`:

```
prepare_sim_input.R -> split_pairs.R -> run_power_simulation.R -> compute_power.R -> summarize_power.R
```

Shared code lives in `bin/lib/` (`cli.R`, `sim_input.R`, `pert_input.R`, `simulate.R`,
`sceptre_io.R`), sourced by a script-relative path so it works both standalone and staged onto
`PATH` by a workflow engine.

### Measured improvements

| | Before | After | How verified |
|---|---|---|---|
| Input files required | 4 | **1** | the other three duplicate `@discovery_pairs_with_info`, `@grna_target_data_frame`, `@discovery_result` |
| Per-task read | 427 MB / 3.3s | **59 MB / 1.1s** | timed `readRDS` on both |
| Per-task RAM | 1,774 MB+ | ~250 MB | `object.size` / peak heap |
| Dead code | — | **694 lines removed** | 6 files, none reachable |
| Split imbalance | 1.18× | **1.000** | measured on 2,798 targets |
| Reproducibility | none | seeded and layout-invariant | identical p-values across 1×4 vs 2×2 replicate chunks, and across split layouts |
| Environment | 273-line pinned linux-64 conda | 11 direct specs, 2 platforms | `pixi.lock` solves for `linux-64` + `osx-arm64` |

### Correctness fixes

| Issue | Consequence | Status on `sample1` |
|---|---|---|
| Dispersions in a list column with `NULL` holes; `unlist()` shortened the vector and `rnbinom` recycled it | every gene after the first gap simulated with another gene's dispersion | latent — 273 values where 292 expected, but all 239 *tested* genes had precomputations, so it never fired. Now a hard error |
| `max()` over possibly-empty significant discoveries → `-Inf` | zero power reported for every pair, silently | latent — 252 significant discoveries present. Now validated |
| `group_by(rep) %>% group_by(...)` | first grouping discarded | dead code, removed |
| Discovery-results file required but produced by no rule | undeclared fourth manual input | removed; read from the object |
| Sampled-controls path reordered covariates independently of matrix columns | covariates misaligned with counts | latent — the path never ran. Fixed |
| No `set.seed()` anywhere | no run reproducible | fixed; seeds derive from `(seed, target, rep, effect_size)` |
| Per-task seeding | `num_batches` silently changed results | fixed and verified invariant |
| `effect_label` trimmed trailing zeros unconditionally | `0.2` → column `PowerAtEffectSize2` | fixed; also renamed to `power_at_effect_size_20` |
| `sparseMatrix` local shadowing `Matrix::sparseMatrix` | fragile | gone with the old script |
| Size factors permuted across cells before simulating | each cell's library size was paired with a different cell's perturbation status, since `effect_size[i, j]` is indexed by cell | removed; each cell keeps its own size factor. This changes the RNG stream, so a seed no longer reproduces pre-refactor output draw for draw — see [Methods]({{ site.baseurl }}{% link methods.md %}) |

### Two of my own claims that measurement overturned

- **Control-cell sampling is not a usable speed lever.** It costs 21–60 % of power. Off by default.
  Table in [Methods]({{ site.baseurl }}{% link methods.md %}).
- **The speedup from smaller matrices was ~7×, not the ~109× the size reduction suggested.**
  Per-replicate cost is sub-linear in cell count; ~0.3s per replicate is fixed overhead.

### Infrastructure

- `pixi.toml` / `pixi.lock` — three environments (runtime / build / dev), sceptre pinned to
  `v0.10.3` (`9faa4373`), guarded by `pixi run check-api`
- `.githooks/pre-commit` — auto-fixes line endings, tabs (2 spaces for R/YAML, 4 otherwise),
  trailing whitespace, final newline; rejects files >512 KB, non-snake_case filenames, camelCase R
  identifiers. `pixi run lint` is the read-only repo-wide version
- Documentation published to GitHub Pages
- `config/config.yml` — parameters only; dead `reps` key removed

---

## Reference numbers

Needed to size cluster runs. Measured, not estimated. The right-hand column is Sherlock, on
`day0_grna20_no_shuffle`, and supersedes the laptop figures wherever they differ.

| Quantity | Laptop (`sample1`) | **Sherlock** |
|---|---|---|
| Genes × cells | 292 × 586,309 | same dataset shape |
| Response matrix | 1,095 MB in RAM, 55.9 % dense | |
| `cells_in_use` | 567,690 of 586,309 | |
| QC-passing pairs / targets | 32,386 / 2,798 | **34,886 / 3,026** (median 9 pairs/target, max 36) |
| Perturbed cells per pair | median 396 (IQR 195–542, max 1,927) | |
| Discovery p-value threshold | 7.71404 × 10⁻⁴ | |
| One `run_discovery_analysis()` call | 3.4–4.8s | **~9.7s** |
| Cost model, per (target, replicate) | — | **4.91s + 0.459s × pairs** (R² = 0.968, 36 targets) |
| **Total at 100 replicates** | ~295 CPU-h per effect size | **~858 CPU-h per effect size** |
| Peak RAM, simulation task | 1.5 GB → request 4 GB | ~1.9 GB heap → 8 GB requested |
| Peak RAM, `prepare_sim_input.R` | 7.7 GB → request 12 GB | 16 GB requested, ~70s |

**The ~295 CPU-hour figure was a laptop measurement and is roughly 3× optimistic.** Sherlock's
cores are slower per-thread than an M-series laptop, and sceptre's per-call cost is untouched by the
refactor, so it dominates. Sized against the old pipeline's own `sacct` record for its 2026-05-14
run, the refactor is still ahead on both axes:

| | Old pipeline | Refactored |
|---|---|---|
| Tasks | 1,004 | 1,000 |
| Mean per task | 78.2 min (range 0.8–134.5) | ~50 min |
| CPU-hours per effect size | 1,308 | **~814–858** |
| Wall clock | **14.9 h** (Snakemake, throttled) | **~50 min** (unthrottled array) |

Both ran 100 replicates on the same dataset, so this is like for like. The CPU win is 1.6×; the
wall-clock win is mostly the unthrottled array rather than the code.

**48 % of the cost is the per-target term** — one `run_discovery_analysis()` call carries the full
586,309-cell bookkeeping however few gene pairs ride along, and targets cannot be merged because
perturbation status differs per target. Consequences worth internalising before optimising:

- Halving the *pairs* saves ~26 %, not 50 %. Halving the *replicates* saves exactly 50 %.
- Filtering pairs for a second pass is much dearer than it looks: the 4,322 pairs whose interval
  straddles power 0.8 (12.4 % of pairs) are spread across 1,820 of 3,026 targets, so a top-up run
  pays 60 % of the per-target overhead to redo 12 % of the pairs — 82 % of its cost is overhead.
- Effect sizes cost the same as each other. A six-point sweep is six times the table above
  (~5,100 CPU-h), but the 2,900-task QOS budget still allows one wave: 480 splits × 6 effect sizes
  = 2,880 tasks, ~1h50m each, so ~2 h wall clock.

---

## Left to do

### Step 7 — Nextflow + SLURM profile

Not started. `config/config.yml` already holds the parameters it will consume
(`n_splits`, `reps_per_chunk`, `effect_sizes`, `num_replicates`, `seed`, …).

Planned shape:

```
samplesheet -> PREPARE_SIM_INPUT (per sample)
            -> SPLIT_PAIRS -> flatten
            -> combine(effect_sizes) x combine(rep_chunks)
            -> POWER_SIMULATION (fan-out)
            -> collectFile by (sample, effect_size)
            -> COMPUTE_POWER -> SUMMARIZE_POWER
```

Notes for whoever builds it:

- The combine step needs no script: `collectFile(keepHeader: true, skip: 1)` replaces the deleted
  `combine_sceptre_power_analysis.R`.
- Resources from the table above, not the old Snakemake guesses (which asked 8 GB for simulation and
  64 GB for prepare).
- `cpus 1` on the simulation process — sceptre runs with `parallel = FALSE` and parallelism comes
  from the fan-out.
- Removing `Snakefile`, `rules/` and `R/` belongs to this step. That also turns `pixi run lint`
  green: it currently flags exactly those legacy files.

### Step 8 — tests and comparisons

1. **Unit tests** (`testthat`, `pixi run test`) for the pure functions: `center_effect_size_matrix`
   centring, `create_effect_size_matrix` shape, split binning conservation, `wilson_interval`
   against known values, `build_dispersion_vector` erroring on a missing gene, `effect_label`
   across `0.15 / 0.2 / 0.5 / 0.125`.
2. **Synthetic test data generator** (`bin/make_test_data.R`) so CI can run without lab data.
3. **Old-vs-new comparison** — the main outstanding item; see below.
4. **Two-stage replicate allocation** — documented in
   [Choosing num_replicates]({{ site.baseurl }}{% link choosing-num-replicates.md %}) but not
   orchestrated. The scripts already support it through `--rep-offset`.

### Step 10 — re-derive the environment

The dependency list was written before the code was finished. Re-audit the finished scripts'
`library()` and `::` calls and trim `pixi.toml` to match.

---

## Running the comparisons on the cluster

### Setup

```sh
git clone https://github.com/EngreitzLab/element-gene-power-analysis.git
cd element-gene-power-analysis
pixi install
pixi run setup        # compiles sceptre; needs a compiler and network
pixi run check-api    # must print "All checks passed"
```

`pixi run setup` installs sceptre into `.pixi/rlibs` inside the checkout, so the whole environment is
self-contained and needs no module system. It does require outbound network access to GitHub — if the
compute nodes are offline, run `setup` on a login node first.

### Smoke test before committing CPU-days

Three targets, 2 replicates — should finish in a couple of minutes and proves the environment works.

```sh
mkdir -p prepared splits sim
Rscript bin/prepare_sim_input.R \
  --sceptre-object results/sample1/sceptre_object.rds --outdir prepared/

Rscript bin/split_pairs.R --pairs prepared/pairs.tsv --n-splits 280 --outdir splits/

Rscript bin/run_power_simulation.R \
  --sim-input prepared/sim_input.rds --sceptre-template prepared/sceptre_template.rds \
  --pairs splits/split_001.tsv --grna-targets prepared/grna_targets.tsv \
  --effect-size 0.15 --reps 2 --seed 20250812 --out sim/smoke.tsv
```

### Full run as a SLURM array

`prepare_sim_input.R` once, then one array task per (split, effect size).

```sh
#!/usr/bin/env bash
#SBATCH --job-name=power-prepare
#SBATCH --mem=12G
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
pixi run Rscript bin/prepare_sim_input.R \
    --sceptre-object results/sample1/sceptre_object.rds \
    --outdir prepared/
pixi run Rscript bin/split_pairs.R \
    --pairs prepared/pairs.tsv --n-splits 280 --outdir splits/
```

```sh
#!/usr/bin/env bash
#SBATCH --job-name=power-sim
#SBATCH --array=1-280
#SBATCH --mem=4G
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/sim_%A_%a.out

# Selected by position rather than by formatting the index: split_pairs.R pads the filename to
# the width of --n-splits, so 280 splits give split_001.tsv but 1000 would give split_0001.tsv.
SPLIT=$(ls splits/split_*.tsv | sed -n "${SLURM_ARRAY_TASK_ID}p")
for ES in 0.15 0.2; do
    pixi run Rscript bin/run_power_simulation.R \
        --sim-input prepared/sim_input.rds \
        --sceptre-template prepared/sceptre_template.rds \
        --pairs "$SPLIT" \
        --grna-targets prepared/grna_targets.tsv \
        --effect-size "$ES" --reps 100 --seed 20250812 \
        --out "sim/split_${SLURM_ARRAY_TASK_ID}_es${ES}.tsv"
done
```

Then per effect size:

```sh
for ES in 0.15 0.2; do
    pixi run Rscript bin/compute_power.R \
        --simulations "$(ls sim/*_es${ES}.tsv | paste -sd, -)" \
        --threshold-file prepared/discovery_threshold.txt \
        --out "power_es${ES}.tsv"
done

pixi run Rscript bin/summarize_power.R \
    --power power_es0.15.tsv,power_es0.2.tsv --out power_summary.tsv
```

Expect roughly an hour per array task at 10 targets × 100 replicates × 2 effect sizes. Time limit is
set to 3h for headroom — per-target cost varies about 1.5× with pair count.

### Comparison 1 — old vs new pipeline (the important one)

**Why it is outstanding:** the upstream statistics are proven bit-identical (size factors, normalised
means, raw means: max absolute difference exactly 0), and dropping `@grna_matrix` was proven to leave
discovery results unchanged. What has *not* been checked is the simulation itself, because seeding
moved inside the replicate loop, which changes the RNG stream. So the comparison has to be
distributional rather than exact.

The old code needs the old environment, which is why the cluster is the right place: the old
`envs/power_analysis.yml` is a fully-pinned `linux-64` solve. It cannot be solved on macOS, but
should solve on a Linux compute node.

```sh
# 1. check out the pre-refactor tree alongside the current one.
#    b3cf7ab, not c7448a7: it is the tip of the old `low-mem` branch and includes the
#    "Improve memory usage in power simulations" commit, so it is the version actually
#    being run in practice.
git worktree add ../egpa-old b3cf7ab
cd ../egpa-old

# 2. the old pipeline expects four inputs under results/<sample>/ with these exact names
mkdir -p results/sample1
cp /path/to/sceptre_object.rds        results/sample1/final_sceptre_object.rds
cp /path/to/grna_groups_table.rds     results/sample1/
cp /path/to/gene_grna_group_pairs.rds results/sample1/

# the fourth input has to be an .rds, so convert the discovery-results CSV:
Rscript -e 'saveRDS(read.csv("/path/to/sceptre_discovery_results.csv"),
                    "results/sample1/results_run_discovery_analysis.rds")'

# 3. run the old pipeline for the power results, NOT `all`
#    (`all` stops at the concatenated per-replicate file)
snakemake --use-conda --conda-frontend conda --cores 8 \
    results/sample1/power_analysis/power_analysis_results_es_0.15.tsv
```

Set `num_replicates` in that tree's `config/config.yml` to whatever you use for the new run so the
two are comparable, and remember the old `config/config.yml` also has the dead `reps: 20` key — it is
`num_replicates` that takes effect.

Then compare against the new `power_es0.15.tsv`. Suggested acceptance criteria:

| Check | Expectation |
|---|---|
| Pairs present | identical set |
| Per-pair power difference | centred on 0; no systematic offset (paired sign test not significant) |
| Mean power | agrees within the Monte-Carlo error of the replicate count used |
| Correlation of per-pair power | high (> 0.95 at 100 replicates) |

A systematic offset in either direction is a real discrepancy and worth chasing. Scatter across
individual pairs is expected — at 100 replicates the standard error of a single estimate is up to
0.05, so differences of ±0.1 on individual pairs are normal.

Use the same `--reps` in both, and note the old pipeline cannot be seeded, so it produces a different
draw every run. Running it two or three times gives a sense of its own run-to-run spread, which is
the yardstick for judging the old-vs-new difference.

### Comparison 2 — monotonicity at full replicate count

Power must not decrease as effect size increases. At 12 replicates on 33 pairs this held (mean power
0.356 → 0.604 → 0.838 for 0.15 / 0.25 / 0.5) with one violation of 0.08 that a two-proportion test
could not distinguish from noise. Worth repeating at 100 replicates across all 32,386 pairs, where
noise-driven violations should nearly vanish. A violation that survives is a bug.

### Comparison 3 — timing and memory at scale

Run with `--with-report`-style accounting (or just `sacct`) and record actual per-task wall time and
`MaxRSS`, so the resource requests above can be tightened. The 4 GB figure comes from a laptop run of
three targets; confirm it holds for splits containing the largest targets (up to 36 pairs).

---

## Decisions already made

Recorded so they are not relitigated:

| Decision | Rationale |
|---|---|
| Power-analysis only; no differential-expression step | the three orphan diffex scripts were unwired duplicates |
| Only the sceptre object as input | the other three files duplicate data inside it |
| No Bioconductor | `SingleCellExperiment` was used purely as a container |
| sceptre pinned to a commit, installed by `pixi run setup` | not on conda-forge/bioconda, and the pipeline reads unexported slots |
| `n_control_cells` off by default | measured 21–60 % power loss |
| Threshold derived from `@discovery_result`, not `alpha` | reflects the correction actually applied |
| Wilson intervals reported per pair | the normal approximation gives `[0, 0]` for zero successes |
| Prefer more `n_splits` over more `reps_per_chunk` | chunking replicates re-pays the per-target setup |
| local + slurm profiles only | no cloud executor needed |

## The analysis question this is all for

Recorded because it determines which columns matter, and therefore how many replicates to buy.

The goal is **false-negative triage**: for a pair CRISPRi did not call, is the link absent,
or could the experiment not have seen one? Three consequences follow, and
[Output]({{ site.baseurl }}{% link output.md %}#interpreting-negatives) explains each in full.

**1. Threshold `power_ci_low`, never `power`.** "This negative is biological" asserts that power was
*at least* 0.8. Certifying that needs 88/100 successes — or 29/30, which is why 30 replicates is
unusable for this question regardless of what the precision tables say. At 0.15 and 100 replicates:
37.7 % of pairs certified, 12.4 % ambiguous, 49.9 % clearly underpowered (that last half is not a
replicate problem — it needs more perturbed cells or a larger effect).

**2. Minimum detectable effect size is the deliverable, not six power columns.** `summarize_power.R`
now emits `min_detectable_effect_size` plus a `_ci_low` / `_ci_high` bracket, where `_ci_high` is
derived from `power_ci_low`: the interval inverts, as power rises with effect size. It turns the
sweep into one sentence per pair: *"we would have caught a knockdown of ≥ 25 %, so the absence of a
call rules out effects that large."* Qualifying requires clearing the threshold at that effect size
**and every larger one tested**; taking the first effect size that clears lets noise bias
every pair towards looking more detectable than it is.

**3. Per-element power needs reframing before it can be reported.** The tempting version — "this
element is well powered whatever gene you pair it with" — is not what the data say. At 0.15 only
**20.9 %** of the variance in per-pair power is between elements; 79 % is gene to gene within an
element, the mean within-element SD is **0.34**, and element mean power correlates with
`mean_pert_cells` at only 0.38. A mean over the pairs an element happens to have tested is also
incomparable between elements: it inherits the expression levels of whichever genes sit nearby.

The defensible framing is conditional: per-pair power answers *"could we have detected this link?"*,
per-element power answers *"how well did we perturb this element?"* — whose sufficient statistic is
the perturbed-cell count, not anything about genes. Report it at a reference gene: "for a gene at
median expression, element E has 0.7 power at a 15 % knockdown."

**Not implemented.** The fit is cheap and needs no new simulation — `power_summary.tsv` already
carries `mean_pert_cells` and `average_expression_all_cells` per pair, so fitting
`power ~ f(pert_cells, expression)` and evaluating at reference expression is post-processing.
Simulating a reference-gene panel instead would cost ~450 CPU-h, because the per-target term
is paid whether a target carries one gene or thirty. Whoever picks this up should decide where it
lives: a new `bin/` script, or a section of `summarize_power.R`.

What *does* hold at element level is a floor: **21.9 % of elements (663/3,026) have no tested pair
that could have reached power 0.8**, so they cannot support a "regulates nothing" claim under any
reading, and they should be excluded from biological interpretation rather than reported as null
results. Conversely only **5 of 3,026** elements have every tested pair certified — element-wide
negative claims are essentially never assertable at 100 replicates and a 15 % knockdown.

## Open questions

- **Two-stage replicate allocation**: worth wiring in? Costed properly it is ~1,768 CPU-h per effect
  size against 3,430 for a uniform 400 and 858 for a uniform 100 — a 1.9× saving over uniform 400,
  not the 3.7× an earlier version of this document claimed. That claim came from counting replicate
  equivalents, which ignores the per-target term; see
  [Choosing num_replicates]({{ site.baseurl }}{% link choosing-num-replicates.md %}). It only pays
  off if per-pair cutoff decisions are being made — which, per the section above, they are.
- **Is the 5 % effect size worth 100 replicates?** Nearly every pair is underpowered there, so the
  replicates buy a precise estimate of a number that is 0. Sweeping 5 % at 30 reps and spending
  the difference on 20–25 % would carry more information. Untested.
- **`--target-overhead` in `split_pairs.R`** defaults to 0, balancing purely on pairs. It should
  now be set from the measured cost model: the per-target term is 4.91s against 0.459s per pair, so
  `--target-overhead` ≈ 10.7 pair-equivalents. Splits are currently balanced on pairs alone, which
  is why per-task times still vary with how many targets a split happens to collect.
- **ODM / out-of-core support** is designed for but not implemented. The seam is
  `get_response_matrix()` in `bin/lib/sceptre_io.R`, which currently errors explicitly on an `odm`.
  The hard part is that poscounts size factors need a per-cell median over genes — a column-wise
  reduction — while an `odm` is row-accessible.
