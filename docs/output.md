---
title: Output
nav_order: 3
---

# Output

## `power_summary.tsv` — the main result

One row per element–gene pair, one set of columns per effect size.

| Column | Meaning |
|---|---|
| `grna_target` | Perturbation target (element), e.g. `chr2:201986246-201986547`. |
| `response_id` | Gene. |
| `mean_pert_cells` | Mean number of perturbed cells across replicates. |
| `average_expression_all_cells` | Raw mean expression of the gene across all cells. Not size-factor normalised — see the note below. |
| `PowerAtEffectSize15` | Power at a 15% knockdown. One column per effect size; the suffix is `effect_size × 100`. |
| `PowerAtEffectSize15_ci_low`, `_ci_high` | 95% Wilson interval for that estimate. |
| `PowerAtEffectSize15_n_reps` | Replicates contributing to it. |
| `min_detectable_effect_size` | Smallest **tested** effect size reaching `--power-threshold` (default 0.8). `NA` means no tested effect size did. |
| `max_effect_size_tested` | The largest effect size in the run, so `NA` above can be interpreted. |

`NA` in `min_detectable_effect_size` is a statement about the effect sizes you ran, **not** evidence
that a pair is undetectable. If you tested 0.15 and 0.2 and a pair needs 0.4, it will be `NA`.
That is why `max_effect_size_tested` sits next to it.

## `power_es<effect_size>.tsv` — per effect size

| Column | Meaning |
|---|---|
| `grna_target`, `response_id` | The pair. |
| `power` | Fraction of replicates in which sceptre would have called the association. |
| `power_ci_low`, `power_ci_high` | 95% Wilson score interval. |
| `n_reps` | Replicates contributing. Can be below `--reps` if any replicate produced no fold-change estimate. |
| `mean_log_2_fold_change` | Mean simulated log₂ fold change across replicates. |
| `mean_pert_cells` | Mean perturbed cells. |
| `average_expression_all_cells` | Raw mean expression. |
| `effect_size` | The effect size simulated. |

### What `power` actually counts

A replicate counts as a success when **both** hold:

```
p_value < threshold   AND   log_2_fold_change < 0
```

The second condition makes this one-sided: a replicate only counts if the simulated perturbation
*reduced* expression. `threshold` is the largest nominal p-value that survived multiple-testing
correction in the real discovery analysis, so "detected" means the same thing here as it did in your
actual results.

### Reading the confidence intervals

Always read `power` together with its interval, especially at the boundaries. `power = 0` does not
mean the pair is undetectable — at 100 replicates the Wilson upper bound is 0.037, so the data are
consistent with a true power of a few percent. Likewise `power = 1` has a lower bound of 0.963.

The interval is Wilson rather than `p̂ ± 1.96·SE` precisely because the normal approximation returns
`[0, 0]` for 0 successes, asserting certainty the data do not support. See
[Choosing num_replicates]({{ site.baseurl }}{% link choosing-num-replicates.md %}).

## `sim_*.tsv` — per-replicate detail

One row per (pair, replicate). Mostly useful for debugging or for re-deriving power with a
different threshold. Columns are sceptre's discovery-result columns (`response_id`, `grna_target`,
`n_nonzero_trt`, `n_nonzero_cntrl`, `pass_qc`, `p_value`, `fold_change`, `se_fold_change`,
`log_2_fold_change`, `significant`) plus:

| Column | Meaning |
|---|---|
| `num_pert_cells` | Perturbed cells for this target. |
| `rep` | Replicate index, unique across chunks thanks to `--rep-offset`. |
| `effect_size` | The effect size simulated. |
| `average_expression_all_cells` | Raw mean expression of the gene. |

Note that `n_nonzero_trt`, `n_nonzero_cntrl` and `significant` are carried over from the **real**
discovery pairs, not recomputed from the simulated data — they describe the observed experiment, and
are kept because sceptre requires them on the pair table.

## Intermediates

| File | Contents |
|---|---|
| `sim_input.rds` | Per-gene `mean`, `dispersion`, `average_expression_all_cells`; per-cell `size_factors` and categorical covariates; the gRNA and target perturbation matrices. No count matrix — the simulation draws counts rather than reading them. |
| `sceptre_template.rds` | The sceptre object with `@response_matrix` and `@grna_matrix` emptied. Neither is read by the discovery analysis once gRNAs are assigned. |
| `pairs.tsv` | `grna_target`, `response_id` for QC-passing pairs only. |
| `grna_targets.tsv` | `grna_id`, `grna_target`. |
| `discovery_threshold.txt` | A single number: the p-value a replicate must beat. |
| `split_*.tsv` | Subsets of `pairs.tsv`, balanced by pair count. |

## Two means, deliberately

`sim_input.rds` carries two per-gene means and they are not interchangeable:

- **`mean`** — size-factor normalised. This is what the simulation draws counts from.
- **`average_expression_all_cells`** — raw. Reported in the output so power can be related to
  expression level.

Only the raw one appears in the output tables.
