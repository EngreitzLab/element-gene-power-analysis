---
title: Choosing num_replicates
nav_order: 4
---

# Choosing `num_replicates`

`num_replicates` (`--reps` on `run_power_simulation.R`) is the number of Monte-Carlo simulation
replicates run per element–gene pair. It is the **only** pipeline parameter with statistical
content — `n_splits` and `reps_per_chunk` only decide how the work is divided across tasks.

This page explains what the number controls, how to pick it for your question, and how to read the
confidence intervals in the output.

---

## What is actually being estimated

For one pair and one effect size, each replicate simulates a fresh count matrix under that effect
size and asks sceptre whether it would have called the association. A replicate is a **success**
when

```
p_value < threshold   AND   log_2_fold_change < 0
```

The second condition makes this a one-sided criterion: a replicate only counts if the simulated
perturbation *reduced* expression, which is the direction a repressive element is expected to act
in. The reported power is the success fraction:

```
power = successes / n_reps
```

so it is a binomial proportion with `n_reps` trials, and every replicate is an independent draw.
`num_replicates` sets the precision of that proportion — nothing else. It does not make the
simulation more realistic, and it does not change the quantity being estimated.

Two consequences follow directly.

**Resolution.** The estimate can only take multiples of `1/n_reps`. At 20 replicates the only
achievable values are 0, 0.05, 0.10, …; a reported power of 0.85 is not representable at all.

**Monte-Carlo error.** The standard error of the estimate is `sqrt(p(1-p)/n_reps)`, largest when
the true power is near 0.5 and shrinking to zero at the extremes:

| `num_replicates` | SE at p=0.5 | 95 % CI width at p=0.5 | SE at p=0.8 | resolution |
|---:|---:|---:|---:|---:|
| 20 | 0.112 | ±0.219 | 0.089 | 0.050 |
| 50 | 0.071 | ±0.139 | 0.057 | 0.020 |
| 100 | 0.050 | ±0.098 | 0.040 | 0.010 |
| 200 | 0.035 | ±0.069 | 0.028 | 0.005 |
| 400 | 0.025 | ±0.049 | 0.020 | 0.003 |
| 1000 | 0.016 | ±0.031 | 0.013 | 0.001 |

---

## How many replicates do you need?

There is no single right answer, because it depends on what you do with the number. Find your case
below.

### Aggregate summaries — 100 is generous, even 20 would do

If you report a mean power across pairs, or an expected number of discoveries, the Monte-Carlo
error averages down across pairs as `1/sqrt(num_replicates × n_pairs)`:

| pairs summarised | SE of mean power at `num_replicates` = 100 |
|---:|---:|
| 1 | 0.050 |
| 100 | 0.005 |
| 32,386 | **0.0003** |

For a dataset-level number over tens of thousands of pairs, replicate count is not your limiting
source of uncertainty — the effect-size assumption is.

### Per-pair power in a results table — 100 is reasonable

At `num_replicates = 100` a per-pair estimate carries a 95 % interval of roughly ±0.10 in the worst
case. That is honest and usually adequate, provided the interval is reported alongside it (it is —
see below). At 20 replicates the interval is ±0.22, which is too coarse to publish per pair.

### Per-pair decisions at a cutoff — you need ~400

This is the case where 100 replicates is genuinely weak. If you filter pairs on something like
"power ≥ 0.8", the standard error at p=0.8 is 0.040, so the boundary is blurred by about ±0.08:
pairs whose true power is 0.72 and 0.88 are routinely swapped across the cutoff. Raising to 400
halves that to ±0.039.

Before reaching for 400 uniformly, read the next section — you almost certainly do not need it for
every pair.

### Two-stage allocation: 400-level precision at roughly a quarter of the cost

The variance `p(1-p)/n_reps` collapses at the extremes, and most pairs *are* at the extremes.
Measured on the `sample1` dataset (33 pairs, 30 replicates, effect size 0.15):

- **55 %** of pairs came out at exactly 0 or 1
- only **21 %** fell in the uncertain band (0.1, 0.9)

Extra replicates spent on a pair that is 0/30 or 30/30 buy almost nothing. So:

1. **Stage 1** — run `num_replicates = 30` for every pair.
2. **Stage 2** — top up to 400 only for pairs whose stage-1 estimate lies in (0.1, 0.9).

Average cost is about `30 + 0.21 × 370 ≈ 108` replicate-equivalents per pair, versus 400 for the
uniform design: **~3.7× cheaper for the same precision where precision matters.**

Pooling the two stages is unbiased — replicates are exchangeable draws, so combining them is fine
as long as `n_reps` is reported per pair, which the output does. The one thing to avoid is
comparing raw power values between pairs while ignoring that they rest on different replicate
counts; use the intervals for that.

To run it, invoke the second stage with `--rep-offset 30` on the subset of pairs that need it, so
the `rep` column stays unique and the two sets of results can simply be concatenated.

---

## Reading the confidence intervals

`compute_power.R` reports, per pair:

| column | meaning |
|---|---|
| `power` | success fraction, `successes / n_reps` |
| `n_reps` | replicates actually contributing (replicates with a missing `log_2_fold_change` are excluded, so this can be below `num_replicates`) |
| `power_ci_low`, `power_ci_high` | 95 % **Wilson** score interval |

The interval is Wilson rather than the textbook `p̂ ± 1.96·SE` for a specific reason: the normal
approximation breaks down exactly where these estimates live. A pair with 0 successes out of 100
gets the interval `[0, 0]` under the normal approximation — asserting with certainty that the power
is zero. The Wilson interval gives `[0, 0.037]`, which is what the data actually support: 100
replicates without a single success is consistent with a true power of a few percent.

So **do not read `power = 0` as "this pair cannot be detected"**. Read it together with
`power_ci_high`. The same applies at the top: `power = 1` at 100 replicates has a Wilson lower
bound of 0.963, not 1.

Worked examples at `num_replicates = 100`:

| successes | `power` | Wilson 95 % interval |
|---:|---:|---|
| 0 | 0.00 | [0.000, 0.037] |
| 5 | 0.05 | [0.022, 0.112] |
| 50 | 0.50 | [0.404, 0.596] |
| 80 | 0.80 | [0.711, 0.867] |
| 100 | 1.00 | [0.963, 1.000] |

Note how the 0.8 row overlaps a cutoff of 0.75 *and* is consistent with 0.86 — that is the
cutoff-blurring described above, made explicit.

---

## Cost

Replicate count is the parameter you pay for linearly. Measured on `sample1` (2,798 targets,
32,386 pairs, all control cells, one `run_discovery_analysis()` call ≈ 3.4–4.8s):

| `num_replicates` | CPU-hours per effect size |
|---:|---:|
| 30 | ~89 |
| 100 | ~295 |
| 400 | ~1,180 |
| two-stage to 400 | ~320 |

Doubling replicates doubles cost and only improves precision by `sqrt(2)` — which is the whole
argument for spending them where the variance is, rather than uniformly.

---

## Reproducibility

Replicate `k` of a given pair always receives the same random draws, because the seed is derived
from `(seed, target, rep, effect_size)` rather than from the task. Two practical consequences:

- Results are **invariant to `n_splits` and `reps_per_chunk`**. Re-running the same analysis with a
  different task layout reproduces the same numbers.
- Increasing `num_replicates` is **incremental**: replicates 1–100 keep exactly the values they had,
  and 101–400 are new. You can extend an existing run rather than redo it.

Both properties depend on passing the same `--seed`. It is required, not optional — the estimates
are stochastic, and an unseeded run cannot be reproduced or extended.

## See also

- `docs/methods.md` — how the effect size and guide-to-guide variability are parameterised
- `docs/output.md` — every output column
- `docs/usage.md` — `n_splits` / `reps_per_chunk` and SLURM sizing
