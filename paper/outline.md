# Outline

The flow, section by section. Placeholders only — no drafted prose. Figure placeholders are
`<!-- FIGURE n -->` and are itemised in `figures.md`.

---

## §0. Positioning — this is retrospective power, not experimental design

**State this in the abstract and again in the first paragraph of §2.** It is the single sentence that
separates this work from PerturbPlan (Niu et al. 2026), and a reader who misses it will think the two
papers do the same thing.

| | PerturbPlan | This work |
|---|---|---|
| When it is used | **before** the experiment | **after** the experiment |
| Question | "how should I build the screen?" | "given the screen that ran, what does this uncalled pair mean?" |
| Granularity | aggregate — power averaged over pairs of interest | **per element–gene pair** |
| Inputs | design parameters chosen by the user (cells, MOI, reads, gRNAs per target) | quantities **measured** in the completed screen (perturbed cells, expression, dispersion) |
| Output | power for a candidate design | minimum detectable effect size per pair, with an interval |
| Effect size | an assumption supplied as input ("minimum FC of interest") | the axis being solved for |

The two are complements, not competitors: one chooses the experiment, the other interprets it. Say so
explicitly and cite them as such — it is a stronger and more accurate framing than staking a claim to
novelty they would dispute, and their own paper leaves retrospective analysis out of scope.

### The objection this framing invites, and the answer

"Post hoc power" is a term with a bad reputation in statistics, and for a good reason: computing power
at the *observed* effect size is circular, because it is a deterministic function of the p-value and
therefore adds no information (Hoenig & Heisey 2001, *The Abuse of Power*). A reviewer will raise this.
**[TODO]** get that citation properly and confront it head-on in §2 rather than waiting for review.

The answer is that this is not observed power:

- Power is evaluated at **pre-specified** knockdown magnitudes — a fixed grid of effect sizes chosen
  before looking at any pair's result — not at whatever effect the data happened to show.
- The quantity reported is a **minimum detectable effect size**, which is a statement about the
  experiment's resolution for that pair, not a restatement of its p-value.
- It is therefore a *sensitivity* or *detectability* analysis in the design-analysis sense, computed
  after the fact because that is when the covariates are known — not an attempt to rescue a null result
  by recomputing its power.

Prefer the language of **detectability** and **minimum detectable effect size** over "post hoc power"
throughout, while naming the term once so readers can locate the work in the literature.

---

## 1. Why power matters for element–gene screens

**The setup.** Single-cell CRISPRi screens test element–gene pairs and return a p-value per pair.
A pair that is not called is recorded as a negative.

**The problem.** A negative has two possible causes that the screen output does not distinguish:

- the element does not regulate the gene, or
- the experiment could not have seen it if it did.

**Why that matters, concretely — two audiences:**

1. **Models trained on CRISPR data.** scE2G and ENCODE-rE2G learn from both positives and negatives.
   An underpowered pair labelled negative is a mislabelled training example. Worse than noise: the
   mislabelling is *not random* — it concentrates in low-expression genes and poorly perturbed
   elements, which is exactly the regime a regulatory model is trying to learn.
   **[NEEDS DATA]** quantify: what fraction of the negatives in a published training set are pairs
   that could not have reached power 0.8? On `day0_grna20_no_shuffle` at a 15 % knockdown, 49.9 % of
   tested pairs are clearly underpowered and 21.9 % of elements (663/3,026) have no tested pair that
   could have reached 0.8 at all (`docs/status.md`). If those proportions hold in a training set, the
   negative class is substantially contaminated.
2. **Any analysis conditioning on null results.** Comparing regulatory architecture across cell types,
   or claiming an element is inert, silently conditions on power. Where power correlates with a
   biological covariate — and it does, through expression — power becomes a confounder.

<!-- FIGURE 1: the problem. (a) distribution of per-pair power at a 15% knockdown, partitioned into
     certified / ambiguous / clearly underpowered. (b) the same pairs' power against gene expression,
     showing that power is not random with respect to biology. (c) schematic: identical screen output,
     two different underlying truths. -->

**Framing to keep.** The question is *false-negative triage*: for a pair CRISPRi did not call, is the
link absent, or could the experiment not have seen one?

---

## 2. What exists today

**Search done — see `literature.md`. The claim as originally framed does not survive, and §2 must be
rewritten around the narrower, true gap:** methods exist for *prospective, aggregate* power; none give
*retrospective, per-pair* detectability with quantified uncertainty.

Open with the positioning table from §0, then cover:

- **scPower** (Schmid et al.) — design-stage power for single-cell DE/eQTL studies. The comparison
  reviewers will demand. Establish precisely what it does and does not do: our reading is that it
  sizes an experiment prospectively rather than certifying an individual pair's negative
  retrospectively, but this needs verifying against the paper.
- sceptre's own treatment of calibration and power (Barry et al.).
- Power analyses reported in the primary screen papers — Gasperini et al., Ray et al. (DC-TAP), and
  the Perturb-seq lineage — which are typically dataset-specific and not released as reusable methods.
- Anything doing per-pair retrospective power for CRISPR screens. If it exists, the paper's
  contribution shifts from "first method" to "first *scalable and predictive* method", which is still
  a contribution but must be framed honestly.

Findings go in `literature.md`.

---

## 3. The method: simulation-based per-pair power

**What it does.** For a pair, impose a known knockdown on the gene in the element's perturbed cells,
re-run the real testing procedure, and count how often it is detected. Repeat over replicates.

Points to make:

- The simulation calls the **same** test the discovery analysis used, so the power estimate inherits
  the discovery analysis's multiple-testing correction rather than a nominal alpha.
- Counts are drawn per cell from that cell's own size factor and the gene's fitted dispersion.
- Seeds derive from `(seed, target, replicate, effect_size)`, so estimates are invariant to how the
  work is partitioned across a cluster. This is what makes the paired comparisons in §7 exact.

**Report the interval, not the point estimate.** Asserting "this negative is biological" asserts power
was *at least* 0.8, which is a claim about a lower bound. Certifying it needs 88/100 replicates — or
29/30, which is why 30 replicates cannot support the claim regardless of what a precision table says.
`power_ci_low` is the quantity to threshold; Wilson intervals, because the normal approximation gives
`[0, 0]` at zero successes.

**The deliverable is minimum detectable effect size, not a row of power columns.** Per pair: the
smallest tested effect size at which `power_ci_low` clears the threshold *and* every larger effect size
does too. Taking the first that clears lets noise bias every pair toward looking more detectable.

<!-- FIGURE 2: method schematic. sceptre object -> impose knockdown on perturbed cells -> simulate
     counts -> re-run the real test -> replicate -> power with interval -> minimum detectable effect
     size. Annotate where the null model enters, which §7 shows matters. -->

**Cost.** **[NEEDS DATA]** final numbers from the 100-replicate array. Current measured model per
(target, replicate): `5.51s + 0.891s x pairs` under fitted null models. Emphasise the structural point
— roughly half the cost is a per-call term paid regardless of how few pairs ride along — because it is
what makes §4 and §5 necessary rather than merely nice.

---

## 4. Making it affordable: three effect sizes instead of six

**The curve is derived, not fitted.** For a Wald-type test at a fixed threshold,

```
power(e) = Phi(beta / SE - z),   beta = -log(1 - e),   z = qnorm(1 - alpha)
```

`beta` is the effect on the scale the test works on, `SE` collects everything pair-specific, `z` is
fixed by the discovery threshold and shared across pairs. On the probit scale this is a straight line
through `-z` with slope `1/SE`: **one free parameter per pair**. Three effect sizes therefore leave two
degrees of freedom to *check* the fit rather than just enough to force it.

Fitted as a binomial GLM with a probit link, no intercept, `-z` as an offset — not least squares on
`qnorm(power)`, because 0/100 and 100/100 still carry information about the slope where `qnorm` is
infinite.

**Held-out validation.** Fit on three effect sizes, predict the other three, compare to measured.
Best grid on `dc_tap_paper_wtc11_no_shuf`: 0.05 / 0.1 / 0.25 → MAE 0.0167, p90 error 0.056, and the
six-point minimum detectable effect size reproduced exactly for 89.0 % of pairs.

**Honesty required here.** `deviance/df` is 2.4–2.7 (median), so the straight line is *not* the true
curve at the saturated ends. It interpolates the transition to within ~1.5x Monte-Carlo noise. Frame as
interpolation, not extrapolation.

<!-- FIGURE 3: the curve. (a) probit-scale linearity for example pairs, six measured points with the
     three-point fit overlaid. (b) held-out predicted vs measured power, all pairs, with the noise
     floor as a reference band. (c) MAE by fit grid. (d) deviance/df distribution, marking where the
     straight line fails. -->

---

## 5. How the three points are chosen

**The finding.** Grid *placement* matters more than the number of points. A one-parameter sigmoid is
only determined by points away from 0 and 1, so a grid that brackets the transition beats one that
samples it. Placement is dataset-specific: `wtc11` transitions between 5 % and 15 %, `day0` mostly
between 5 % and 25 %.

Projecting each `day0` pair's curve from its measured 0.15 power:

| Grid | `day0` pairs with ≥1 point where power is 0.1–0.9 |
|---|---:|
| 5 / 25 / 50 | 26.3 % |
| 5 / 15 / 50 | 60.4 % |
| **10 / 20 / 35** | **100.0 %** |
| 5 / 15 / 25 / 50 | 72.3 % |

**The recipe.** Pick the grid from a cheap pilot — one effect size at low replicate count locates the
transition distribution — rather than reusing another dataset's grid.

**[NEEDS DATA]** the pilot recipe is reasoned, not yet demonstrated end to end. Show: run the pilot,
choose the grid from it, run that grid, and confirm the resulting fit is as good as the best grid
chosen with hindsight.

<!-- FIGURE 4: grid choice. (a) transition-location distribution for two datasets, showing they differ.
     (b) fraction of pairs with an informative point, by grid. (c) pilot-chosen grid vs
     hindsight-optimal grid, held-out MAE. -->

---

## 6. The main result: power for pairs never simulated

**Why it matters.** sceptre can test an element against genes on other chromosomes (`trans`), which is
how a null distribution for element–gene testing is obtained — but simulating power pair by pair at
that scale is infeasible. If power can be predicted from covariates the pipeline already reports, it
becomes available wherever those covariates are, including pairs and experiments never simulated.

**The model.** Theory says `SE^2 ~ (1/n_pert_cells)(1/mu + 1/theta)`, so `k = 1/SE` should be
predictable. On 6,030 `wtc11` pairs:

```
log(k) = -0.605 + 0.502 * log(perturbed cells) + 0.351 * log(expression)
                        theory: 0.500              theory: 0 to 0.5
R^2 = 0.867     residual sd of log k = 0.243  ->  k to within x1.28
```

The perturbed-cell exponent landing on 0.502 against a theoretical 0.500 confirms the `sqrt(n)` law on
real data.

**The limitation, stated up front rather than buried.** x1.28 in `k` is ≈ ±0.20 in power mid-transition.
That is *model* error, so it does not shrink with more replicates, and it is unlikely to be spread
evenly across genes. A point prediction is therefore not a substitute for measurement where per-pair
certification is concerned.

**The fix, and a contribution in its own right.** Calibrate the residual quantiles on a measured subset
and widen each predicted interval by them. Honest intervals, wider than measured ones — and the width
quantifies what a measured sweep buys.

**[NEEDS DATA] — this section is the paper's load-bearing claim and is the least validated:**

- Cross-dataset transfer of the *coefficients*, not just per-dataset fit quality. Needs
  `dc_tap_paper_k562_no_shuf` and Gasperini et al.
- Where the residuals concentrate (expression decile? dispersion? cell count?), since non-random error
  is what harms a downstream classifier.
- Perturbed-cell count explains only 2.3 % of variance *here* because it barely varies; expression
  alone explains 83.7 %. A dataset with real variation in cell count is needed to show this is not
  effectively a one-covariate model.
- A direct `trans` demonstration: simulate a tractable subset of `trans` pairs and compare against the
  prediction. Without it, §6 extrapolates from `cis` to an unmeasured regime.

<!-- FIGURE 5: prediction. (a) fitted k against predicted k, with the sqrt(n) exponent annotated.
     (b) partial dependence on perturbed cells and on expression. (c) residuals by expression decile
     -- the non-randomness check. (d) calibrated prediction intervals vs measured power, held out. -->

<!-- FIGURE 6: the trans demonstration. predicted vs directly simulated power for a subset of trans
     pairs, with calibrated intervals. The panel that turns section 6 from argument into evidence. -->

---

## 7. What power means at different levels of aggregation

**Per pair** — "could we have detected this link?" This is the well-defined level and the one the
deliverable reports.

**Per element** — the tempting claim, "this element is well powered whatever gene you pair it with",
is not what the data say. At a 15 % knockdown only **20.9 %** of the variance in per-pair power is
between elements; ~79 % is gene-to-gene *within* an element. Mean within-element SD is **0.34**, and
element mean power correlates with `mean_pert_cells` at only 0.38. A mean over the pairs an element
happens to have tested is also not comparable between elements, since it inherits the expression levels
of whichever genes sit nearby.

The defensible version is conditional: per-element power answers "how well did we perturb this
element?", whose sufficient statistic is the perturbed-cell count, not anything about genes. Report it
at a reference gene: *"for a gene at median expression, element E has 0.7 power to detect a 15 %
knockdown."* Note this is exactly §6's model evaluated at a fixed expression — the aggregation
question and the prediction question have the same answer.

**What does hold at element level is a floor.** 21.9 % of elements have no tested pair that could have
reached power 0.8, so they cannot support a "regulates nothing" claim under any reading. Conversely
only 5 of 3,026 elements have every tested pair certified — element-wide negative claims are
essentially never assertable at 100 replicates and a 15 % knockdown.

<!-- FIGURE 7: aggregation. (a) variance decomposition, between- vs within-element. (b) within-element
     power spread for example elements, ordered by mean. (c) element power at a reference gene vs naive
     mean over tested pairs -- showing the naive version is confounded by neighbouring genes. -->

---

## 8. Validation and tests

The "solid tests" section. What must be in it:

- **Reproducibility.** Identical results across partitionings — 1x4 vs 2x2 replicate chunks, and across
  split layouts. Already verified; state it as a property of the seeding scheme.
- **Monotonicity.** Power must not decrease as effect size increases. **[NEEDS DATA]** at 100
  replicates across all pairs; on the six-point `wtc11` sweep only 220 of 32,870 consecutive
  comparisons decrease (0.67 %), largest 0.060 ≈ one standard error.
- **Held-out curve prediction** (§4), repeated on a second dataset. **[NEEDS DATA]**
- **Cross-dataset coefficient transfer** (§6). **[NEEDS DATA]**
- **A null-model trap worth its own subsection.** Building a power simulation on a sceptre object
  silently inherits `@response_precomputations` from the real discovery analysis, so simulated counts
  are tested against coefficients fitted to *real* counts. Measured: 7 of 265 discovery calls flip
  against a faithful refit, **all 7 in the same direction**, so power is systematically understated.
  Refitting inside every call is correct but 4.3x dearer; fitting once per (gene, replicate) on a null
  simulation is equivalent to refitting (0 flips of 265) at +1.7 %. Anyone reimplementing this will hit
  the same trap, which is why it belongs in the paper and not only in the repo.
- **Comparison against a reference implementation.** **[NEEDS DATA]** old vs new pipeline,
  distributional rather than exact since the RNG stream differs.

<!-- FIGURE 8 (supplementary): the null-model trap. (a) p-value agreement looks fine in aggregate
     (Spearman 0.999). (b) the same data at the call level, showing one-directional flips. The
     pedagogical point: correlation hides the thing that matters. -->

<!-- FIGURE 9 (supplementary): cost. measured runtime per call against pairs per call, the three null
     model configurations, and the implied CPU-hours per effect size. -->

---

## 9. Discussion

- What the method licenses: per-pair negative claims with a stated effect-size bound; power-aware
  training labels; power as a covariate rather than a confounder.
- What it does not: element-wide "regulates nothing" claims, and per-pair certification from prediction
  alone without calibrated intervals.
- Practical guidance: report `power_ci_low`; report minimum detectable effect size; choose the grid from
  a pilot; never inherit a null model.
- Limits: validated on N datasets (**[NEEDS DATA]** — currently 1); `trans` demonstrated on a subset;
  the linear-probit form fails at the saturated ends.

---

## Open questions to resolve before drafting

1. Literature search (§2) — does the "no satisfying method" claim survive?
2. Is the paper's headline §6 (prediction) or §3–4 (a cheap correct method)? Currently written as §6,
   which is also the least validated. If the validation does not land, §3–4 is a complete paper on its
   own and §6 becomes a discussion point.
3. Which datasets ship: `wtc11` plus `k562` is within-DC-TAP; Gasperini et al. is the cross-lab test.
4. Is the software a contribution? A tool paper needs the workflow engine and a test suite; a methods
   paper does not.
