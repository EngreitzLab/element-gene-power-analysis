# Literature

Prior art for the claim that no satisfying method exists to quantify power for element–gene pairs.

**Headline: the claim as written in `outline.md` §2 does not survive, and the reason is close to home.**
Two tools already address power for single-cell CRISPR screens, and one of them is from this lab. The
gap is narrower and more specific than "there is no method" — but it is real, and it is explicitly
signposted as future work by the paper that closes the wider gap.

---

## PerturbPlan — the paper that must be cited first

Niu Z, He Y, Galante J, Gschwind AR, Ray J, Steinmetz LM, **Engreitz JM**, Katsevich E.
*PerturbPlan: An analytical framework for designing Perturb-seq experiments.* bioRxiv, 23 May 2026.
PMID 42239142, PMC13228452, [DOI](https://doi.org/10.64898/2026.05.22.727199).
(Retrieved from PubMed.)

Note the author list: **Jesse Engreitz and Judhajeet Ray are co-authors**, with Eugene Katsevich, who
wrote sceptre. This is not an outside competitor — it is the neighbouring paper, and coordination is
possible rather than adversarial.

**What it does.** Derives an analytical formula for the power to detect perturbation–expression
associations, reproducing a simulation-based calculator's estimates while cutting runtime by up to
**seven orders of magnitude**, and wraps it in an interactive web app addressing 11 design questions
for Perturb-seq and TAP-seq.

**Formula inputs.** Sequenced reads per cell, number of cells, number of targeted elements, gRNAs per
target, perturbed cell count (`Nt = K × N × MOI / (KL + K0)`), plus baseline relative expression `α_j`
and dispersion `θ_j` from negative binomial fits to reference data.

**What it explicitly does *not* do** — quoting the paper:

> "While PerturbPlan is a tool for prospective power analysis and experimental design, it is sometimes
> of interest to analyze power retrospectively."

and, on using retrospective power for "curating perturbation-gene associations (and non-associations)
for predictive model training":

> "beyond the scope of the current work."

That second quote is *our* motivation, named and set aside. It is the strongest available evidence that
the gap is real, and it should be quoted in the introduction rather than paraphrased.

**Also not covered:** power is **aggregate**, not per pair — "Overall power is computed by averaging
the detection power over all element-gene pairs of interest". No confidence intervals on power. No
treatment of `trans` pairs or pairs that were never tested.

### Consequences for our framing

| Our claim | Status after PerturbPlan |
|---|---|
| "No satisfying method exists to quantify power" | **Dead as written.** Must become: no method gives *retrospective, per-pair* power with quantified uncertainty. |
| §4 — the curve is derived, one free parameter per pair, three points suffice | **Overlaps.** Their formula is analytical too, and from the same statistical structure. Our contribution is not "an analytical form exists" but that a *per-pair* curve can be fitted from few measured points and validated against a full sweep. |
| §6 — power predictable from perturbed cells and expression | **Substantially overlaps.** Their formula already takes perturbed-cell count, baseline expression and dispersion. Our `log(k) ~ log(pert_cells) + log(expression)` regression, R² 0.867, may simply be a weaker empirical version of a formula that already exists in closed form. |
| §6 — reaches pairs never simulated, including `trans` | **Still open.** They do not address `trans` or untested pairs. |
| §3 — per-pair power with a Wilson lower bound, minimum detectable effect size as the deliverable | **Still open.** They take a "minimum FC of interest" as an input, not an output, and report no intervals. |
| §7 — aggregation: per-pair vs per-element variance decomposition | **Still open**, and sharpened: they *average over pairs*, which is precisely the aggregation our §7 argues is not comparable across elements. |
| §8 — the inherited-null-model trap | **Untouched by anyone.** Genuinely new, and it bears on the simulation-based reference their formula is calibrated against. |

### The decision this forces

Our §6 fits a regression to predict `k`. PerturbPlan derives `k` analytically from the same quantities.
Fitting our own regression alongside an existing closed-form solution from our own lab is hard to
defend. Two better options:

1. **Apply their formula per pair and validate it retrospectively**, with calibrated intervals from
   measured residuals. This turns the overlap into a contribution: they showed the formula reproduces
   *aggregate* simulation power; nobody has shown it is trustworthy *per pair*, which is what a negative
   claim about a specific element–gene link requires. Extending it to `trans` follows naturally.
2. Keep the regression only as a diagnostic — where the analytical formula and the data disagree —
   rather than as the predictive engine.

Option 1 is stronger, cheaper, and cooperative. It also reframes the whole paper favourably: *the
analytical formula makes design tractable; we establish what it takes to trust it per pair, and use
that to triage negatives.* **[DECISION NEEDED]**

---

## Ray et al. — the simulation-based calculator PerturbPlan replaces

PerturbPlan names exactly one prior tool, its reference [9]: Ray J. et al., *An unbiased survey of
distal element-gene regulatory interactions with direct-capture targeted Perturb-seq.*

This is the DC-TAP paper, and it is the simulation-based power calculator described as "computationally
costly and requiring high-performance computing".

**Confirmed by EM: the pipeline in this repository is *not* that calculator.** They are separate
implementations. This matters in both directions and should be stated explicitly in §2:

- We cannot claim to be the tool PerturbPlan benchmarked against, so PerturbPlan's demonstration that
  its formula "recapitulates power estimates from the simulation-based tool" says nothing directly
  about agreement with *this* pipeline. That agreement is an open question and worth measuring — if our
  simulation and their formula disagree per pair, that is a result in itself.
- Conversely, the "computationally costly, requires HPC" criticism is levelled at a different tool, so
  our cost figures stand on their own and are not the ones PerturbPlan improved on by 10^7.

**[TODO]** get the full Ray et al. citation, and check what its supplement says about how power was
computed and at what granularity — it is the closest methodological relative even though it is not this
code.

**[TODO]** get the full citation, and check what its supplement says about how power was computed and at
what granularity.

---

## sceptre — the test being emulated

- Barry T, Wang X, Morris JA, Roeder K, Katsevich E. *SCEPTRE improves calibration and sensitivity in
  single-cell CRISPR screen analysis.* Genome Biology, 2021. High-MOI conditional resampling.
- Barry T, Mason K, Roeder K, Katsevich E. *Robust differential expression testing for single-cell
  CRISPR screens at low multiplicity of infection.* Genome Biology 25:124, 2024.
- Barry T, Roeder K, Katsevich E. *Exponential family measurement error models for single-cell CRISPR
  screens.* Biostatistics, 2024. Cite if the 'mixture' gRNA assignment is used.

**Relevant detail:** sceptre's own `run_power_check()` is *not* a power analysis in our sense — it
applies sceptre to **positive control pairs** to confirm the method can detect known associations. It
validates the method, not the experiment's power for a given pair. Worth one sentence in §2 to
forestall the "sceptre already does this" objection.

Note also that sceptre's test is a conditional-resampling test with a skew-normal approximation, not a
Wald test. Our §4 curve assumes a Wald-type form, which is why fitted `z` (median 2.045 on `wtc11`)
does not exactly match `qnorm(1 - threshold)`. PerturbPlan presumably confronts the same issue —
**[TODO]** check how they handle it, since it is the main theoretical soft spot in §4.

---

## scPower — adjacent, and the comparison reviewers will still ask for

Schmid KT, Höllbacher B, Cruceanu C, Böttcher A, Lickert H, Binder EB, Theis FJ, Heinig M.
*scPower accelerates and optimizes the design of multi-sample single cell transcriptomic studies.*
Nature Communications 12:6625, 2021.
[Nature](https://www.nature.com/articles/s41467-021-26779-7) ·
[GitHub](https://github.com/heiniglab/scPower)

Design-stage power for multi-sample scRNA-seq **DE and eQTL** studies: negative binomial regression on
pseudo-bulk counts, plus a model for the probability of detecting cell-type-specific expression as a
function of cells per type and reads per cell. Optimises sample size / cells / depth under a budget.

**Why it is not our competitor:** different inferential target (inter-individual DE/eQTL, not
perturbation–response association), prospective, and no notion of a per-pair negative. One sentence in
§2 is enough — but it must be there, because it is the best-known power tool in single-cell.

---

## Still to check

- **Ray et al. supplement** — how power was computed, at what granularity. Highest priority, since it
  is the direct ancestor.
- Gasperini M, et al. *A Genome-wide Framework for Mapping Gene Regulation via Cellular Genetic
  Screens.* Cell, 2019 — what power analysis, if any, it reports. Needed anyway as a validation dataset.
- scE2G and ENCODE-rE2G papers — how they actually treat negatives, and whether any power filter is
  already applied. The motivation in §1 depends on this. If they already filter on power, the
  contribution shifts from "they should" to "here is how to do it properly".
- Whether PerturbPlan's formula has been applied per pair anywhere since May 2026.
- CRISPR screen power literature outside single-cell readouts (bulk FACS-based screens) — likely not
  transferable, but worth one search to be able to say so.

## Searches run

- `scPower power analysis single-cell RNA-seq design Schmid` (web)
- `power analysis single-cell CRISPR screen Perturb-seq sample size tool` (web) → found PerturbPlan
- `sceptre power analysis CRISPR screen Barry Katsevich` (web)
- `statistical power CRISPRi single-cell screen enhancer gene pair` (PubMed) → 0 results
- PubMed metadata for PMID 42239142; PMC13228452 full text

bioRxiv blocks direct fetching (403); PMC is the working route.
