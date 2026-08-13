# Paper

Skeleton for the methods paper on quantifying power in single-cell CRISPRi element–gene screens.

**Nothing here is written prose yet.** `outline.md` is the flow with placeholders; `figures.md` is the
figure list with, for each panel, the data and script that would produce it and what is still missing.

## What the paper claims

In one sentence, as currently understood:

> Power for element–gene pairs in single-cell CRISPRi screens can be quantified per pair, cheaply
> enough to run at screen scale, and **predicted for pairs that were never simulated** — including
> `trans` pairs, where direct simulation is infeasible.

The three supporting claims, in the order they appear:

1. Distinguishing an underpowered pair from a real negative is required by anything trained on CRISPR
   negatives (scE2G, ENCODE-rE2G) and by any analysis that conditions on a screen's null results.
2. A full effect-size sweep is not necessary: the power curve is derived, has one free parameter per
   pair, and three well-placed effect sizes determine it.
3. The curve's parameter is itself predictable from perturbed-cell count and gene expression, which
   extends power estimates to unmeasured pairs.

## Conventions

- Markdown, one file per section once drafting starts. `outline.md` stays as the map.
- Figure placeholders are written as `<!-- FIGURE n: short description -->` so they can be grepped.
- Every quoted number carries its provenance inline — dataset, replicate count, and the job or script
  that produced it. Numbers without provenance are marked `[unverified]`.
- Claims not yet supported by data are marked `**[NEEDS DATA]**` and cross-referenced to
  `docs/status.md`, which holds the measurements and the open gaps.

## Status

| Piece | State |
|---|---|
| Flow / outline | drafted here, from EM's framing |
| Figure list | drafted, mostly placeholders |
| Literature check | **not started** — see `outline.md` §2 |
| Curve validation, 2nd dataset | **not run** — only `dc_tap_paper_wtc11_no_shuf` has a six-point sweep |
| Covariate model, cross-dataset | **not run** |
| `trans` demonstration | **not run** |

The measurements, cost models and known traps live in
[`docs/status.md`](../docs/status.md) — this directory should cite it rather than restate it.
