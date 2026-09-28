"""Drawing the counts.

One function, and the only place randomness touches expression. Everything it
needs has been decided elsewhere: the baseline by `baseline.py`, the per-cell
multiplier by `perturbation.py`, the dispersion by `gene_model.py`.

    mu[i, j] = baseline[i, j] * effect_size[i, j]
    count[i, j] ~ NegBinomial(mean = mu[i, j], size = theta[i])

`size` is theta, the same parameterisation sceptre's own model uses, so the
counts this draws are counts from the model the test assumes -- which is the
point of taking the baseline from that model too (see `baseline.py`).
"""

from __future__ import annotations

import numpy as np

# int16 holds counts to 32,767, which covers almost every draw. Not all of them: sceptre's fitted
# mean for a very highly expressed gene can reach far past anything observed in the cells with
# extreme covariates -- HBA2 in DC-TAP K562 has an observed maximum of 6,583 and a fitted mean of
# 20,104 in its most extreme cell, and a draw from NB(mean 20,104, theta 2.1) passes 32,767 about
# one time in seven. Those draws are legitimate draws from the model the test assumes, so a draw
# that does not fit is promoted to the next integer width rather than refused.
_COUNT_DTYPE = np.int16
_WIDER = (np.int16, np.int32, np.int64)


def draw_counts(
    baseline: np.ndarray,
    effect_size: np.ndarray,
    theta: np.ndarray,
    rng: np.random.Generator,
    *,
    dtype: np.dtype | None = _COUNT_DTYPE,
) -> np.ndarray:
    """Simulated counts, `(n_genes, n_cells)`.

    `theta` is the NB size, one per gene, broadcast down the rows.

    numpy parameterises the negative binomial as `(n, p)` with mean
    `n(1-p)/p`, so `n = theta` and `p = theta / (theta + mu)` give mean `mu`
    and variance `mu + mu^2/theta` -- R's `rnbinom(mu=, size=)`.

    Returned as `int16` by default, or the narrowest of int16/int32/int64 at
    least as wide as `dtype` that holds the largest draw. These are counts;
    holding them as float64 costs four times the memory for no information, and
    the simulation's whole shape depends on how many replicates fit in one
    process. A draw is never wrapped or clipped.
    """
    baseline = np.asarray(baseline, dtype=float)
    if baseline.shape != effect_size.shape:
        raise ValueError(f"baseline is {baseline.shape} but effect_size is {effect_size.shape}")
    theta = np.asarray(theta, dtype=float)
    if theta.shape != (baseline.shape[0],):
        raise ValueError(f"theta is {theta.shape}, expected ({baseline.shape[0]},)")
    if np.any(theta <= 0) or not np.all(np.isfinite(theta)):
        raise ValueError("theta must be finite and positive")

    mu = baseline * effect_size
    size = theta[:, None]
    # p == 1 exactly where mu == 0, which numpy accepts and which draws 0
    # every time -- the right answer for a gene knocked all the way down.
    p = size / (size + mu)
    counts = rng.negative_binomial(size, p)

    if dtype is None:
        return counts
    largest = counts.max(initial=0)
    for candidate in _WIDER:
        if np.dtype(candidate).itemsize >= np.dtype(dtype).itemsize and largest <= np.iinfo(
            candidate
        ).max:
            return counts.astype(candidate)
    return counts
