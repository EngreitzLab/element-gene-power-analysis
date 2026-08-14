source_cli()

test_that("derive_seed is deterministic", {
  expect_equal(
    derive_seed(20250812, "chr1:100-200", 7, 0.15),
    derive_seed(20250812, "chr1:100-200", 7, 0.15)
  )
})

test_that("derive_seed separates every component of the key", {
  base <- derive_seed(20250812, "chr1:100-200", 7, 0.15)
  expect_false(base == derive_seed(20250813, "chr1:100-200", 7, 0.15))  # base seed
  expect_false(base == derive_seed(20250812, "chr1:100-300", 7, 0.15))  # target
  expect_false(base == derive_seed(20250812, "chr1:100-200", 8, 0.15))  # replicate
  expect_false(base == derive_seed(20250812, "chr1:100-200", 7, 0.20))  # effect size
})

test_that("derive_seed does not depend on how work is partitioned", {
  # This is the property the whole seeding scheme exists for. Replicates 1-4 must get the same
  # seeds whether they are run as one chunk of 4 or two chunks of 2, or --n-splits and
  # --reps-per-chunk -- purely computational knobs -- would change the reported power.
  one_chunk <- vapply(1:4, function(r) derive_seed(20250812, "chr1:100-200", r, 0.15), integer(1))
  two_chunks <- c(
    vapply(1:2, function(r) derive_seed(20250812, "chr1:100-200", r, 0.15), integer(1)),
    vapply(3:4, function(r) derive_seed(20250812, "chr1:100-200", r, 0.15), integer(1))
  )
  expect_equal(one_chunk, two_chunks)
})

test_that("derive_seed returns a valid, in-range integer seed", {
  seeds <- vapply(
    c("chr1:1-2", "chr19:10936677-10936987", "non-targeting", NULL_FIT_TARGET_KEY),
    function(t) derive_seed(20250812, t, 1, 0.15),
    integer(1)
  )
  expect_type(seeds, "integer")
  expect_true(all(is.finite(seeds)))
  expect_true(all(seeds >= 0))
})

test_that("the null-fit key cannot collide with a real target", {
  # Targets are genomic coordinates or control labels; the sentinel is deliberately unusable as one.
  expect_true(grepl("^__", NULL_FIT_TARGET_KEY))
  expect_false(grepl("^__", "chr19:10936677-10936987"))
})

test_that("stable_hash spreads distinct keys across distinct seeds", {
  # Not a claim about cryptographic quality -- only that ~3,000 targets x 100 replicates do not
  # pile onto a handful of seeds, which would correlate simulations that must be independent.
  keys <- as.vector(outer(
    sprintf("chr1:%d-%d", seq_len(200) * 1000, seq_len(200) * 1000 + 300),
    1:20,
    function(t, r) paste(20250812, t, r, 0.15, sep = "|")
  ))
  hashes <- vapply(keys, stable_hash, integer(1))
  expect_gt(length(unique(hashes)) / length(hashes), 0.999)
})

test_that("init_seed accepts what derive_seed produces, and says which seed it used", {
  # The message is deliberate, not noise: without it a run's log gives no way to reproduce it.
  seed <- derive_seed(20250812, "chr1:100-200", 1, 0.15)
  expect_message(init_seed(seed), "RNG seed")
  expect_equal(suppressMessages(init_seed(seed)), seed)
})

test_that("init_seed refuses to run unseeded", {
  # Before this, nothing called set.seed() anywhere and no run was reproducible.
  expect_error(init_seed(NULL), "A seed is required")
  expect_error(init_seed(NA), "A seed is required")
})

test_that("init_seed actually sets the stream", {
  init_seed(42)
  a <- suppressMessages({ init_seed(42); runif(3) })
  b <- suppressMessages({ init_seed(42); runif(3) })
  expect_equal(a, b)
})
