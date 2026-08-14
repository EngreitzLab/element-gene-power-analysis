test_that("wilson_interval matches published values at the boundaries", {
  # The boundaries are the whole reason this is not the normal approximation, so they are what the
  # test pins. Reference values are the standard Wilson score interval at 95 %.
  zero <- wilson_interval(0, 100)
  expect_equal(zero$low, 0)
  expect_equal(zero$high, 0.0370, tolerance = 1e-3)

  full <- wilson_interval(100, 100)
  expect_equal(full$low, 0.9630, tolerance = 1e-3)
  expect_equal(full$high, 1)

  half <- wilson_interval(50, 100)
  expect_equal(half$low, 0.4038, tolerance = 1e-3)
  expect_equal(half$high, 0.5962, tolerance = 1e-3)
})

test_that("wilson_interval does not collapse where the normal approximation does", {
  # p +/- z*sqrt(p(1-p)/n) gives [0, 0] for 0 successes, which would assert that power is certainly
  # zero and let a pair be called underpowered with no uncertainty at all.
  expect_gt(wilson_interval(0, 100)$high, 0)
  expect_lt(wilson_interval(100, 100)$low, 1)
})

test_that("wilson_interval is vectorised and stays within [0, 1]", {
  res <- wilson_interval(c(0, 1, 50, 99, 100), 100)
  expect_length(res$low, 5)
  expect_true(all(res$low >= 0 & res$low <= 1))
  expect_true(all(res$high >= 0 & res$high <= 1))
  expect_true(all(res$low <= res$high))
})

test_that("wilson_interval widens as the confidence level rises", {
  narrow <- wilson_interval(50, 100, conf_level = 0.80)
  wide   <- wilson_interval(50, 100, conf_level = 0.99)
  expect_lt(wide$low, narrow$low)
  expect_gt(wide$high, narrow$high)
})

test_that("wilson_interval narrows as replicates increase", {
  expect_lt(
    diff(unlist(wilson_interval(500, 1000)[c("low", "high")])),
    diff(unlist(wilson_interval(50, 100)[c("low", "high")]))
  )
})

test_that("effect_label keeps whole percentages whole", {
  # The bug this pins: trimming trailing zeros unconditionally turned 0.2 into "2" and 0.5 into
  # "5", so a 20 % knockdown was reported in a column named power_at_effect_size_2.
  expect_equal(effect_label(0.2), "20")
  expect_equal(effect_label(0.5), "50")
  expect_equal(effect_label(0.15), "15")
  expect_equal(effect_label(0.05), "5")
  expect_equal(effect_label(0.1), "10")
})

test_that("effect_label keeps a real fractional part", {
  expect_equal(effect_label(0.125), "12_5")
  expect_equal(effect_label(0.025), "2_5")
})

test_that("effect_label emits no dots, which downstream readers treat as separators", {
  for (es in c(0.05, 0.1, 0.125, 0.15, 0.2, 0.25, 0.5)) {
    expect_false(grepl(".", effect_label(es), fixed = TRUE), info = paste("effect size", es))
  }
})

test_that("effect_label is distinct across a realistic sweep", {
  sweep <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.5)
  expect_equal(length(unique(vapply(sweep, effect_label, character(1)))), length(sweep))
})
