#!/usr/bin/env Rscript
#
# Compare two null-model bundles on the replicates they share.
#
# Used to check the Nextflow runner against the sbatch runner. The two need not cover the same
# number of replicates -- a short test run covers 2 where a production run covers 100 -- so the
# comparison is over the intersection of their replicate keys.
#
# The replicates SHOULD agree exactly. A null model is fitted on a null simulation seeded from
# (seed, rep) alone, deliberately not from (seed, target, rep, effect_size), so replicate 1 is the
# same object no matter which runner produced it, how many replicates it produced, or how they were
# chunked across tasks. Anything other than an exact match means the two runners are not doing the
# same arithmetic, and the whole point of keeping both is that they should be.
#
# Usage:
#   compare_null_models.R --a bundle_a.rds --b bundle_b.rds

local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here <- if (length(file_arg) == 1) {
    dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else {
    normalizePath(".")
  }
  for (dir in unique(c(file.path(dirname(here), "lib"), file.path(here, "lib"), here))) {
    if (file.exists(file.path(dir, "cli.R"))) {
      source(file.path(dir, "cli.R"))
      return(invisible(NULL))
    }
  }
  stop("Cannot find lib/cli.R relative to ", here, call. = FALSE)
})

option_list <- list(
  make_option("--a", type = "character", default = NULL, dest = "a",
              help = "First null-model bundle (.rds)."),
  make_option("--b", type = "character", default = NULL, dest = "b",
              help = "Second null-model bundle (.rds).")
)
opts <- parse_args(OptionParser(option_list = option_list))
require_options(opts, c("a", "b"))

say <- function(...) cat(..., "\n", sep = "")

a <- readRDS(opts$a)
b <- readRDS(opts$b)

say("A: ", opts$a)
say("   seed ", a$seed, ", ", length(a$precomputations), " replicate(s)")
say("B: ", opts$b)
say("   seed ", b$seed, ", ", length(b$precomputations), " replicate(s)")
say("")

if (!identical(a$seed, b$seed)) {
  stop("Seeds differ (", a$seed, " vs ", b$seed, "); the bundles are not comparable.",
       call. = FALSE)
}

shared <- intersect(names(a$precomputations), names(b$precomputations))
if (length(shared) == 0) {
  stop("The bundles share no replicate.", call. = FALSE)
}
say("shared replicates: ", length(shared), " (", paste(head(shared, 10), collapse = ", "),
    if (length(shared) > 10) ", ..." else "", ")")
say("")

failures <- 0L
for (rep_id in shared) {
  ra <- a$precomputations[[rep_id]]
  rb <- b$precomputations[[rep_id]]

  genes_a <- names(ra)
  genes_b <- names(rb)
  if (!setequal(genes_a, genes_b)) {
    say("  rep ", rep_id, ": FAIL -- gene sets differ (",
        length(setdiff(genes_a, genes_b)), " only in A, ",
        length(setdiff(genes_b, genes_a)), " only in B)")
    failures <- failures + 1L
    next
  }

  if (identical(ra[sort(genes_a)], rb[sort(genes_a)])) {
    say("  rep ", rep_id, ": identical (", length(genes_a), " genes)")
  } else {
    # Not identical -- quantify, so the report says how far off rather than just "no".
    worst <- 0
    for (g in genes_a) {
      d <- tryCatch(max(abs(unlist(ra[[g]]) - unlist(rb[[g]]))), error = function(e) NA_real_)
      if (!is.na(d) && d > worst) worst <- d
    }
    say("  rep ", rep_id, ": FAIL -- max abs coefficient difference ", format(worst, digits = 6))
    failures <- failures + 1L
  }
}

say("")
if (failures > 0L) {
  say("VERDICT: FAIL -- ", failures, " of ", length(shared), " shared replicate(s) differ.")
  quit(status = 1L)
}
say("VERDICT: PASS -- all ", length(shared), " shared replicate(s) identical.")
