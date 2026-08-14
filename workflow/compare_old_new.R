#!/usr/bin/env Rscript
#
# Compare the refactored pipeline's power estimates against the old pipeline's.
#
# This lives in workflow/ rather than src/ deliberately: src/ is the pipeline surface that the
# workflow engine stages onto PATH, and a one-off validation tool does not belong there.
#
# WHAT IS AND IS NOT COMPARABLE
#
# The comparison is distributional, not exact, and cannot be made exact. Seeding moved inside the
# replicate loop during the refactor, which changes the RNG stream, and removing the size-factor
# shuffle removed a sample() call from it as well. Two runs of the *old* pipeline do not agree
# with each other either -- it has no set.seed() anywhere. So the yardstick for "do these agree"
# is the old pipeline's own run-to-run spread, and the Monte-Carlo error of the replicate count.
#
# What *is* directly comparable is anything that does not depend on the RNG: the set of pairs
# tested, and the number of perturbed cells per pair. Those are checked exactly, because a
# difference there is a wiring bug rather than noise.
#
# Usage:
#   compare_old_new.R --old  results/day0_grna20_no_shuffle/power_analysis/power_analysis_results_es_0.15.tsv \
#                     --new  results/refactor/power/power_es0.15.tsv \
#                     --outdir results/refactor/comparison

## ARGUMENTS =======================================================================================

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) {
    if (is.null(default)) {
      stop("Missing required argument ", flag, call. = FALSE)
    }
    return(default)
  }
  if (idx == length(args)) stop(flag, " needs a value.", call. = FALSE)
  args[idx + 1]
}

old_path <- get_arg("--old")
new_path <- get_arg("--new")
outdir <- get_arg("--outdir", ".")
# Correlation below this is treated as a failure. From docs/status.md: > 0.95 at 100 replicates.
min_correlation <- as.numeric(get_arg("--min-correlation", "0.95"))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

report_path <- file.path(outdir, "comparison_report.txt")
pairs_path <- file.path(outdir, "comparison_pairs.tsv")

# Everything printed also lands in the report file, so the run is self-documenting.
report_lines <- character(0)
say <- function(...) {
  line <- paste0(...)
  report_lines <<- c(report_lines, line)
  cat(line, "\n", sep = "")
}

## LOAD ============================================================================================

read_tsv <- function(path) {
  if (!file.exists(path)) stop("No such file: ", path, call. = FALSE)
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

old <- read_tsv(old_path)
new <- read_tsv(new_path)

for (required in c("grna_target", "response_id", "power")) {
  if (!required %in% colnames(old)) stop("Old table lacks column ", required, call. = FALSE)
  if (!required %in% colnames(new)) stop("New table lacks column ", required, call. = FALSE)
}

say("OLD VS NEW POWER COMPARISON")
say("===========================")
say("")
say("old:  ", old_path)
say("      ", nrow(old), " pairs")
say("new:  ", new_path)
say("      ", nrow(new), " pairs")
say("")

failures <- character(0)
fail <- function(...) failures <<- c(failures, paste0(...))

## CHECK 1 -- PAIR SETS ============================================================================
#
# Exact check. The two pipelines derive the pair list differently (the old one from a separate
# gene_grna_group_pairs.rds, the new one from @discovery_pairs_with_info inside the object), so
# this is a genuine test that they agree on what is being measured, not a formality.

say("CHECK 1 -- pair sets")
say("--------------------")

old_key <- paste(old$grna_target, old$response_id, sep = "\r")
new_key <- paste(new$grna_target, new$response_id, sep = "\r")

if (anyDuplicated(old_key)) fail("The old table has ", sum(duplicated(old_key)), " duplicated pairs.")
if (anyDuplicated(new_key)) fail("The new table has ", sum(duplicated(new_key)), " duplicated pairs.")

only_old <- setdiff(old_key, new_key)
only_new <- setdiff(new_key, old_key)
shared <- intersect(old_key, new_key)

say("  shared pairs:      ", length(shared))
say("  only in old:       ", length(only_old))
say("  only in new:       ", length(only_new))

show_examples <- function(keys, label) {
  if (length(keys) == 0) return(invisible(NULL))
  parts <- do.call(rbind, strsplit(utils::head(keys, 5), "\r", fixed = TRUE))
  say("  first ", min(5, length(keys)), " ", label, ":")
  for (i in seq_len(nrow(parts))) say("    ", parts[i, 1], "  ", parts[i, 2])
}
show_examples(only_old, "only in old")
show_examples(only_new, "only in new")

if (length(only_old) > 0 || length(only_new) > 0) {
  fail("Pair sets differ: ", length(only_old), " only in old, ", length(only_new), " only in new.")
} else {
  say("  PASS -- identical pair sets")
}
say("")

if (length(shared) == 0) {
  say("No shared pairs; nothing further to compare.")
  writeLines(report_lines, report_path)
  stop("The two tables share no pairs.", call. = FALSE)
}

## JOIN ============================================================================================

merged <- merge(
  old[, intersect(c("grna_target", "response_id", "power", "mean_log_2_fold_change",
                    "mean_pert_cells"), colnames(old))],
  new[, intersect(c("grna_target", "response_id", "power", "mean_log_2_fold_change",
                    "mean_pert_cells", "power_ci_low", "power_ci_high", "n_reps"),
                  colnames(new))],
  by = c("grna_target", "response_id"),
  suffixes = c("_old", "_new")
)
merged$power_diff <- merged$power_new - merged$power_old

## CHECK 2 -- PERTURBED CELL COUNTS ================================================================
#
# Does not depend on the RNG: the same cells carry the same guides in both pipelines. Any
# disagreement is a wiring difference in how perturbation status is derived, which would make the
# power comparison meaningless -- so it is worth checking before looking at power at all.

say("CHECK 2 -- perturbed cells per pair (RNG-independent)")
say("----------------------------------------------------")

if (all(c("mean_pert_cells_old", "mean_pert_cells_new") %in% colnames(merged))) {
  cell_diff <- merged$mean_pert_cells_new - merged$mean_pert_cells_old
  n_differ <- sum(abs(cell_diff) > 1e-6, na.rm = TRUE)
  say("  pairs differing:   ", n_differ, " of ", nrow(merged))
  say("  max |difference|:  ", format(max(abs(cell_diff), na.rm = TRUE), digits = 6))
  if (n_differ > 0) {
    say("  mean difference:   ", format(mean(cell_diff, na.rm = TRUE), digits = 6))
    worst <- merged[order(-abs(cell_diff)), ][seq_len(min(5, nrow(merged))), ]
    say("  largest differences:")
    for (i in seq_len(nrow(worst))) {
      say("    ", worst$grna_target[i], "  ", worst$response_id[i],
          "  old ", worst$mean_pert_cells_old[i], "  new ", worst$mean_pert_cells_new[i])
    }
    fail(n_differ, " pair(s) have different perturbed-cell counts. This does not depend on the ",
         "RNG, so it indicates a real difference in how perturbation status is derived.")
  } else {
    say("  PASS -- identical")
  }
} else {
  say("  SKIPPED -- mean_pert_cells not present in both tables")
}
say("")

## CHECK 3 -- MEAN POWER ===========================================================================

say("CHECK 3 -- mean power")
say("---------------------")

mean_old <- mean(merged$power_old, na.rm = TRUE)
mean_new <- mean(merged$power_new, na.rm = TRUE)
mean_diff <- mean_new - mean_old

# Monte-Carlo standard error of the *mean over pairs*. Each pair's power is a binomial proportion
# over n_reps replicates, so its standard error is sqrt(p(1-p)/n_reps); averaging over N pairs
# divides the variance by N. Both runs carry this error, hence the sqrt(2).
n_reps <- if ("n_reps" %in% colnames(merged)) stats::median(merged$n_reps, na.rm = TRUE) else 100
per_pair_var <- merged$power_new * (1 - merged$power_new) / n_reps
se_mean <- sqrt(2 * sum(per_pair_var, na.rm = TRUE)) / nrow(merged)

say("  old:               ", format(mean_old, digits = 5))
say("  new:               ", format(mean_new, digits = 5))
say("  difference:        ", format(mean_diff, digits = 5))
say("  Monte-Carlo SE:    ", format(se_mean, digits = 4), "  (at ", n_reps, " replicates)")
say("  difference in SE:  ", format(mean_diff / se_mean, digits = 4), " sigma")

# Deliberately loose: the point is to catch a systematic shift, not to fail on a real but tiny
# difference. Anything past 5 sigma on ~35,000 pairs is a shift, not noise.
if (abs(mean_diff) > 5 * se_mean) {
  fail("Mean power differs by ", format(mean_diff, digits = 4), " = ",
       format(abs(mean_diff / se_mean), digits = 3), " Monte-Carlo standard errors. ",
       "This is a systematic shift, not sampling noise.")
} else {
  say("  PASS -- within 5 Monte-Carlo standard errors")
}
say("")

## CHECK 4 -- PAIRED SIGN TEST =====================================================================
#
# Asks whether the new pipeline is systematically higher or lower, ignoring magnitude. Ties
# (identical power, common when both are 0 or 1) carry no directional information and are dropped,
# which is the standard treatment.

say("CHECK 4 -- paired sign test on per-pair differences")
say("---------------------------------------------------")

n_up <- sum(merged$power_diff > 0, na.rm = TRUE)
n_down <- sum(merged$power_diff < 0, na.rm = TRUE)
n_tied <- sum(merged$power_diff == 0, na.rm = TRUE)

say("  new > old:         ", n_up)
say("  new < old:         ", n_down)
say("  tied:              ", n_tied)

if (n_up + n_down > 0) {
  sign_p <- stats::binom.test(n_up, n_up + n_down, p = 0.5)$p.value
  say("  sign test p:       ", format.pval(sign_p, digits = 4))
  say("  median difference: ", format(stats::median(merged$power_diff, na.rm = TRUE), digits = 4))
  say("  mean difference:   ", format(mean(merged$power_diff, na.rm = TRUE), digits = 4))

  # With ~35,000 pairs the sign test detects offsets far too small to matter, so a bare p-value
  # is not a useful gate on its own. It is reported, and paired with the effect size: a
  # significant p AND a median offset past 0.02 is what counts as a real discrepancy.
  median_diff <- stats::median(merged$power_diff, na.rm = TRUE)
  if (sign_p < 0.05 && abs(median_diff) > 0.02) {
    fail("Directional bias: p = ", format.pval(sign_p, digits = 3), " with a median offset of ",
         format(median_diff, digits = 3), ".")
  } else if (sign_p < 0.05) {
    say("  NOTE -- statistically significant but the median offset is under 0.02, which at this")
    say("          number of pairs is expected. Not treated as a failure.")
  } else {
    say("  PASS -- no directional bias")
  }
} else {
  say("  All differences are exactly zero.")
}
say("")

## CHECK 5 -- CORRELATION ==========================================================================

say("CHECK 5 -- per-pair agreement")
say("----------------------------")

usable <- stats::complete.cases(merged$power_old, merged$power_new)
correlation <- if (stats::sd(merged$power_old[usable]) == 0 ||
                   stats::sd(merged$power_new[usable]) == 0) {
  NA_real_
} else {
  stats::cor(merged$power_old[usable], merged$power_new[usable])
}

abs_diff <- abs(merged$power_diff)
say("  Pearson r:         ", format(correlation, digits = 5))
say("  mean |difference|: ", format(mean(abs_diff, na.rm = TRUE), digits = 4))
say("  median |diff|:     ", format(stats::median(abs_diff, na.rm = TRUE), digits = 4))
say("  90th percentile:   ", format(stats::quantile(abs_diff, 0.9, na.rm = TRUE), digits = 4))
say("  max |difference|:  ", format(max(abs_diff, na.rm = TRUE), digits = 4))
say("  pairs differing by more than 0.1:  ",
    sum(abs_diff > 0.1, na.rm = TRUE), " of ", nrow(merged))
say("  pairs differing by more than 0.2:  ",
    sum(abs_diff > 0.2, na.rm = TRUE), " of ", nrow(merged))

if (is.na(correlation)) {
  say("  NOTE -- correlation undefined (one side has zero variance)")
} else if (correlation < min_correlation) {
  fail("Per-pair correlation is ", format(correlation, digits = 4), ", below the ",
       min_correlation, " expected at ", n_reps, " replicates.")
} else {
  say("  PASS -- correlation at or above ", min_correlation)
}
say("")

# Where the two disagree most, for follow-up. Scatter of +/-0.1 on individual pairs is expected
# at 100 replicates; a cluster of large differences sharing a target is not.
worst <- merged[order(-abs_diff), ][seq_len(min(10, nrow(merged))), ]
say("  10 largest per-pair differences:")
say(sprintf("    %-34s %-14s %8s %8s %8s", "grna_target", "response_id", "old", "new", "diff"))
for (i in seq_len(nrow(worst))) {
  say(sprintf("    %-34s %-14s %8.3f %8.3f %8.3f",
              worst$grna_target[i], worst$response_id[i],
              worst$power_old[i], worst$power_new[i], worst$power_diff[i]))
}
say("")

## CHECK 6 -- CONFIDENCE INTERVAL COVERAGE =========================================================
#
# The new pipeline reports a Wilson interval per pair. If the two pipelines agree, the old
# estimate should fall inside it about conf_level of the time. Systematically low coverage means
# the difference is larger than replicate noise explains.

if (all(c("power_ci_low", "power_ci_high") %in% colnames(merged))) {
  say("CHECK 6 -- old estimate within the new Wilson interval")
  say("------------------------------------------------------")
  inside <- merged$power_old >= merged$power_ci_low & merged$power_old <= merged$power_ci_high
  coverage <- mean(inside, na.rm = TRUE)
  say("  coverage:          ", format(100 * coverage, digits = 4), "%")
  say("  (both estimates carry replicate noise, so this understates true agreement;")
  say("   informational, not a pass/fail gate)")
  say("")
}

## WRITE ===========================================================================================

utils::write.table(merged[order(-abs(merged$power_diff)), ], pairs_path,
                   sep = "\t", quote = FALSE, row.names = FALSE)

say("VERDICT")
say("=======")
if (length(failures) == 0) {
  say("  PASS -- the refactored pipeline agrees with the old one on every check.")
} else {
  say("  FAIL -- ", length(failures), " check(s) did not pass:")
  for (f in failures) say("    - ", f)
}
say("")
say("per-pair detail: ", pairs_path)

writeLines(report_lines, report_path)
cat("\nreport written to", report_path, "\n")

if (length(failures) > 0) quit(status = 1)
