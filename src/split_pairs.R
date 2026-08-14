#!/usr/bin/env Rscript
#
# Split the QC-passing discovery pairs into balanced chunks, one per parallel task.
#
# Whole targets stay together, because the simulation works target by target: it builds the
# perturbation status and guide-level effect sizes once per target and then loops over reps.
# Splitting a target across tasks would duplicate that setup.
#
# Two changes from the original split_target_response_pairs.R:
#
#  * Balance on pairs, not on target count. The original counted unique targets per bin
#    (`sapply(splits, n_distinct)`), so a target with 36 pairs weighed the same as one with 1.
#    Measured on sample1 this was worth ~18% on the straggler bin -- real but modest.
#  * O(n log n) instead of O(n^2). The original recomputed n_distinct() across every growing bin
#    for each of the 2,798 targets in turn.
#
# It also writes a header row. The original wrote col_names = FALSE and the reader re-labelled
# the columns positionally as c("grna_group", "response_id") while the writer had written
# grna_target -- correct only by accident of column order.
#
# Usage:
#   split_pairs.R --pairs pairs.tsv --n-splits 10 --outdir splits/

local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  lib <- if (length(file_arg) == 1) {
    file.path(dirname(dirname(normalizePath(sub("^--file=", "", file_arg)))), "lib")
  } else {
    "lib"
  }
  source(file.path(lib, "cli.R"))
})

## ARGUMENTS =======================================================================================

option_list <- list(
  make_option("--pairs", type = "character", default = NULL, dest = "pairs",
              help = "TSV of discovery pairs with columns grna_target and response_id."),
  make_option("--n-splits", type = "integer", default = NULL, dest = "n_splits",
              help = "Number of splits to produce."),
  make_option("--outdir", type = "character", default = ".", dest = "outdir",
              help = "Directory to write the splits into [default %default]."),
  make_option("--prefix", type = "character", default = "split_", dest = "prefix",
              help = "Filename prefix for each split [default %default]."),
  make_option("--target-overhead", type = "double", default = 0, dest = "target_overhead",
              help = paste("Extra weight per target, in pair-equivalents, to account for the",
                           "fixed per-target cost of a sceptre call. 0 balances purely on pair",
                           "count [default %default]. Raise it once the per-target overhead has",
                           "been measured; see docs/development.md."))
)

opts <- parse_args(OptionParser(
  option_list = option_list,
  description = "Split discovery pairs into balanced per-task chunks."
))
require_options(opts, c("pairs", "n_splits"))

if (opts$n_splits < 1) {
  stop("--n-splits must be at least 1.", call. = FALSE)
}
if (opts$target_overhead < 0) {
  stop("--target-overhead must not be negative.", call. = FALSE)
}

## SPLIT ===========================================================================================

pairs <- read_tsv_file(opts$pairs, required_columns = c("grna_target", "response_id"))
log_step("Read ", nrow(pairs), " pairs from ", opts$pairs)

if (nrow(pairs) == 0) {
  stop(opts$pairs, " contains no pairs.", call. = FALSE)
}

pairs_per_target <- table(pairs$grna_target)
n_target <- length(pairs_per_target)
log_step("Targets: ", n_target, " (", min(pairs_per_target), "-", max(pairs_per_target),
         " pairs each, median ", stats::median(pairs_per_target), ")")

if (opts$n_splits > n_target) {
  stop("--n-splits (", opts$n_splits, ") exceeds the number of targets (", n_target, "). ",
       "Every split must contain at least one target; lower --n-splits.", call. = FALSE)
}

# Longest-processing-time-first bin packing: sort targets by weight descending, then repeatedly
# place the next target into whichever split is currently lightest. Simple, deterministic, and
# within a small constant factor of optimal for this shape of problem.
weights <- as.numeric(pairs_per_target) + opts$target_overhead
names(weights) <- names(pairs_per_target)
weights <- sort(weights, decreasing = TRUE)

split_load <- numeric(opts$n_splits)
split_of_target <- integer(length(weights))
for (i in seq_along(weights)) {
  lightest <- which.min(split_load)
  split_of_target[i] <- lightest
  split_load[lightest] <- split_load[lightest] + weights[[i]]
}
names(split_of_target) <- names(weights)

## WRITE ===========================================================================================

if (!dir.exists(opts$outdir)) {
  dir.create(opts$outdir, recursive = TRUE)
}

# Zero-padded so lexicographic order matches numeric order, which keeps Nextflow channel
# ordering and any manual `ls` predictable.
width <- max(2, nchar(as.character(opts$n_splits)))
pairs$split <- split_of_target[pairs$grna_target]

written <- character(opts$n_splits)
for (k in seq_len(opts$n_splits)) {
  chunk <- pairs[pairs$split == k, c("grna_target", "response_id"), drop = FALSE]
  path <- file.path(opts$outdir,
                    sprintf("%s%0*d.tsv", opts$prefix, width, k))
  write_tsv_file(chunk[order(chunk$grna_target, chunk$response_id), ], path)
  written[k] <- path
}

pairs_written <- sum(vapply(written, function(p) nrow(read_tsv_file(p)), numeric(1)))
if (pairs_written != nrow(pairs)) {
  stop("Wrote ", pairs_written, " pairs but read ", nrow(pairs),
       "; the split lost or duplicated rows.", call. = FALSE)
}

targets_per_split <- tabulate(split_of_target, nbins = opts$n_splits)
pairs_per_split <- vapply(seq_len(opts$n_splits),
                          function(k) sum(pairs$split == k), numeric(1))

log_step("Wrote ", opts$n_splits, " splits to ", opts$outdir)
message("  split  targets   pairs")
for (k in seq_len(opts$n_splits)) {
  message(sprintf("  %5d  %7d  %6d", k, targets_per_split[k], pairs_per_split[k]))
}
message(sprintf("  imbalance (max/min pairs): %.3f", max(pairs_per_split) / min(pairs_per_split)))
