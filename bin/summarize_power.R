#!/usr/bin/env Rscript
#
# Merge the per-effect-size power tables into one wide table, one row per element-gene pair.
#
# Without this the pipeline leaves you with a separate file per effect size and no combined view.
# This recovers the useful half of the old format_sceptre_output.R -- a single table with one
# power column per effect size -- without that script's dependencies on `distances`,
# `guide_targets` and `features` files, none of which the pipeline produces.
#
# It also derives the quantity the analysis is usually actually after: the smallest tested effect
# size at which a pair reaches a target power. See --power-threshold.
#
# Usage:
#   summarize_power.R --power power_es0.15.tsv,power_es0.2.tsv --out power_summary.tsv

local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  lib <- if (length(file_arg) == 1) {
    file.path(dirname(normalizePath(sub("^--file=", "", file_arg))), "lib")
  } else {
    "bin/lib"
  }
  source(file.path(lib, "cli.R"))
})

## ARGUMENTS =======================================================================================

option_list <- list(
  make_option("--power", type = "character", default = NULL, dest = "power",
              help = "Comma-separated compute_power.R outputs, one per effect size."),
  make_option("--power-threshold", type = "double", default = 0.8, dest = "power_threshold",
              help = paste("Power level used to derive min_detectable_effect_size: the smallest",
                           "tested effect size at which a pair reaches it [default %default].")),
  make_option("--out", type = "character", default = NULL, dest = "out",
              help = "Output TSV, one row per pair.")
)

opts <- parse_args(OptionParser(
  option_list = option_list,
  description = "Merge per-effect-size power tables into one wide table."
))
require_options(opts, c("power", "out"))

if (opts$power_threshold <= 0 || opts$power_threshold > 1) {
  stop("--power-threshold must be in (0, 1].", call. = FALSE)
}

## LOAD ============================================================================================

paths <- trimws(strsplit(opts$power, ",", fixed = TRUE)[[1]])
required <- c("grna_target", "response_id", "power", "effect_size")

tables <- lapply(paths, function(path) {
  df <- read_tsv_file(path, required_columns = required)
  effect_sizes <- unique(df$effect_size)
  if (length(effect_sizes) != 1) {
    stop(path, " contains ", length(effect_sizes), " effect sizes; expected exactly one.",
         call. = FALSE)
  }
  df
})
effect_sizes <- vapply(tables, function(df) df$effect_size[1], numeric(1))

if (anyDuplicated(effect_sizes)) {
  stop("Two or more inputs share an effect size (",
       paste(effect_sizes[duplicated(effect_sizes)], collapse = ", "), ").", call. = FALSE)
}

# Ascending, so min_detectable_effect_size below can take the first that clears the threshold and
# the columns read left to right from weakest to strongest perturbation.
ordering <- order(effect_sizes)
tables <- tables[ordering]
effect_sizes <- effect_sizes[ordering]
log_step("Merging ", length(tables), " effect size(s): ", paste(effect_sizes, collapse = ", "))

#' Column suffix for an effect size, as a percentage: 0.15 -> "15", giving the column name
#' power_at_effect_size_15.
#'
#' Trailing zeros are trimmed only *after a decimal point*, so 0.125 -> "12.5" rather than
#' "12.50". Trimming unconditionally is wrong and quietly so: it would turn 0.2 into "2" and 0.5
#' into "5", i.e. a column called power_at_effect_size_2 for a 20 % knockdown.
#'
#' The decimal point then becomes an underscore (12.5 -> "12_5") to keep every column name
#' snake_case and free of dots, which several downstream readers treat as name separators.
effect_label <- function(effect_size) {
  formatted <- format(effect_size * 100, trim = TRUE, scientific = FALSE)
  if (grepl(".", formatted, fixed = TRUE)) {
    formatted <- sub("0+$", "", formatted)
    formatted <- sub("\\.$", "", formatted)
  }
  gsub(".", "_", formatted, fixed = TRUE)
}

## MERGE ===========================================================================================

# Outer join on the pair, so a pair tested at only some effect sizes still appears (with NA
# elsewhere) rather than being silently dropped.
pairs <- unique(do.call(rbind, lapply(tables, function(df) df[, c("grna_target", "response_id")])))
pairs <- pairs[order(pairs$grna_target, pairs$response_id), , drop = FALSE]
rownames(pairs) <- NULL
key_of <- function(df) paste(df$grna_target, df$response_id, sep = "\r")
pair_keys <- key_of(pairs)

out <- pairs

# Columns shared across effect sizes: taken from the first table that has the pair. Flagged if they
# disagree between effect sizes, which would mean the inputs do not describe the same experiment.
shared <- c("mean_pert_cells", "average_expression_all_cells")
for (column in shared) {
  if (!any(vapply(tables, function(df) column %in% colnames(df), logical(1)))) next
  values <- matrix(NA_real_, nrow = nrow(pairs), ncol = length(tables))
  for (i in seq_along(tables)) {
    if (!column %in% colnames(tables[[i]])) next
    values[, i] <- tables[[i]][[column]][match(pair_keys, key_of(tables[[i]]))]
  }
  spread <- apply(values, 1, function(v) {
    v <- v[!is.na(v)]
    if (length(v) < 2) 0 else diff(range(v)) / max(abs(v), 1e-12)
  })
  if (any(spread > 1e-6, na.rm = TRUE)) {
    message("  note: '", column, "' differs between effect sizes for ",
            sum(spread > 1e-6, na.rm = TRUE), " pair(s); reporting the first non-missing value.")
  }
  out[[column]] <- apply(values, 1, function(v) if (all(is.na(v))) NA_real_ else v[!is.na(v)][1])
}

# Per-effect-size columns.
power_columns <- character(0)
for (i in seq_along(tables)) {
  df <- tables[[i]]
  idx <- match(pair_keys, key_of(df))
  label <- effect_label(effect_sizes[i])
  base <- paste0("power_at_effect_size_", label)
  power_columns <- c(power_columns, base)

  out[[base]] <- df$power[idx]
  for (from in c("power_ci_low", "power_ci_high", "n_reps")) {
    if (from %in% colnames(df)) {
      suffix <- switch(from, power_ci_low = "_ci_low", power_ci_high = "_ci_high",
                       n_reps = "_n_reps")
      out[[paste0(base, suffix)]] <- df[[from]][idx]
    }
  }
}

## DERIVED =========================================================================================

# Smallest tested effect size reaching --power-threshold. NA means no tested effect size did, which
# is a statement about the effect sizes you ran, not proof the pair is undetectable -- so the
# largest effect size tested is reported alongside it for context.
power_matrix <- as.matrix(out[, power_columns, drop = FALSE])
out$min_detectable_effect_size <- apply(power_matrix, 1, function(row) {
  reached <- which(!is.na(row) & row >= opts$power_threshold)
  if (length(reached) == 0) NA_real_ else effect_sizes[min(reached)]
})
out$max_effect_size_tested <- max(effect_sizes)

out <- out[order(out$min_detectable_effect_size, na.last = TRUE,
                 decreasing = FALSE), , drop = FALSE]
rownames(out) <- NULL

## WRITE ===========================================================================================

write_tsv_file(out, opts$out)
log_step("Wrote ", nrow(out), " pairs x ", length(effect_sizes), " effect size(s) to ", opts$out)

message(sprintf("  reaching power >= %.2f at some tested effect size: %d of %d (%.0f%%)",
                opts$power_threshold, sum(!is.na(out$min_detectable_effect_size)), nrow(out),
                100 * mean(!is.na(out$min_detectable_effect_size))))
for (i in seq_along(effect_sizes)) {
  column <- power_columns[i]
  values <- out[[column]]
  message(sprintf("  %-24s mean power %.3f | >= %.2f in %d of %d pairs",
                  column, mean(values, na.rm = TRUE), opts$power_threshold,
                  sum(values >= opts$power_threshold, na.rm = TRUE), sum(!is.na(values))))
}
