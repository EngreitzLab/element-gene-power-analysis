## Simulate Perturb-seq counts with a specified effect size, including guide-to-guide
## variability in that effect size.
##
## Ported from R/power_simulations_fun.R. fit_negbinom_deseq2() is not carried over: nothing
## called it, and it was the only reason the pipeline depended on DESeq2.
##
## The numerical behaviour follows the original with one deliberate exception: the original drew
## a random permutation of the size factors before simulating, decoupling each cell's library
## size from its identity. That shuffle is wrong -- the size factor belongs to the cell whose
## counts are being simulated -- so it has been removed and each cell now keeps its own size
## factor. This also removes a sample() call from the RNG stream, so a given seed no longer
## reproduces the pre-refactor output draw for draw.

suppressPackageStartupMessages(library(Matrix))

#' Simulate counts for one perturbation.
#'
#' @param x sim_input, already subset to the genes being tested
#' @param effect_size_mat genes x cells multiplier, columns in the same order as x$cells
#' @return a dense genes x cells matrix of counts
draw_counts <- function(x, effect_size_mat) {
  gene_means <- x$row_data$mean
  gene_dispersions <- x$row_data$dispersion
  size_factors <- x$col_data$size_factors

  n_gene <- length(gene_means)
  n_cell <- length(size_factors)

  if (!identical(dim(effect_size_mat), c(n_gene, n_cell))) {
    stop("effect_size_mat is ", paste(dim(effect_size_mat), collapse = " x "),
         " but ", n_gene, " x ", n_cell, " was expected.", call. = FALSE)
  }

  # Cell-to-cell variability. Each cell keeps its own size factor: the effect-size matrix is
  # indexed by cell, so shuffling would pair one cell's perturbation status with another cell's
  # library size.
  #
  # mu[i, j] = gene_means[i] * size_factor[j] * effect_size[i, j].
  # outer() replaces the original matrix(rep(...)) + sweep() pair: same values, one allocation
  # instead of three.
  mu <- outer(gene_means, size_factors) * effect_size_mat

  # size = theta = 1 / dispersion. mu is consumed column-major, and `size` recycles over the
  # gene index, which is why gene_dispersions must be exactly n_gene long -- see
  # build_dispersion_vector() for the guard that guarantees it.
  counts <- rnbinom(n_cell * n_gene, mu = mu, size = 1 / gene_dispersions)

  matrix(counts, nrow = n_gene, ncol = n_cell, dimnames = list(x$genes, x$cells))
}

#' Convert simulated counts into the matrix class sceptre requires.
#'
#' sceptre needs CSR (`dgRMatrix`) here and there is no way around it: this pipeline assigns
#' `@response_matrix` directly, bypassing the normalisation `import_data()` would have applied,
#' and the consumers are unforgiving. `load_row()` (sceptre's s4_helpers.R) dispatches on `odm`
#' and `dgRMatrix` with *no else branch*, so handing it a dense matrix returns NULL silently
#' rather than raising an error; and with `run_permutations = FALSE` the CRT path reads
#' `response_matrix@j/@p/@x` directly, which only a dgRMatrix has.
#'
#' So density does not decide the representation here -- but it does decide the cost, and it is
#' worth logging. The original built a dgCMatrix via `Matrix(..., sparse = TRUE)` and then
#' re-built it as dgRMatrix with `as(..., "RsparseMatrix")`; this constructs the dgRMatrix in a
#' single pass instead. Rows are iterated because the gene count is small (single or double
#' digits) and iterating rows yields CSR's required within-row column ordering for free.
#'
#' @param counts dense genes x cells matrix from draw_counts()
#' @param report_density log the realised nonzero fraction
as_sceptre_response_matrix <- function(counts, report_density = TRUE) {
  n_gene <- nrow(counts)
  n_cell <- ncol(counts)

  nonzero_cols <- lapply(seq_len(n_gene), function(g) which(counts[g, ] != 0))
  per_row <- lengths(nonzero_cols)

  if (isTRUE(report_density)) {
    density <- sum(per_row) / (as.numeric(n_gene) * n_cell)
    message(sprintf(
      "Simulated matrix %d x %d, density %.1f%% (%s storage would be smaller; sceptre requires CSR)",
      n_gene, n_cell, 100 * density, if (density > 0.5) "dense" else "sparse"
    ))
  }

  values <- unlist(
    lapply(seq_len(n_gene), function(g) counts[g, nonzero_cols[[g]]]),
    use.names = FALSE
  )

  new("dgRMatrix",
      j = as.integer(unlist(nonzero_cols, use.names = FALSE)) - 1L,
      p = c(0L, as.integer(cumsum(per_row))),
      x = as.numeric(values),
      Dim = c(n_gene, n_cell),
      Dimnames = list(rownames(counts), colnames(counts)))
}

#' Build the per-gene dispersion vector from sceptre's cached precomputations.
#'
#' The original wrote `rowData(sce)$dispersion <- dispersion_values[rownames(sce)]`, which
#' produced a *list* column with a NULL entry for any gene sceptre had no precomputation for.
#' `unlist()` then dropped those NULLs, shortening the vector, and rnbinom() recycled it -- so
#' every gene after the first gap would have been simulated with another gene's dispersion, with
#' no warning. On sample1 this never fired (all 239 tested genes have precomputations, though 19
#' of the 292 genes in the response matrix do not), but the failure is silent, so it is worth a
#' hard guard rather than a comment.
#'
#' @param precomputations sceptre_object@response_precomputations
#' @param genes genes to build the vector for
#' @return named numeric vector, one finite value per gene, in the order given
build_dispersion_vector <- function(precomputations, genes) {
  missing_genes <- setdiff(genes, names(precomputations))
  if (length(missing_genes) > 0) {
    stop(length(missing_genes), " gene(s) have no entry in @response_precomputations and so no ",
         "dispersion estimate, including: ",
         paste(utils::head(missing_genes, 5), collapse = ", "),
         ". Re-run sceptre's precomputation, or drop these genes from the discovery pairs.",
         call. = FALSE)
  }

  theta <- vapply(precomputations[genes], function(p) p$theta, numeric(1))
  dispersion <- 1 / theta

  if (any(!is.finite(dispersion))) {
    bad <- genes[!is.finite(dispersion)]
    stop(length(bad), " gene(s) have a non-finite dispersion (theta of 0 or NA), including: ",
         paste(utils::head(bad, 5), collapse = ", "), ".", call. = FALSE)
  }

  stats::setNames(dispersion, genes)
}

## GUIDE-LEVEL VARIABILITY =========================================================================

#' Pick one expressed guide at random per cell, from a CSC perturbation matrix.
#'
#' Returns 0-based row indices, or -1 for a cell with no guide (the caller adds 1, mapping it to
#' the no-effect row).
sample_guide <- function(pert_status) {
  num_cols <- dim(pert_status)[[2]]
  return_vector <- integer(num_cols)

  for (col_idx in seq_len(num_cols)) {
    start_idx <- pert_status@p[col_idx] + 1
    end_idx <- pert_status@p[col_idx + 1]
    column_non_zeros <- end_idx - (start_idx - 1)

    if (column_non_zeros > 0) {
      selected_i <- pert_status@i[start_idx + sample(0:(column_non_zeros - 1), 1)]
      return_vector[col_idx] <- selected_i
    } else {
      return_vector[col_idx] <- -1L
    }
  }

  return_vector
}

#' Collapse a perturbation matrix to one status per cell, choosing randomly among multiples.
convert_pert_mat_to_vector <- function(pert_mat) {
  sample_guide(pert_mat) + 1
}

#' The cell ordering create_guide_pert_status() produces: perturbed cells then control cells.
#'
#' `es_mat[, order(cell_order(pert_status))]` restores x$cells order. Computed once per target
#' rather than once per rep -- the original re-derived it by name inside the rep loop.
cell_order <- function(pert_status) {
  c(which(pert_status == 1), which(pert_status == 0))
}

#' Per-cell gRNA perturbation status: which guide (if any) each cell carries.
#'
#' Control-cell statuses are offset past the targeting guides so that every guide, targeting or
#' not, has a distinct index into the effect-size table.
#'
#' IMPORTANT: the result is ordered perturbed-cells-first-then-control-cells (each block in
#' ascending cell position), which is *not* the order of x$cells. Use cell_order() to get the
#' permutation and reorder the effect-size matrix back. The original relied on cell barcodes as
#' names for this; positions are cheaper and cannot silently mismatch.
create_guide_pert_status <- function(pert_status, grna_perts, pert_guides) {
  grnas_pert_cells <- grna_perts[pert_guides, pert_status == 1]
  if (!is.null(nrow(grnas_pert_cells))) {
    grnas_pert_cells <- convert_pert_mat_to_vector(grnas_pert_cells)
  }

  grnas_ctrl_cells <- grna_perts[!rownames(grna_perts) %in% pert_guides, pert_status == 0]
  if (!is.null(nrow(grnas_ctrl_cells))) {
    grnas_ctrl_cells <- convert_pert_mat_to_vector(grnas_ctrl_cells)
  }

  ctrl_perts <- grnas_ctrl_cells > 0
  grnas_ctrl_cells[ctrl_perts] <- grnas_ctrl_cells[ctrl_perts] + max(grnas_pert_cells)

  c(grnas_pert_cells, grnas_ctrl_cells)
}

#' Effect-size matrix with guide-to-guide variability.
#'
#' Each guide gets its own effect size drawn around the target effect size (targeting guides) or
#' around 1 (control guides), with standard deviation guide_sd. Negative draws are clamped to 0.
create_effect_size_matrix <- function(grna_pert_status, pert_guides, gene_effect_sizes, guide_sd) {
  n_pert_guides <- length(pert_guides)
  n_ctrl_guides <- max(grna_pert_status) - n_pert_guides

  guide_effect_sizes_pert <- vapply(gene_effect_sizes, FUN = rnorm, n = n_pert_guides,
                                    sd = guide_sd, FUN.VALUE = numeric(n_pert_guides))
  guide_effect_sizes_ctrl <- vapply(rep(1, length(gene_effect_sizes)), FUN = rnorm,
                                    n = n_ctrl_guides, sd = guide_sd,
                                    FUN.VALUE = numeric(n_ctrl_guides))
  guide_effect_sizes <- rbind(guide_effect_sizes_pert, guide_effect_sizes_ctrl)
  guide_effect_sizes[guide_effect_sizes < 0] <- 0

  # Row 1 is the no-effect row, used by cells carrying no guide.
  guide_effect_sizes <- rbind(1, guide_effect_sizes)

  t(guide_effect_sizes[grna_pert_status + 1, ])
}

#' Shift the effect-size matrix so each gene's mean effect equals the requested effect size.
#'
#' Drawing per-guide effect sizes and clamping negatives to 0 biases the mean upwards, so the
#' perturbed and control blocks are each re-centred: perturbed on gene_effect_sizes, control on 1.
center_effect_size_matrix <- function(effect_size_mat, pert_status, gene_effect_sizes) {
  is_pert <- pert_status == 1
  is_ctrl <- pert_status == 0

  mean_es_pert <- rowMeans(effect_size_mat[, is_pert, drop = FALSE])
  mean_es_ctrl <- rowMeans(effect_size_mat[, is_ctrl, drop = FALSE])

  effect_size_mat[, is_pert] <- effect_size_mat[, is_pert] + (gene_effect_sizes - mean_es_pert)
  effect_size_mat[, is_ctrl] <- effect_size_mat[, is_ctrl] + (1 - mean_es_ctrl)

  effect_size_mat[effect_size_mat < 0] <- 0
  effect_size_mat
}
