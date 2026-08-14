## Build the per-perturbation input for one power simulation.
##
## Ported from R/differential_expression_fun.R (pert_input at :444, pert_input_sampled at
## :464), which were the only two of that file's sixteen functions the pipeline ever called.
## The other fourteen -- MAST/DEsingle/LFC differential expression, normalisation, cell
## filtering, distance-based gene filtering and output annotation -- were unreachable and are
## not carried over.
##
## Behaviour is unchanged apart from operating on a sim_input list instead of a
## SingleCellExperiment.

#' Use every non-perturbed cell as the control group.
#'
#' @param x sim_input
#' @param target perturbation target (a row of x$perts[[level]])
#' @param level which perturbation level to read the status from
pert_input <- function(x, target, level = "cre_perts") {
  message("Creating input for perturbation '", target, "'.")
  add_pert_status(x, target, level)
}

#' Sample a fixed number of control cells instead of using all of them.
#'
#' Much cheaper than pert_input() when most cells are unperturbed, since every downstream
#' matrix is (genes x cells). Opt-in and experimental: the caller must also repair the
#' sceptre object's cell indexing, which pert_input() does not require. See
#' run_power_simulation.R.
#'
#' @param n_ctrl number of control cells to draw
#' @param cell_batches name of a x$col_data column to stratify the draw by, or NULL
pert_input_sampled <- function(x, target, level = "cre_perts", n_ctrl, cell_batches = NULL) {
  message("Creating input for perturbation '", target, "' with ", n_ctrl, " sampled control cells.")

  pert_matrix <- x$perts[[level]]
  row <- match(target, rownames(pert_matrix))
  if (is.na(row)) {
    stop("Target '", target, "' is not a row of perts$", level, ".", call. = FALSE)
  }
  # Columns are positionally aligned to x$cells, so this is already in cell order.
  is_pert <- as.vector(pert_matrix[row, ]) > 0

  n_available_ctrl <- sum(!is_pert)
  if (n_ctrl > n_available_ctrl) {
    stop("Requested ", n_ctrl, " control cells for target '", target, "' but only ",
         n_available_ctrl, " non-perturbed cells exist.", call. = FALSE)
  }

  if (!is.null(cell_batches)) {
    if (!cell_batches %in% colnames(x$col_data)) {
      stop("cell_batches column '", cell_batches, "' is not present in col_data. Available: ",
           paste(colnames(x$col_data), collapse = ", "), ".", call. = FALSE)
    }

    # Draw controls with the same batch composition as the perturbed cells. A perturbation
    # localised to one batch can leave too few control cells in that batch to satisfy the
    # quota, so fall back to an unstratified draw rather than failing the whole task.
    ctrl_idx <- tryCatch({
      message("Sampling control cells with equal proportions as perturbed cells across '",
              cell_batches, "'.")
      batch_of_cell <- x$col_data[[cell_batches]]
      pert_batch_prop <- table(batch_of_cell[is_pert]) / sum(is_pert)
      ctrl_per_batch <- round(pert_batch_prop * n_ctrl)

      drawn <- lapply(names(ctrl_per_batch), function(batch) {
        eligible <- which(!is_pert & batch_of_cell == batch)
        sample(eligible, size = ctrl_per_batch[[batch]], replace = FALSE)
      })
      unlist(drawn, use.names = FALSE)
    }, error = function(e) {
      message("Error occurred in cell batch sampling for differential expression.")
      message("Most likely the perturbation being tested is localized to one batch, so there ",
              "aren't enough control cells to sample within that batch.")
      message("The specific error is: ", conditionMessage(e))
      message("Sampling control cells irrespective of batch instead.")
      sample(which(!is_pert), size = n_ctrl, replace = FALSE)
    })
  } else {
    ctrl_idx <- sample(which(!is_pert), size = n_ctrl, replace = FALSE)
  }

  # Perturbed cells first, then the sampled controls. run_power_simulation.R relies on this
  # ordering when it rewrites the sceptre object's grna_group_idxs.
  x <- subset_cells(x, c(which(is_pert), ctrl_idx))
  add_pert_status(x, target, level)
}
