# Script: create_sce_object.R

### SETUP =====================================================================

# Saving image for debugging
if (!file.exists("RDA_objects/create_sce_object")) { dir.create("RDA_objects/create_sce_object", recursive = TRUE) }
save.image(paste0("RDA_objects/create_sce_object/", snakemake@wildcards$sample, ".rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading in packages")
# required packages and functions
suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
})
message("Loading input files")

# load in the final sceptre object
final_sceptre_object <- readRDS(snakemake@input$final_sceptre_object)

make_assignment_sparse_matrix <- function(assignments, n_cols, row_names, col_names,
                                          cell_index = NULL) {
  n_rows <- length(assignments)
  assignment_lengths <- lengths(assignments)
  n_assignments <- sum(assignment_lengths)

  if (n_assignments == 0L) {
    return(sparseMatrix(
      i = integer(0),
      j = integer(0),
      x = numeric(0),
      dims = c(n_rows, n_cols),
      dimnames = list(row_names, col_names)
    ))
  }

  row_idx <- rep.int(seq_len(n_rows), assignment_lengths)
  col_idx <- unlist(assignments, use.names = FALSE)
  if (!is.null(cell_index)) {
    mapped_col_idx <- rep(NA_integer_, length(col_idx))
    valid_source_idx <- !is.na(col_idx) & col_idx >= 1L & col_idx <= length(cell_index)
    mapped_col_idx[valid_source_idx] <- cell_index[col_idx[valid_source_idx]]
    col_idx <- mapped_col_idx
  }

  keep_idx <- !is.na(col_idx) & col_idx >= 1L & col_idx <= n_cols
  if (any(!keep_idx)) {
    message("Dropping ", sum(!keep_idx), " assignment indices outside the SCE column range.")
  }

  output <- sparseMatrix(
    i = row_idx[keep_idx],
    j = as.integer(col_idx[keep_idx]),
    x = rep.int(1, sum(keep_idx)),
    dims = c(n_rows, n_cols),
    dimnames = list(row_names, col_names)
  )
  if (length(output@x) > 0L) {
    output@x[] <- 1
  }
  output
}

estimate_poscounts_size_factors_sparse <- function(counts, locfunc = stats::median) {
  if (!inherits(counts, "sparseMatrix")) {
    stop("Expected a sparse Matrix object for count normalization. ",
         "DESeq2::DESeqDataSetFromMatrix() would coerce this matrix to dense.")
  }

  counts <- as(counts, "CsparseMatrix")
  if (any(counts@x < 0)) {
    stop("Count matrix contains negative values.")
  }
  counts <- drop0(counts)

  log_counts <- counts
  log_counts@x <- log(log_counts@x)
  log_geomeans <- Matrix::rowSums(log_counts) / ncol(counts)
  log_geomeans[Matrix::rowSums(counts) == 0] <- -Inf
  finite_geomeans <- is.finite(log_geomeans)

  size_factors <- vapply(seq_len(ncol(counts)), function(col_idx) {
    start_idx <- counts@p[col_idx] + 1L
    end_idx <- counts@p[col_idx + 1L]
    if (start_idx > end_idx) {
      return(NA_real_)
    }

    entry_idx <- seq.int(start_idx, end_idx)
    row_idx <- counts@i[entry_idx] + 1L
    count_values <- counts@x[entry_idx]
    keep_idx <- finite_geomeans[row_idx] & count_values > 0
    if (!any(keep_idx)) {
      return(NA_real_)
    }

    exp(locfunc(log(count_values[keep_idx]) - log_geomeans[row_idx[keep_idx]]))
  }, numeric(1))

  if (any(!is.finite(size_factors)) || any(size_factors <= 0)) {
    stop("Could not compute finite positive size factors for all cells.")
  }

  size_factors
}

sparse_normalized_row_means <- function(counts, size_factors) {
  counts <- as(counts, "CsparseMatrix")
  counts@x <- counts@x / rep.int(size_factors, diff(counts@p))
  Matrix::rowSums(counts) / ncol(counts)
}


### CREATE SCE OBJECT ========================================================

expr <- final_sceptre_object@response_matrix[[1]]
if (!inherits(expr, "sparseMatrix")) {
  stop("final_sceptre_object@response_matrix[[1]] must be a sparse Matrix object. ",
       "A dense gene-by-cell count matrix is too large for this pipeline.")
}
expr <- as(expr, "CsparseMatrix")

if ("batch" %in% colnames(final_sceptre_object@covariate_data_frame)) {
  cell_metadata <- data.frame(
    cell_barcode = rownames(final_sceptre_object@covariate_data_frame),
    cell_batches = final_sceptre_object@covariate_data_frame$batch
  )
  sce <- SingleCellExperiment(assays = list(counts = expr), colData = cell_metadata)
  colnames(sce) <-  rownames(final_sceptre_object@covariate_data_frame)
} else {
  cell_metadata <- NULL
  sce <- SingleCellExperiment(assays = list(counts = expr))
  colnames(sce) <-  rownames(final_sceptre_object@covariate_data_frame)
}


### ADD SCEPTRE GRNA_PERTS ASSIGNMENTS TO SCE =================================

# Add the individual grna perts first
message("Adding the individual grna_perts")
get_grna_assignments <- function(sceptre_object) {
  if (!sceptre_object@functs_called[["assign_grnas"]]) {
    stop("`assign_grnas()` has not yet been called on the `sceptre_object`.")
  }
  return(sceptre_object@initial_grna_assignment_list)
}

individual_grna_assignments <- get_grna_assignments(final_sceptre_object)

# Number of rows and columns for the matrix
nRows <- length(individual_grna_assignments)
nCols <- length(colnames(sce))

# Build perturbation status matrix in one sparse allocation.
sparseMatrix <- make_assignment_sparse_matrix(
  assignments = individual_grna_assignments,
  n_cols = nCols,
  row_names = names(individual_grna_assignments),
  col_names = colnames(sce)
)

# Add to sce object
message("Adding grna_perts to sce object")
altExp(sce, e = "grna_perts") <- SummarizedExperiment(assays = list(perts = sparseMatrix))


### ADD SCEPTRE CRE_PERTS ASSIGNMENTS TO SCE ==================================

# Now add the individual "cre_perts"
message("Adding the individual cre_perts")
get_cre_assignments <- function(sceptre_object) {
  if (!sceptre_object@functs_called[["assign_grnas"]]) {
    stop("`assign_grnas()` has not yet been called on the `sceptre_object`.")
  }
  return(sceptre_object@grna_assignments$grna_group_idxs)
}

individual_cre_assignments <- get_cre_assignments(final_sceptre_object)

# Number of rows and columns for the matrix
nRows <- length(individual_cre_assignments)
nCols <- length(colnames(sce))

# Indices from individual_cre_assignments point to cells_in_use.
sparseMatrix <- make_assignment_sparse_matrix(
  assignments = individual_cre_assignments,
  n_cols = nCols,
  row_names = names(individual_cre_assignments),
  col_names = colnames(sce),
  cell_index = final_sceptre_object@cells_in_use
)

# Add to sce object
message("Adding cre_prets to sce object")
altExp(sce, e = "cre_perts") <- SummarizedExperiment(assays = list(perts = sparseMatrix))


### ADD SCEPTRE DISPERION ESTIMATES ===========================================

dispersion_values <- lapply(final_sceptre_object@response_precomputations, function(x) 1/x$theta)

# Add to rowData
rowData(sce)$dispersion <- dispersion_values[rownames(sce)]


### CALCULATE SIZE FACTORS ====================================================

# Calculate total_umis and detected_genes for Deseq2 object creation
response_matrix <- assay(sce, "counts")
coldata <- data.frame(
  total_umis = colSums(response_matrix),
  detected_genes = colSums(response_matrix > 0)
)
rowData(sce)[, "average_expression_all_cells"] <- Matrix::rowMeans(response_matrix)

# Compute DESeq2-style poscounts size factors without coercing sparse counts to dense.
size_factors <- estimate_poscounts_size_factors_sparse(response_matrix)
colData(sce)[, "total_umis"] <- coldata$total_umis
colData(sce)[, "detected_genes"] <- coldata$detected_genes
colData(sce)[, "size_factors"] <- size_factors
rowData(sce)[, "mean"] <- sparse_normalized_row_means(response_matrix, size_factors)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
saveRDS(sce, snakemake@output$perturb_sce)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)
