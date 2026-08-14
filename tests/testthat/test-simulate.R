test_that("center_effect_size_matrix puts each gene's perturbed mean on its target effect size", {
  # Drawing per-guide effect sizes and clamping negatives to 0 biases the mean upwards, so the
  # centring step is what makes "effect size 0.15" mean a 15 % knockdown on average rather than
  # something slightly weaker.
  set.seed(1)
  n_genes <- 3
  n_cells <- 200
  gene_effect_sizes <- c(0.85, 0.90, 0.95)
  pert_status <- rep(c(1, 0), each = n_cells / 2)

  mat <- matrix(runif(n_genes * n_cells, 0.5, 1.2), nrow = n_genes)
  centred <- center_effect_size_matrix(mat, pert_status, gene_effect_sizes)

  expect_equal(rowMeans(centred[, pert_status == 1, drop = FALSE]), gene_effect_sizes)
  expect_equal(rowMeans(centred[, pert_status == 0, drop = FALSE]), rep(1, n_genes))
})

test_that("center_effect_size_matrix preserves shape and never goes negative", {
  set.seed(2)
  mat <- matrix(runif(4 * 50, 0, 0.2), nrow = 4)
  pert_status <- rep(c(1, 0), each = 25)
  centred <- center_effect_size_matrix(mat, pert_status, rep(0.1, 4))

  expect_equal(dim(centred), dim(mat))
  expect_true(all(centred >= 0))
})

test_that("create_effect_size_matrix returns genes x cells", {
  # Orientation is genes in rows, cells in columns. run_power_simulation.R relies on it directly --
  # `es_mat[, restore_cell_order]` reorders columns by cell, and center_effect_size_matrix subsets
  # columns by perturbation status. A transposed matrix would pair each cell with another gene's
  # effect, which is the class of bug that produced the size-factor shuffle.
  set.seed(3)
  n_cells <- 40
  gene_effect_sizes <- c(0.8, 0.9)
  # 0 = no guide, 1..2 = targeting guides, 3..4 = control guides
  grna_pert_status <- sample(0:4, n_cells, replace = TRUE)
  pert_guides <- c("g1", "g2")

  mat <- create_effect_size_matrix(grna_pert_status, pert_guides, gene_effect_sizes, guide_sd = 0.13)

  expect_equal(nrow(mat), length(gene_effect_sizes))
  expect_equal(ncol(mat), n_cells)
})

test_that("create_effect_size_matrix leaves guide-free cells unperturbed", {
  set.seed(4)
  # Cells 1 and 2 carry no guide; 3-4 carry targeting guides, 5-6 control guides.
  grna_pert_status <- c(0, 0, 1, 2, 3, 4)
  mat <- create_effect_size_matrix(grna_pert_status, c("g1", "g2"), c(0.8, 0.9), guide_sd = 0.13)

  # Row 1 of the guide table is the no-effect row, so a cell with no guide gets exactly 1 for
  # every gene. Cells are columns.
  expect_equal(unname(mat[, 1]), c(1, 1))
  expect_equal(unname(mat[, 2]), c(1, 1))
  expect_false(all(mat[, 3] == 1))
})

test_that("create_effect_size_matrix never returns a negative effect", {
  # guide_sd large enough that the normal draws go negative and must be clamped.
  set.seed(5)
  mat <- create_effect_size_matrix(sample(0:4, 100, replace = TRUE), c("g1", "g2"),
                                   c(0.5, 0.5), guide_sd = 2)
  expect_true(all(mat >= 0))
})

test_that("build_dispersion_vector inverts theta, in the order asked for", {
  precomp <- list(
    GENE_A = list(theta = 4),
    GENE_B = list(theta = 2),
    GENE_C = list(theta = 10)
  )
  disp <- build_dispersion_vector(precomp, c("GENE_C", "GENE_A"))

  expect_equal(names(disp), c("GENE_C", "GENE_A"))
  expect_equal(unname(disp), c(1 / 10, 1 / 4))
})

test_that("build_dispersion_vector errors rather than recycling when a gene is missing", {
  # The original bug: dispersions came back in a list column with NULL holes, and unlist() silently
  # shortened the vector so rnbinom recycled -- every gene after the first gap was simulated with
  # another gene's dispersion. A hard error is the point of this function.
  precomp <- list(GENE_A = list(theta = 4))
  expect_error(
    build_dispersion_vector(precomp, c("GENE_A", "GENE_MISSING")),
    "no entry in @response_precomputations"
  )
})

test_that("build_dispersion_vector errors on a non-finite dispersion", {
  expect_error(
    build_dispersion_vector(list(GENE_A = list(theta = 0)), "GENE_A"),
    "non-finite dispersion"
  )
  expect_error(
    build_dispersion_vector(list(GENE_A = list(theta = NA_real_)), "GENE_A"),
    "non-finite dispersion"
  )
})
