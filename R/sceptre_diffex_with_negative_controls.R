# Script: sceptre_diffex_with_negative_controls.R

### SETUP =====================================================================

# Saving image for debugging
if (!file.exists("RDA_objects/sceptre_diffex_with_negative_controls")) { dir.create("RDA_objects/sceptre_diffex_with_negative_controls", recursive = TRUE) }
save.image(paste0("RDA_objects/sceptre_diffex_with_negative_controls/", snakemake@wildcards$sample, ".rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

# Download Sceptre
# library(devtools)
# message("Installing Sceptre")
# devtools::install_github("katsevich-lab/sceptre")
# message("Sceptre Installation Complete")

message("Loading in packages")
suppressPackageStartupMessages({
  library(tidyverse)
  library(sceptre)
})

message("Loading input files")
sceptre_object <- readRDS(snakemake@input$sceptre_diffex_input)
gene_gRNA_group_pairs <- readRDS(snakemake@input$gene_gRNA_group_pairs)
gRNA_groups_table <- readRDS(snakemake@input$gRNA_groups_table)
negative_control_genes <- as.vector(snakemake@params$negative_control_genes)
features <- read_tsv(snakemake@input$features, col_names = FALSE) %>% filter(X3 == "Gene Expression")


### ADD NEGATIVE CONTROL PAIRS ================================================

# Which negative control genes aren't in the gene expression matrix
neg_control_genes_not_represented <- negative_control_genes[which(!negative_control_genes %in% features$X2)]
print(neg_control_genes_not_represented)
negative_control_genes <- setdiff(negative_control_genes, neg_control_genes_not_represented)

# For each negative control gene, match the gene up with each target
targets <- setdiff(unique(gRNA_groups_table$grna_target), "non-targeting") # Remove non-targeting
neg_control_gene_target_pairs <- expand.grid(grna_target = targets, response_id = negative_control_genes)

# Then convert the gene names to ensemble ids
neg_control_gene_id_target_pairs <- neg_control_gene_target_pairs %>% 
  left_join(features, by = c("response_id" = "X2")) %>% 
  select(grna_target, X1) %>%
  dplyr::rename(response_id = X1)

# Then add to gene_gRNA_group_pairs and remove duplicates
all_pairs <- rbind(gene_gRNA_group_pairs, neg_control_gene_id_target_pairs)
all_pairs <- unique(all_pairs)


### RUN DIFFEX ================================================================

# Specify all_pairs in Set analysis parameters
sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = all_pairs,
  side = "both",
  grna_integration_strategy = "union",
)

sceptre_object <- sceptre_object %>% 
  assign_grnas(method = "thresholding", threshold = 5) %>%
  run_qc() %>%
  run_calibration_check() %>%
  run_discovery_analysis(parallel = FALSE)


### SAVE OUTPUT ===============================================================

# Create "plots" subdirectory inside dirname(snakemake@output$discovery_results)
plots_dir <- file.path(dirname(snakemake@output$discovery_results), "plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Define ggsave wrapper for consistent usage
save_plot <- function(plot, file_name) {
  file_path <- file.path(plots_dir, paste0(file_name, ".pdf"))
  ggsave(
    filename = file_path,
    plot = plot,
    device = "pdf",
    width = 5,
    height = 5
  )
}

# Save plots with their respective names
message("Saving plots")
save_plot(plot_grna_count_distributions(sceptre_object), "plot_grna_count_distributions")
save_plot(plot_assign_grnas(sceptre_object), "plot_assign_grnas")
save_plot(plot_run_calibration_check(sceptre_object), "plot_run_calibration_check")
save_plot(plot_run_discovery_analysis(sceptre_object), "plot_run_discovery_analysis")
save_plot(plot_covariates(sceptre_object), "plot_covariates")
save_plot(plot_run_qc(sceptre_object), "plot_run_qc")

# Write all outputs to directory
message("Saving output files")
write_outputs_to_directory(sceptre_object, dirname(snakemake@output$discovery_results))

# Save the final sceptre object
saveRDS(sceptre_object, snakemake@output$final_sceptre_object)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)