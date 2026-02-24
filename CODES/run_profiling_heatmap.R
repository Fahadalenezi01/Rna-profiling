# ==============================================================================
# R Script for Body Fluid Profiling using a One-Versus-Rest Heatmap
# ==============================================================================

# --- 1. Automated Package Installation ---
cat("--- Checking and installing necessary R packages ---\n")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
packages <- c("rtracklayer", "dplyr", "tidyr", "pheatmap", "org.Hs.eg.db")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing package:", pkg, "\n"))
    BiocManager::install(pkg)
  }
}
cat("--- All packages are installed. Loading libraries. ---\n")

# --- 2. Load Libraries ---
library(rtracklayer)
library(dplyr)
library(tidyr)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db)

# --- 3. Configuration & Sample Sheet ---
stringtie_dir_home <- "/home/alajdal/stringtie_quantification_results"
base_dir_d <- "/mnt/d"
output_dir <- "body_fluid_profiling_heatmap"
top_n_specific <- 25 # We will select the top 25 specific genes from each group

sample_sheet <- data.frame(
  sample = c("SL1", "SL2", "SL3", "SL3_rep", "SL4", "BL1", "BL2", "BL3", "BL4"),
  condition = c(rep("Saliva", 5), rep("Blood", 4)),
  filepath = c(
    file.path(base_dir_d, "SL1", "stringtie_output_SL1.gtf"),
    file.path(base_dir_d, "SL2", "stringtie_output_SL2.gtf"),
    file.path(base_dir_d, "SL3", "stringtie_output_SL3.gtf"),
    file.path(stringtie_dir_home, "SL3_rep.transcripts.gtf"),
    file.path(stringtie_dir_home, "SL4.transcripts.gtf"),
    file.path(base_dir_d, "BL1", "stringtie_output_BL1.gtf"),
    file.path(base_dir_d, "BL2", "stringtie_output_BL2.gtf"),
    file.path(base_dir_d, "BL3", "stringtie_output_BL3.gtf"),
    file.path(base_dir_d, "BL4", "stringtie_output_BL4.gtf")
  )
)

# --- 4. Create Output Directory ---
if (!dir.exists(output_dir)) dir.create(output_dir)

# ==============================================================================
# PART 1: Import and Prepare the Expression Matrix
# ==============================================================================
cat("--- PART 1: Importing and preparing TPM data ---\n")

all_tpm_data <- list()
for (i in 1:nrow(sample_sheet)) {
  sample_name <- sample_sheet$sample[i]
  file_path <- sample_sheet$filepath[i]
  if (!file.exists(file_path)) {
    warning(paste("File not found for", sample_name, "- Skipping."))
    next
  }
  
  gtf_df <- as.data.frame(rtracklayer::import(file_path))
  
  if ("TPM" %in% names(gtf_df) && "ref_gene_name" %in% names(gtf_df)) {
    gene_tpm <- gtf_df %>%
      filter(type == "transcript", !is.na(ref_gene_name)) %>%
      mutate(TPM = as.numeric(TPM)) %>%
      group_by(ref_gene_name) %>%
      summarise(Total_TPM = sum(TPM))
    
    colnames(gene_tpm)[colnames(gene_tpm) == "Total_TPM"] <- sample_name
    all_tpm_data[[sample_name]] <- gene_tpm
  }
}

# Merge all samples into a single TPM matrix
tpm_matrix <- Reduce(function(x, y) full_join(x, y, by = "ref_gene_name"), all_tpm_data)
rownames(tpm_matrix) <- tpm_matrix$ref_gene_name
tpm_matrix$ref_gene_name <- NULL
tpm_matrix[is.na(tpm_matrix)] <- 0

# UPDATED: Relaxed the filtering cutoff to be more inclusive.
# Keep genes with TPM > 1 in at least 2 samples.
keep <- rowSums(tpm_matrix > 1) >= 2
tpm_matrix_filtered <- tpm_matrix[keep, ]

# --- FINAL ROBUSTNESS CHECK 1 ---
if (nrow(tpm_matrix_filtered) == 0) {
  stop("Analysis stopped: No genes passed the initial expression filter (TPM > 1 in at least 2 samples).")
}

# Log2 transform the data for visualization (log2(TPM + 1))
log2_tpm_matrix <- log2(tpm_matrix_filtered + 1)

# ==============================================================================
# PART 2: One-Versus-Rest Analysis to Find Specific Genes
# ==============================================================================
cat("--- PART 2: Identifying body fluid-specific genes ---\n")

specific_genes <- c()
conditions <- unique(sample_sheet$condition)

for (cond in conditions) {
  cat(paste("  Finding top specific genes for:", cond, "\n"))
  
  in_group_samples <- sample_sheet$sample[sample_sheet$condition == cond]
  out_group_samples <- sample_sheet$sample[sample_sheet$condition != cond]
  
  mean_in_group <- rowMeans(log2_tpm_matrix[, in_group_samples, drop = FALSE])
  mean_out_group <- rowMeans(log2_tpm_matrix[, out_group_samples, drop = FALSE])
  
  specificity_score <- mean_in_group - mean_out_group
  
  top_specific <- names(sort(specificity_score, decreasing = TRUE))[1:top_n_specific]
  specific_genes <- c(specific_genes, top_specific)
}

specific_genes <- unique(specific_genes)

# ==============================================================================
# PART 3: Generate the Profiling Heatmap
# ==============================================================================
cat("--- PART 3: Generating the final profiling heatmap ---\n")

# --- FINAL ROBUSTNESS CHECK 2 ---
if (length(specific_genes) == 0 || all(is.na(specific_genes))) {
  cat("\n!!! NOTE: No specific genes were identified after one-vs-rest analysis. Cannot generate heatmap. !!!\n")
} else {
  
  # Create the final matrix for plotting with only the specific genes
  plot_matrix <- log2_tpm_matrix[specific_genes, ]

  # Convert Ensembl IDs to Gene Symbols for better readability
  ensembl_ids <- gsub("\\..*$", "", rownames(plot_matrix))
  gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  display_names <- ifelse(!is.na(gene_symbols), gene_symbols, rownames(plot_matrix))
  rownames(plot_matrix) <- display_names

  # Create an annotation data frame for the columns (samples)
  annotation_df <- data.frame(
    Condition = sample_sheet$condition,
    row.names = sample_sheet$sample
  )

  # Define colors for the annotation
  ann_colors <- list(
    Condition = c(Saliva = "blue", Blood = "red")
  )

  # --- Create and Save the Heatmap ---
  heatmap_file <- file.path(output_dir, "body_fluid_profiling_heatmap.tiff")

  pheatmap(
    plot_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_col = annotation_df,
    annotation_colors = ann_colors,
    main = "Body Fluid Gene Expression Profile (Top 25 Specific Genes)",
    fontsize_row = 8,
    scale = "row", # This is key: it shows relative expression patterns
    filename = heatmap_file,
    width = 8.5,  # Slightly wider for better layout
    height = 9.69,
    dpi = 600
  )

  cat(paste("  Profiling heatmap saved to:", heatmap_file, "\n"))
}

cat("\n--- Analysis complete! ---\n")
