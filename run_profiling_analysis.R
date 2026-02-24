# ==============================================================================
# R Script for Exploratory Data Analysis based on Top Variable Genes
# (Method adapted from user's successful heatmap_PCA.R script)
# ==============================================================================

# --- 1. Automated Package Installation ---
cat("--- Checking and installing necessary R packages ---\n")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
packages <- c("rtracklayer", "dplyr", "tidyr", "pheatmap", "ggplot2", "ggrepel", "org.Hs.eg.db")
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
library(ggplot2)
library(ggrepel)
library(AnnotationDbi)
library(org.Hs.eg.db)

# --- 3. Configuration & Sample Sheet ---
stringtie_dir_home <- "/home/alajdal/stringtie_quantification_results"
base_dir_d <- "/mnt/d"
output_dir <- "exploratory_analysis_variable_genes"
top_n_variable <- 50 # Number of top variable genes to use

sample_sheet <- data.frame(
  sample = c(
    "SL1", "SL2", "SL3_rep", "SL4",
    "BL1", "BL2", "BL3", "BL4",
    "50SL50BL", "70SL30BL", "30SL70BL"
  ),
  condition = c(
    rep("Saliva", 4),
    rep("Blood", 4),
    rep("Mixed", 3)
  ),
  filepath = c(
    file.path(base_dir_d, "SL1", "stringtie_output_SL1.gtf"),
    file.path(base_dir_d, "SL2", "stringtie_output_SL2.gtf"),
    file.path(stringtie_dir_home, "SL3_rep.transcripts.gtf"),
    file.path(stringtie_dir_home, "SL4.transcripts.gtf"),
    file.path(base_dir_d, "BL1", "stringtie_output_BL1.gtf"),
    file.path(base_dir_d, "BL2", "stringtie_output_BL2.gtf"),
    file.path(base_dir_d, "BL3", "stringtie_output_BL3.gtf"),
    file.path(base_dir_d, "BL4", "stringtie_output_BL4.gtf"),
    file.path(stringtie_dir_home, "50SL50BL.transcripts.gtf"),
    file.path(stringtie_dir_home, "70SL30BL.transcripts.gtf"),
    file.path(stringtie_dir_home, "30SL70BL.transcripts.gtf")
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

# ==============================================================================
# PART 2: Select Top Variable Genes
# ==============================================================================
cat("--- PART 2: Selecting top variable genes ---\n")

# Filter for expressed genes: Keep genes with a mean TPM > 2 across all samples
tpm_matrix_filtered <- tpm_matrix[rowMeans(tpm_matrix) > 2, ]

# Calculate the variance for each gene
tpm_matrix_filtered$variance <- apply(tpm_matrix_filtered, 1, var)

# Get the top N most variable genes
top_variable_genes <- tpm_matrix_filtered %>%
  arrange(desc(variance)) %>%
  head(top_n_variable) %>%
  select(-variance) # Remove the variance column

# Log2 transform the data for visualization
log2_tpm_matrix <- log2(top_variable_genes + 1)

# ==============================================================================
# PART 3: Generate Heatmap and PCA Plot
# ==============================================================================
cat("--- PART 3: Generating visualizations ---\n")

# --- Convert Ensembl IDs to Gene Symbols for better readability ---
ensembl_ids <- gsub("\\..*$", "", rownames(log2_tpm_matrix))
gene_symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
display_names <- ifelse(!is.na(gene_symbols), gene_symbols, rownames(log2_tpm_matrix))
rownames(log2_tpm_matrix) <- display_names

# --- Create Annotation for Plots ---
annotation_df <- data.frame(
  Condition = factor(sample_sheet$condition),
  row.names = sample_sheet$sample
)
ann_colors <- list(
  Condition = c(Saliva = "blue", Blood = "red", Mixed = "darkgreen")
)

# --- Generate Heatmap ---
heatmap_file <- file.path(output_dir, "top_50_variable_genes_heatmap.tiff")
pheatmap(
  log2_tpm_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  main = paste("Top", top_n_variable, "Most Variable Genes"),
  fontsize_row = 8,
  scale = "row",
  filename = heatmap_file,
  width = 8.5,
  height = 9.69,
  dpi = 600
)
cat(paste("  Heatmap saved to:", heatmap_file, "\n"))

# --- Generate PCA Plot ---
pca_results <- prcomp(t(log2_tpm_matrix), scale. = TRUE)
pca_data <- as.data.frame(pca_results$x)
pca_data$condition <- sample_sheet$condition[match(rownames(pca_data), sample_sheet$sample)]
pca_data$sample <- rownames(pca_data)

percentVar <- round(100 * pca_results$sdev^2 / sum(pca_results$sdev^2))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text_repel(aes(label = sample), show.legend = FALSE, size = 4) +
  scale_color_manual(values = c("Saliva" = "blue", "Blood" = "red", "Mixed" = "darkgreen")) +
  labs(
    title = "PCA of Top 50 Most Variable Genes",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "bottom")

pca_plot_file <- file.path(output_dir, "top_50_variable_genes_PCA_plot.tiff")
ggsave(
  pca_plot_file,
  plot = pca_plot,
  device = "tiff",
  width = 8.5,
  height = 8.5,
  units = "in",
  dpi = 600
)
cat(paste("  PCA plot saved to:", pca_plot_file, "\n"))

cat("\n--- Analysis complete! ---\n")
