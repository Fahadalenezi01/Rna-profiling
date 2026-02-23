# ==============================================================================
# R Script for UMAP Analysis of Saliva, Blood, and Mixed Samples
# ==============================================================================

# --- 1. Automated Package Installation ---
# This section checks for required packages and installs them if they are missing.
cat("--- Checking and installing necessary R packages ---\n")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# List of required packages
packages <- c("rtracklayer", "dplyr", "tidyr", "ggplot2", "ggrepel", "uwot", "DESeq2")

# Loop through packages and install if missing
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing package:", pkg, "\n"))
    if (pkg %in% c("rtracklayer", "DESeq2")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}
cat("--- All packages are installed. Loading libraries. ---\n")

# --- 2. Load Libraries ---
library(rtracklayer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(uwot)
library(DESeq2)

# --- 3. Configuration & Sample Sheet ---
stringtie_dir_home <- "/home/alajdal/stringtie_quantification_results"
base_dir_d <- "/mnt/d"
output_dir <- "umap_analysis_results"

# Sample sheet includes all samples for the analysis
sample_sheet <- data.frame(
  sample = c(
    "SL1", "SL2", "SL3_rep", "SL4",             # Pure Saliva
    "BL1", "BL2", "BL3", "BL4",                 # Pure Blood
    "50SL50BL", "70SL30BL", "30SL70BL"           # Mixed Samples
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
cat("--- PART 1: Importing and preparing TPM data for all samples ---\n")

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
      summarise(Total_TPM = sum(TPM, na.rm = TRUE))
    
    colnames(gene_tpm)[colnames(gene_tpm) == "Total_TPM"] <- sample_name
    all_tpm_data[[sample_name]] <- gene_tpm
  }
}

tpm_matrix <- Reduce(function(x, y) full_join(x, y, by = "ref_gene_name"), all_tpm_data)
tpm_matrix <- as.data.frame(tpm_matrix)
rownames(tpm_matrix) <- tpm_matrix$ref_gene_name
tpm_matrix$ref_gene_name <- NULL
tpm_matrix[is.na(tpm_matrix)] <- 0

# Convert to a count matrix for DESeq2 transformation
count_matrix <- round(tpm_matrix)

# Filter for expressed genes
keep <- rowSums(count_matrix >= 10) >= 2
count_matrix_filtered <- count_matrix[keep, ]

# ==============================================================================
# PART 2: Perform Data Transformation and UMAP
# ==============================================================================
cat("--- PART 2: Performing data transformation and UMAP ---\n")

# Create a DESeq2 object for data transformation
dds <- DESeqDataSetFromMatrix(countData = count_matrix_filtered,
                              colData = sample_sheet,
                              design = ~ 1) # No design needed

# Perform variance stabilizing transformation - excellent for dimensionality reduction
vsd <- vst(dds, blind = TRUE)
transformed_data <- t(assay(vsd)) # UMAP needs samples in rows, genes in columns

# Run UMAP
set.seed(123) # for reproducibility
# FIX: Set n_neighbors to be less than the number of samples to prevent error
n_neighbors_value <- min(15, nrow(transformed_data) - 1)
cat(paste("  Running UMAP with n_neighbors =", n_neighbors_value, "\n"))
umap_results <- umap(transformed_data, n_neighbors = n_neighbors_value)
colnames(umap_results) <- c("UMAP1", "UMAP2")

# Combine UMAP results with sample information
umap_plot_data <- as.data.frame(umap_results)
umap_plot_data$sample <- rownames(transformed_data)
umap_plot_data <- left_join(umap_plot_data, sample_sheet, by = "sample")

# ==============================================================================
# PART 3: Generate the UMAP Plot
# ==============================================================================
cat("--- PART 3: Generating the UMAP plot ---\n")

umap_plot <- ggplot(umap_plot_data, aes(x = UMAP1, y = UMAP2, color = condition, shape = condition)) +
  geom_point(size = 5, alpha = 0.8) +
  # Use ggrepel for smarter label placement
  geom_text_repel(aes(label = sample), show.legend = FALSE, box.padding = 0.5) +
  # Set custom colors
  scale_color_manual(values = c("Saliva" = "blue", "Blood" = "red", "Mixed" = "darkgreen")) +
  # Add titles and labels
  labs(
    title = "UMAP of Saliva, Blood, and Mixed Samples",
    subtitle = "Visualization of sample relationships based on gene expression",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    color = "Condition",
    shape = "Condition"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

# --- Save the Plot ---
plot_file <- file.path(output_dir, "sample_UMAP_plot.tiff")

ggsave(
  plot_file,
  plot = umap_plot,
  device = "tiff",
  width = 11,
  height = 9,
  units = "in",
  dpi = 600
)

cat(paste("  UMAP plot saved to:", plot_file, "\n"))
cat("\n--- UMAP analysis complete! ---\n")
