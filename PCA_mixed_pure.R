# ==============================================================================
# R Script for Targeted PCA using Saliva and Blood Marker Genes
# ==============================================================================

# --- 1. Automated Package Installation ---
# This section checks for required packages and installs them if they are missing.
cat("--- Checking and installing necessary R packages ---\n")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# List of required packages
packages <- c("rtracklayer", "dplyr", "tidyr", "ggplot2", "ggrepel", "DESeq2")

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
library(DESeq2)

# --- 3. Configuration & Sample Sheet ---
stringtie_dir_home <- "/home/alajdal/stringtie_quantification_results"
base_dir_d <- "/mnt/d"
output_dir <- "targeted_pca_analysis"

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

# ==============================================================================
# PART 2: Identify Robust Saliva and Blood Marker Genes
# ==============================================================================
cat("--- PART 2: Identifying robust marker genes from pure samples ---\n")

log_tpm_matrix <- log2(tpm_matrix + 1)

saliva_samples_pure <- sample_sheet$sample[sample_sheet$condition == "Saliva"]
blood_samples_pure <- sample_sheet$sample[sample_sheet$condition == "Blood"]

run_wilcox_test <- function(gene_expression_row) {
  saliva_values <- as.numeric(gene_expression_row[saliva_samples_pure])
  blood_values <- as.numeric(gene_expression_row[blood_samples_pure])
  if (var(saliva_values) == 0 && var(blood_values) == 0) return(1)
  test_result <- wilcox.test(saliva_values, blood_values)
  return(test_result$p.value)
}

p_values <- apply(log_tpm_matrix, 1, run_wilcox_test)
log2fc <- rowMeans(log_tpm_matrix[, saliva_samples_pure]) - rowMeans(log_tpm_matrix[, blood_samples_pure])

results_df <- data.frame(gene = names(p_values), p_value = p_values, log2fc = log2fc)

saliva_markers <- results_df %>% filter(log2fc > 1) %>% arrange(p_value) %>% head(25)
blood_markers <- results_df %>% filter(log2fc < -1) %>% arrange(p_value) %>% head(25)
top_50_markers <- c(saliva_markers$gene, blood_markers$gene)

cat(paste("  Identified", length(top_50_markers), "total marker genes for PCA.\n"))

# ==============================================================================
# PART 3: Perform PCA using only the Marker Genes
# ==============================================================================
cat("--- PART 3: Performing PCA using only marker genes ---\n")

# Convert TPMs to integer counts for DESeq2
count_matrix <- round(tpm_matrix)

# Filter the count matrix to include only the top 50 markers
marker_count_matrix <- count_matrix[rownames(count_matrix) %in% top_50_markers, ]

# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = marker_count_matrix,
                              colData = sample_sheet,
                              design = ~ 1)

# Perform variance stabilizing transformation
# FIX: Use varianceStabilizingTransformation directly when number of genes is small
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# Generate PCA data using the built-in DESeq2 function
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# ==============================================================================
# PART 4: Generate the Targeted PCA Plot
# ==============================================================================
cat("--- PART 4: Generating the targeted PCA plot ---\n")

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  geom_point(size = 5, alpha = 0.8) +
  geom_text_repel(aes(label = name), show.legend = FALSE, box.padding = 0.5) +
  scale_color_manual(values = c("Saliva" = "blue", "Blood" = "red", "Mixed" = "darkgreen")) +
  labs(
    title = "Targeted PCA of Saliva, Blood, and Mixed Samples",
    subtitle = "Analysis based on top 50 Saliva & Blood marker genes",
    x = paste0("PC1: ", percentVar[1], "% variance (Tissue Identity)"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
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
plot_file <- file.path(output_dir, "targeted_marker_PCA_plot.tiff")

ggsave(
  plot_file,
  plot = pca_plot,
  device = "tiff",
  width = 11,
  height = 9,
  units = "in",
  dpi = 600
)

cat(paste("  Targeted PCA plot saved to:", plot_file, "\n"))
cat("\n--- Targeted PCA analysis complete! ---\n")
