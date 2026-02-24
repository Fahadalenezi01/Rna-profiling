# ==============================================================================
# R Script for Principal Component Analysis (PCA) of Saliva vs. Blood Samples
# ==============================================================================

# --- 1. Automated Package Installation ---
cat("--- Checking and installing necessary R packages ---\n")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
packages <- c("rtracklayer", "dplyr", "tidyr", "ggplot2", "DESeq2", "org.Hs.eg.db")
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
library(ggplot2)
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)

# --- 3. Configuration & Sample Sheet ---
stringtie_dir_home <- "/home/alajdal/stringtie_quantification_results"
base_dir_d <- "/mnt/d"
output_dir <- "pca_analysis_results"

# UPDATED: The sample sheet now includes the "Mixed" samples and removes "SL3".
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

# Convert to a count matrix (DESeq2 works best with counts)
# We will just round the TPM values for this visualization purpose
count_matrix <- round(tpm_matrix)

# Filter for expressed genes: Keep genes with count > 10 in at least 2 samples
keep <- rowSums(count_matrix >= 10) >= 2
count_matrix_filtered <- count_matrix[keep, ]

# ==============================================================================
# PART 2: Perform PCA
# ==============================================================================
cat("--- PART 2: Performing Principal Component Analysis ---\n")

# Create a DESeq2 object (we won't run the full analysis, just use it for transformation)
dds <- DESeqDataSetFromMatrix(countData = count_matrix_filtered,
                              colData = sample_sheet,
                              design = ~ 1) # No design needed for PCA

# Perform a variance stabilizing transformation - the best input for PCA
vsd <- vst(dds, blind = TRUE)

# ==============================================================================
# PART 3: Generate the PCA Plot
# ==============================================================================
cat("--- PART 3: Generating the PCA plot ---\n")

# Use the built-in DESeq2 function to generate the plot data
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

# Calculate the percentage of variance explained by each principal component
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Create the final, publication-quality plot using ggplot2
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  geom_point(size = 5) +
  # Add labels to the points
  geom_text(aes(label = name), vjust = -1, show.legend = FALSE) +
  # Set custom colors
  scale_color_manual(values = c("Saliva" = "blue", "Blood" = "red", "Mixed" = "darkgreen")) +
  # Add titles and axis labels with variance explained
  labs(
    title = "PCA of Saliva, Blood, and Mixed Samples",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "Condition",
    shape = "Condition"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    plot.margin = unit(c(1, 1, 2, 1), "cm")
  )

# --- Save the Plot ---
pca_plot_file <- file.path(output_dir, "saliva_vs_blood_vs_mixed_PCA_plot.tiff")

ggsave(
  pca_plot_file,
  plot = pca_plot,
  device = "tiff",
  width = 8.5,
  height = 8.5, # Make it square for better proportions
  units = "in",
  dpi = 600
)

cat(paste("  PCA plot saved to:", pca_plot_file, "\n"))
cat("\n--- PCA analysis complete! ---\n")
