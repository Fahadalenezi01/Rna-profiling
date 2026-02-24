# ==============================================================================
# R Script for Mixture Analysis using Marker Gene Signatures
# ==============================================================================

# --- 1. Automated Package Installation ---
# This section checks for required packages and installs them if they are missing.
cat("--- Checking and installing necessary R packages ---\n")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# List of required packages
packages <- c("rtracklayer", "dplyr", "tidyr", "ggplot2", "tibble")

# Loop through packages and install if missing
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing package:", pkg, "\n"))
    if (pkg == "rtracklayer") {
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
library(tibble)

# --- 3. Configuration & Sample Sheet ---
stringtie_dir_home <- "/home/alajdal/stringtie_quantification_results"
base_dir_d <- "/mnt/d"
output_dir <- "mixture_analysis_results"

# Sample sheet now includes ALL samples to identify markers and then analyze mixtures
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

# Get pure sample names for marker identification
saliva_samples_pure <- sample_sheet$sample[sample_sheet$condition == "Saliva"]
blood_samples_pure <- sample_sheet$sample[sample_sheet$condition == "Blood"]

# Perform Wilcoxon test to find differentially expressed genes between pure samples
run_wilcox_test <- function(gene_expression_row) {
  saliva_values <- as.numeric(gene_expression_row[saliva_samples_pure])
  blood_values <- as.numeric(gene_expression_row[blood_samples_pure])
  # Add error handling for cases with no variance
  if (var(saliva_values) == 0 && var(blood_values) == 0) return(1)
  test_result <- wilcox.test(saliva_values, blood_values)
  return(test_result$p.value)
}

p_values <- apply(log_tpm_matrix, 1, run_wilcox_test)

# Calculate log2 fold change between pure samples
log2fc <- rowMeans(log_tpm_matrix[, saliva_samples_pure]) - rowMeans(log_tpm_matrix[, blood_samples_pure])

# Create a results data frame
results_df <- data.frame(
  gene = names(p_values),
  p_value = p_values,
  log2fc = log2fc
)

# Select top 25 Saliva markers and top 25 Blood markers
saliva_markers <- results_df %>% filter(log2fc > 0) %>% arrange(p_value) %>% head(25)
blood_markers <- results_df %>% filter(log2fc < 0) %>% arrange(p_value) %>% head(25)

cat(paste("  Identified", nrow(saliva_markers), "saliva markers and", nrow(blood_markers), "blood markers.\n"))

# ==============================================================================
# PART 3: Calculate and Normalize Signature Scores
# ==============================================================================
cat("--- PART 3: Calculating and normalizing signature scores ---\n")

# --- Step 3a: Calculate reference scores from PURE samples ---
pure_saliva_matrix <- log_tpm_matrix[ , saliva_samples_pure]
pure_blood_matrix <- log_tpm_matrix[ , blood_samples_pure]

# The reference score is the average expression of the signature in its own tissue type
saliva_ref_score <- mean(as.matrix(pure_saliva_matrix[saliva_markers$gene, ]))
blood_ref_score <- mean(as.matrix(pure_blood_matrix[blood_markers$gene, ]))

cat(paste("  Saliva Reference Score:", round(saliva_ref_score, 2), "\n"))
cat(paste("  Blood Reference Score:", round(blood_ref_score, 2), "\n"))

# --- Step 3b: Calculate scores for MIXED samples and normalize them ---
mixed_samples <- sample_sheet$sample[sample_sheet$condition == "Mixed"]
mixed_matrix <- log_tpm_matrix[ , mixed_samples]

# Calculate raw scores for mixed samples
saliva_score_raw <- colMeans(mixed_matrix[saliva_markers$gene, ])
blood_score_raw <- colMeans(mixed_matrix[blood_markers$gene, ])

# Normalize the raw scores by the reference scores
saliva_score_norm <- saliva_score_raw / saliva_ref_score
blood_score_norm <- blood_score_raw / blood_ref_score

# Combine scores into a data frame
scores_df <- data.frame(
  sample = mixed_samples,
  Saliva_Score = saliva_score_norm,
  Blood_Score = blood_score_norm
)

# Reshape data for ggplot2
plot_data <- scores_df %>%
  pivot_longer(cols = c("Saliva_Score", "Blood_Score"),
               names_to = "Marker_Type",
               values_to = "Score")

# Order the samples logically for the plot
plot_data$sample <- factor(plot_data$sample, levels = c("70SL30BL", "50SL50BL", "30SL70BL"))

# ==============================================================================
# PART 4: Generate 100% Stacked Bar Plot
# ==============================================================================
cat("--- PART 4: Generating the proportional mixture analysis plot ---\n")

mixture_plot <- ggplot(plot_data, aes(x = sample, y = Score, fill = Marker_Type)) +
  # Use position = "fill" to create a 100% stacked bar plot from the NORMALIZED scores
  geom_bar(stat = "identity", position = "fill") +
  
  # Customize colors and labels
  scale_fill_manual(values = c("Saliva_Score" = "blue", "Blood_Score" = "red"),
                    name = "Marker Type") +
  # Change y-axis to represent percentage
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proportional Analysis of Mixed Samples",
    subtitle = "Relative contribution of Saliva and Blood marker signatures (Normalized)",
    x = "Mixed Sample",
    y = "Proportion of Marker Signature"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# --- Save the Plot ---
plot_file <- file.path(output_dir, "mixture_analysis_normalized_proportional_plot.tiff")

ggsave(
  plot_file,
  plot = mixture_plot,
  device = "tiff",
  width = 10,
  height = 8,
  units = "in",
  dpi = 600
)

cat(paste("  Normalized proportional mixture plot saved to:", plot_file, "\n"))
cat("\n--- Analysis complete! ---\n")
