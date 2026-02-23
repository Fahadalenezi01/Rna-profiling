# ==============================================================================
# R Script for Identifying Group-Specific Marker Genes
# and Heatmap Visualization
# ==============================================================================

# --- 1. Automated Package Installation ---
# This section checks for required packages and installs them if they are missing.
cat("--- Checking and installing necessary R packages ---\n")

# List of packages from CRAN
cran_packages <- c("dplyr", "pheatmap", "ggplot2", "tibble")
# List of packages from Bioconductor
bioc_packages <- c("rtracklayer")

# Install BiocManager if it's not already installed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Loop through CRAN packages and install if missing
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing CRAN package:", pkg, "\n"))
    install.packages(pkg)
  }
}

# Loop through Bioconductor packages and install if missing
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing Bioconductor package:", pkg, "\n"))
    BiocManager::install(pkg)
  }
}

cat("--- All packages are installed. Loading libraries. ---\n")


# --- 2. Load Necessary Libraries ---
library(rtracklayer)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(tibble)

# --- 3. Configuration & Sample Sheet ---
# This section defines the samples and their metadata.

# List your sample names
sample_names <- c(
  "SL1", "SL2", "SL3", "SL3_rep", "SL4",
  "BL1", "BL2", "BL3", "BL4"
)

# Define the condition (group) for each sample
conditions <- c(
  "Saliva", "Saliva", "Saliva", "Saliva", "Saliva",
  "Blood", "Blood", "Blood", "Blood"
)

# Create the main sample sheet data frame
sample_sheet <- data.frame(
  sample = sample_names,
  condition = factor(conditions) # 'factor' is important for annotations
)

# Define the full paths to the StringTie GTF files.
# Please double-check these paths to ensure they are correct.
base_dir_d <- "/mnt/d"
stringtie_dir_home <- "/home/alajdal/stringtie_quantification_results"

file_paths <- c(
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

# Add the file paths to our sample sheet
sample_sheet$filepath <- file_paths

# Define the output directory for this analysis
output_dir <- "group_specific_marker_analysis"

# --- 4. Create Output Directory ---
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# ==============================================================================
# PART 1: Import Data and Create TPM Matrix from GTF Files
# ==============================================================================
cat("--- PART 1: Importing StringTie quantification data to create TPM matrix ---\n")

# Initialize an empty list to store TPM data for each sample
all_tpms <- list()

# Loop through each file to extract TPM values
for (i in 1:nrow(sample_sheet)) {
  sample_name <- sample_sheet$sample[i]
  file_path <- sample_sheet$filepath[i]
  cat(paste("  Processing file:", basename(file_path), "\n"))

  if (!file.exists(file_path)) {
    warning(paste("File not found for sample", sample_name, "at path:", file_path, "- Skipping."))
    next
  }

  gtf_data <- tryCatch({
    rtracklayer::import(file_path)
  }, error = function(e) {
    warning(paste("Error reading GTF for sample", sample_name, ":", e$message, "- Skipping."))
    return(NULL)
  })

  if (is.null(gtf_data)) {
    next
  }
  
  gtf_df <- as.data.frame(gtf_data)

  if ("TPM" %in% names(gtf_df) && "gene_id" %in% names(gtf_df)) {
    gene_tpms <- gtf_df %>%
      filter(type == "transcript") %>%
      group_by(gene_id) %>%
      summarise(total_tpm = sum(as.numeric(TPM), na.rm = TRUE))
    
    colnames(gene_tpms)[colnames(gene_tpms) == "total_tpm"] <- sample_name
    all_tpms[[sample_name]] <- gene_tpms
  } else {
      warning(paste("Could not find 'TPM' or 'gene_id' columns in", basename(file_path)))
  }
}

tpm_matrix <- Reduce(function(x, y) full_join(x, y, by = "gene_id"), all_tpms)
tpm_matrix <- as.data.frame(tpm_matrix)
rownames(tpm_matrix) <- tpm_matrix$gene_id
tpm_matrix$gene_id <- NULL
tpm_matrix[is.na(tpm_matrix)] <- 0

processed_samples <- colnames(tpm_matrix)
sample_sheet <- sample_sheet[sample_sheet$sample %in% processed_samples, ]
tpm_matrix <- tpm_matrix[, sample_sheet$sample]
rownames(sample_sheet) <- sample_sheet$sample

gtf_map_data <- rtracklayer::import(sample_sheet$filepath[1])
gene_id_map <- as.data.frame(gtf_map_data) %>%
  filter(!is.na(ref_gene_name) & !is.na(gene_id)) %>%
  select(gene_id, ref_gene_name) %>%
  distinct()

# ==============================================================================
# PART 2: Filter Genes and Identify Group-Specific Markers
# ==============================================================================
cat("--- PART 2: Filtering genes and identifying group-specific markers ---\n")

log_tpm_matrix <- log2(tpm_matrix + 1)

# --- Gene Filtering Step ---
cat("  Filtering out housekeeping, mitochondrial, and ribosomal RNA genes...\n")

mt_pattern <- "^MT-"
rrna_pattern <- "^RN"
housekeeping_genes <- c("ACTB", "GAPDH", "B2M", "UBC", "HPRT1", "TBP", "RPL13A", "RPS18")

gene_names_df <- data.frame(gene_id = rownames(log_tpm_matrix)) %>%
  left_join(gene_id_map, by = "gene_id")

genes_to_remove <- gene_names_df$ref_gene_name[
  grepl(mt_pattern, gene_names_df$ref_gene_name, ignore.case = TRUE) |
  grepl(rrna_pattern, gene_names_df$ref_gene_name, ignore.case = TRUE) |
  gene_names_df$ref_gene_name %in% housekeeping_genes
]

gene_ids_to_remove <- gene_names_df$gene_id[gene_names_df$ref_gene_name %in% genes_to_remove]

genes_to_keep <- !(rownames(log_tpm_matrix) %in% gene_ids_to_remove)
original_gene_count <- nrow(log_tpm_matrix)
log_tpm_matrix_filtered <- log_tpm_matrix[genes_to_keep, ]
filtered_gene_count <- nrow(log_tpm_matrix_filtered)

cat(paste("  Removed", original_gene_count - filtered_gene_count, "genes.\n"))
cat(paste("  Continuing analysis with", filtered_gene_count, "genes.\n"))

# --- End of Filtering ---


# --- Robust Marker Selection using Wilcoxon Test ---
cat("  Performing Wilcoxon rank-sum test for each gene...\n")
saliva_samples <- sample_sheet$sample[sample_sheet$condition == "Saliva"]
blood_samples <- sample_sheet$sample[sample_sheet$condition == "Blood"]

# This function performs the test for a single gene (a row in the matrix)
run_wilcox_test <- function(gene_expression_row) {
  saliva_values <- as.numeric(gene_expression_row[saliva_samples])
  blood_values <- as.numeric(gene_expression_row[blood_samples])
  
  # Perform the test
  test_result <- wilcox.test(saliva_values, blood_values)
  
  # Calculate log2 fold change of medians for robustness
  log2_median_fc <- log2(median(saliva_values) + 1) - log2(median(blood_values) + 1)
  
  return(c(p_value = test_result$p.value, log2_median_fc = log2_median_fc))
}

# Apply the function to every gene (row) in the filtered matrix
test_results <- t(apply(log_tpm_matrix_filtered, 1, run_wilcox_test))
test_results_df <- as.data.frame(test_results)
test_results_df$gene_id <- rownames(test_results_df)

# Add readable gene names
test_results_df$gene_name <- gene_id_map$ref_gene_name[match(test_results_df$gene_id, gene_id_map$gene_id)]
test_results_df$gene_name[is.na(test_results_df$gene_name)] <- test_results_df$gene_id[is.na(test_results_df$gene_name)]


# Identify Saliva-specific markers (up-regulated in Saliva, significant p-value)
saliva_markers <- test_results_df %>%
  filter(log2_median_fc > 0) %>%
  arrange(p_value) %>%
  head(25)

# Identify Blood-specific markers (up-regulated in Blood, significant p-value)
blood_markers <- test_results_df %>%
  filter(log2_median_fc < 0) %>%
  arrange(p_value) %>%
  head(25)

# Combine the top markers for plotting
top_markers_df <- rbind(saliva_markers, blood_markers)

cat(paste("  Identified top 25 Saliva and 25 Blood specific markers based on statistical significance.\n"))

# Save the marker list to a CSV file
write.csv(
  top_markers_df,
  file = file.path(output_dir, "top_group_specific_markers_wilcox_filtered.csv"),
  row.names = FALSE
)
cat("  Full list of marker genes saved.\n")


# ==============================================================================
# PART 3: Generate the Heatmap of Top Marker Genes
# ==============================================================================
cat("--- PART 3: Generating heatmap of top group-specific markers ---\n")

# Order samples by condition to group them in the heatmap
sample_sheet <- sample_sheet[order(sample_sheet$condition), ]
ordered_tpm_matrix <- log_tpm_matrix_filtered[, sample_sheet$sample]

if (nrow(top_markers_df) > 0) {
    # Get the gene IDs and names for labeling, in the desired order
    top_gene_ids <- top_markers_df$gene_id
    top_gene_names <- top_markers_df$gene_name
    
    # Create the matrix for the plot using the ordered matrix and top genes
    plot_matrix <- ordered_tpm_matrix[top_gene_ids, ]
    rownames(plot_matrix) <- top_gene_names
    
    # Create annotation for rows to show which group they belong to
    row_annotation <- data.frame(
      Marker_Type = rep(c("Saliva Specific", "Blood Specific"), each = 25)
    )
    rownames(row_annotation) <- top_gene_names
    
    cat(paste("  Plotting heatmap for the top", nrow(plot_matrix), "marker genes.\n"))

    heatmap_file <- file.path(output_dir, "top_specific_markers_heatmap_wilcox.tiff")

    # Generate the heatmap
    pheatmap(
      plot_matrix,
      cluster_rows = FALSE, # Disable row clustering to keep Saliva/Blood markers grouped
      cluster_cols = FALSE, # Disable column clustering to keep sample order
      show_rownames = TRUE,
      show_colnames = TRUE,
      annotation_col = sample_sheet[, "condition", drop = FALSE],
      annotation_row = row_annotation,
      main = "Top 25 Saliva vs. Blood Specific Markers\n(log2(TPM+1), Filtered, Wilcoxon Test)",
      fontsize_row = 8,
      scale = "row", # Scale rows to see pattern instead of absolute expression
      filename = heatmap_file,
      width = 8,
      height = 10,
      dpi = 300
    )

    cat(paste("  Heatmap saved to:", heatmap_file, "\n"))
} else {
    cat("\n!!! NOTE: No marker genes were found, so no heatmap can be generated. !!!\n")
}

cat("\n--- Marker identification analysis complete! ---\n")
