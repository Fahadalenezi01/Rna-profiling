# ==============================================================================
# R Script to Generate Thesis-Quality Log-Scale TPM Profile Plots
# ==============================================================================

# --- 1. Automated Package Installation ---
# This section ensures all necessary packages, including the human gene
# annotation database, are installed.
cat("--- Checking and installing necessary R packages ---\n")

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
packages <- c("rtracklayer", "dplyr", "ggplot2", "stringr", "org.Hs.eg.db")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(paste("Installing Bioconductor package:", pkg, "\n"))
    BiocManager::install(pkg)
  }
}
cat("--- All packages are installed. Loading libraries. ---\n")


# --- 2. Load Necessary Libraries ---
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(stringr)
library(AnnotationDbi) # Added for clarity
library(org.Hs.eg.db) # The human gene annotation database

# --- 3. Configuration ---

# Directory where your StringTie GTF files are located
stringtie_dir_home <- "/home/alajdal/stringtie_quantification_results"
base_dir_d <- "/mnt/d"

# A new, separate directory to save the final thesis-quality plots
output_dir <- "log_tpm_profile_plots_thesis_quality"

# Number of top genes to display in each plot
top_n <- 50

# --- 4. Sample Sheet Configuration ---
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

# --- 5. Create Output Directory ---
if (!dir.exists(output_dir)) dir.create(output_dir)

# ==============================================================================
# Main Loop to Process Each Sample
# ==============================================================================

for (i in 1:nrow(sample_sheet)) {
  
  sample_name <- sample_sheet$sample[i]
  sample_condition <- sample_sheet$condition[i]
  gtf_file_path <- sample_sheet$filepath[i]
  
  if (!file.exists(gtf_file_path)) {
    warning(paste("GTF file not found for sample:", sample_name, "- Skipping."))
    next
  }
  
  cat(paste("--- Processing sample:", sample_name, "---\n"))
  
  # --- Read and Process the GTF File ---
  gtf_df <- as.data.frame(rtracklayer::import(gtf_file_path))
  
  if ("TPM" %in% names(gtf_df) && "ref_gene_name" %in% names(gtf_df)) {
    
    gene_tpm <- gtf_df %>%
      filter(type == "transcript", !is.na(ref_gene_name)) %>%
      filter(seqnames != "chrM") %>%
      # CORRECTED: Explicitly use dplyr's select function to avoid conflict
      dplyr::select(ref_gene_name, TPM) %>% 
      mutate(TPM = as.numeric(TPM)) %>%
      group_by(ref_gene_name) %>%
      summarise(Total_TPM = sum(TPM)) %>%
      filter(!str_starts(ref_gene_name, "MT-"), !str_starts(ref_gene_name, "RN")) %>%
      mutate(Log10_TPM = log10(Total_TPM + 1)) %>%
      arrange(desc(Total_TPM))
      
    # --- NEW: Convert Ensembl IDs to Gene Symbols ---
    # Strip version numbers from Ensembl IDs (e.g., ENSG000.1 -> ENSG000)
    gene_tpm$ensembl_id_clean <- gsub("\\..*$", "", gene_tpm$ref_gene_name)
    
    # Use the annotation database to map IDs to symbols
    gene_symbols <- mapIds(org.Hs.eg.db,
                           keys = gene_tpm$ensembl_id_clean,
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")
                           
    # Create the final display name: use the symbol if found, otherwise use the original ID
    gene_tpm$display_name <- ifelse(!is.na(gene_symbols), gene_symbols, gene_tpm$ref_gene_name)
    
    # Get the top N genes
    top_genes <- head(gene_tpm, top_n)
    
    # --- Create the Bar Plot ---
    cat("  Generating plot...\n")
    
    plot_color <- ifelse(sample_condition == "Saliva", "blue", "red")
    plot_title <- str_wrap(paste("Top 50 Expressed Genes in Sample:", sample_name), width = 60)
    
    log_tpm_plot <- ggplot(top_genes, aes(x = reorder(display_name, Total_TPM), y = Log10_TPM)) +
      geom_bar(stat = "identity", fill = plot_color) +
      coord_flip() +
      labs(
        title = plot_title,
        subtitle = "Expression shown on a Log10 scale",
        x = "Gene",
        y = "Log10 (TPM + 1)"
      ) +
      # UPDATED: Using a cleaner theme and smaller text for better thesis presentation
      theme_classic(base_size = 10) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 8), # Smaller font for gene names
        plot.margin = unit(c(1, 1, 2, 1), "cm") # Added bottom margin for caption
      )
      
    # --- Save the Plot ---
    output_filename <- file.path(output_dir, paste0(sample_name, "_log_tpm_profile_thesis.tiff"))
    cat(paste("  Saving plot to:", output_filename, "\n"))
    
    ggsave(
      output_filename,
      plot = log_tpm_plot,
      device = "tiff",
      width = 6.27, # A4 width with 1-inch margins
      height = 8.5, # Reduced height for better proportion and caption space
      units = "in",
      dpi = 600
    )
    
  } else {
    warning(paste("TPM or ref_gene_name attribute not found in GTF for sample:", sample_name))
  }
}

cat("\n--- All samples processed successfully! ---\n")
