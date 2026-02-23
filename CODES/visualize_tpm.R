# ==============================================================================
# R Script to Visualize Top Expressed Genes from StringTie Output
# ==============================================================================

# --- 1. Load Necessary Libraries ---
# Make sure you have installed these packages first.
library(rtracklayer) # To read GTF files
library(dplyr)       # For data manipulation (filtering, arranging)
library(ggplot2)     # For creating the plots
library(stringr)     # For string manipulation

# --- 2. Configuration ---

# Directory where your StringTie GTF files are located
stringtie_dir <- "/home/alajdal/stringtie_quantification_results"

# List of your sample names
sample_names <- c(
  "SL3_rep",
  "SL4",
  "50SL50BL",
  "70SL30BL",
  "30SL70BL"
)

# A new directory to save the output plots
output_dir <- "top_50_known_genes_plots"

# Number of top genes to display in each plot
top_n <- 50

# --- 3. Create Output Directory ---
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# --- 4. Main Loop to Process Each Sample ---

for (sample in sample_names) {
  
  # Construct the full path to the input GTF file
  gtf_file_path <- file.path(stringtie_dir, paste0(sample, ".transcripts.gtf"))
  
  # Check if the file exists before proceeding
  if (!file.exists(gtf_file_path)) {
    warning(paste("GTF file not found for sample:", sample, "- Skipping."))
    next # Skip to the next sample
  }
  
  cat(paste("--- Processing sample:", sample, "---\n"))
  
  # --- Read and Process the GTF File ---
  
  # Import the GTF file
  gtf_data <- rtracklayer::import(gtf_file_path)
  
  # Convert to a data frame for easier manipulation
  gtf_df <- as.data.frame(gtf_data)
  
  # Check if the necessary columns exist
  if ("TPM" %in% names(gtf_df) && "ref_gene_name" %in% names(gtf_df)) {
    
    transcript_data <- gtf_df %>%
      filter(type == "transcript") %>%
      # UPDATED: This is the key change. We now filter to ONLY include transcripts
      # that have a known reference gene name, removing all "STRG" novel transcripts.
      filter(!is.na(ref_gene_name)) %>%
      select(ref_gene_name, TPM) %>% 
      mutate(TPM = as.numeric(TPM))
      
    # Aggregate TPM values by the official gene name
    gene_tpm <- transcript_data %>%
      group_by(ref_gene_name) %>%
      summarise(Total_TPM = sum(TPM)) %>%
      # Filter out mitochondrial and ribosomal RNA
      filter(!str_starts(ref_gene_name, "MT-"),
             !str_starts(ref_gene_name, "RN")) %>%
      arrange(desc(Total_TPM))
      
    # Get the top N genes AFTER all filtering
    top_genes <- head(gene_tpm, top_n)
    
    # --- Create the Bar Plot ---
    
    cat("  Generating plot...\n")
    
    top_genes_plot <- ggplot(top_genes, aes(x = reorder(ref_gene_name, Total_TPM), y = Total_TPM)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(
        title = paste("Top", top_n, "Known Expressed Genes in Sample:", sample),
        x = "Gene Name",
        y = "TPM (Transcripts Per Million)"
      ) +
      theme_bw(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5))
      
    # --- Save the Plot ---
    
    output_filename <- file.path(output_dir, paste0(sample, "_top_known_genes.tiff"))
    
    cat(paste("  Saving plot to:", output_filename, "\n"))
    
    ggsave(
      output_filename,
      plot = top_genes_plot,
      device = "tiff",
      width = 6.27,
      height = 9.69,
      units = "in",
      dpi = 600
    )
    
  } else {
    warning(paste("TPM or ref_gene_name attribute not found in GTF for sample:", sample))
  }
}

cat("\n--- All samples processed successfully! ---\n")
