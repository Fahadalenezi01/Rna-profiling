#!/usr/bin/env Rscript

# One-vs-Rest TPM Limma Pipeline (Full DE Output)
# ===============================================
# 1) Create GeneID -> GeneName mapping
# 2) Merge StringTie TPM by GeneID
# 3) Fit limma one-vs-rest contrasts on log2(TPM+1)
# 4) Extract and save FULL Differential Expression results (NEW)
# 5) Extract DE results for marker identification
# 6) Apply TPM specificity filter
# 7) Enforce exclusivity, map IDs to Names
# 8) Write detailed marker tables (TSV) and simple lists
# 9) Generate heatmap with Gene Names

# 0. Install/load dependencies
# ---------------------------
cran_pkgs <- c("dplyr", "purrr", "tibble", "matrixStats", "pheatmap", "limma")
for(pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(matrixStats)
  library(pheatmap)
  library(limma)
})

# --- Configuration ---
# Define TPM file paths
tpm_paths <- list(
  BL1 = "/home/datavis/Desktop/BL1/filtered_stringtie_abundance_BL1.txt",
  BL2 = "/home/datavis/Desktop/BL2/filtered_stringtie_abundance_BL2.txt",
  BL3 = "/home/datavis/Desktop/BL3/filtered_stringtie_abundance_BL3.txt",
  BL4 = "/home/datavis/Desktop/BL4/filtered_stringtie_abundance_BL4.txt",
  SL1 = "/home/datavis/Desktop/SL1/filtered_stringtie_abundance_SL1.txt",
  SL2 = "/home/datavis/Desktop/SL2/filtered_stringtie_abundance_SL2.txt",
  SL3 = "/home/datavis/Desktop/SL3/filtered_stringtie_abundance_SL3.txt"
)
# Output directory
output_dir <- "/home/datavis/Desktop"

# Filtering Thresholds (ADJUST THESE TO GET MORE/FEWER GENES)
# ---------------------------------------------------------
# For initial DE list and marker filtering
lfc_cut <- 0.15        # Minimum log2 Fold Change
fdr_cut <- 0.25        # Maximum Adjusted P-value (FDR)
# For TPM specificity filter (applied *after* initial DE filter)
min_tpm_target <- 1.0  # Minimum mean TPM in the target group
max_tpm_other  <- 1.5  # Maximum mean TPM allowed in the other group
# ---------------------------------------------------------

# Heatmap settings
heatmap_max_genes <- 50
# --- End Configuration ---


# 1. Create Gene ID -> Gene Name Mapping
# --------------------------------------
extract_gene_info <- function(file_path) {
  df <- tryCatch({
      read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE, nrows = 1e5)
  }, error = function(e) {
      message("Reading ", basename(file_path), " as tab-delimited failed, trying comma.")
      tryCatch({
          read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE, nrows = 1e5)
      }, error = function(e2) {
          stop("Failed to read file: ", file_path, ". Error: ", e2$message)
      })
  })

  gid_col <- grep("gene[_ ]?id", names(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(gid_col)) stop("No Gene ID column found in: ", file_path)
  gname_col <- grep("^(gene[ _]?name|symbol)$", names(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(gname_col)) stop("No Gene Name/Symbol column found in: ", file_path, ". Needed for mapping.")

  df %>%
    select(GeneID = all_of(gid_col), GeneName = all_of(gname_col)) %>%
    filter(!is.na(GeneID) & GeneID != "" & !is.na(GeneName) & GeneName != "") %>%
    distinct(GeneID, .keep_all = TRUE)
}

all_gene_info <- map_dfr(tpm_paths, extract_gene_info)
gene_map <- all_gene_info %>% distinct(GeneID, .keep_all = TRUE)

if (nrow(gene_map) == 0) {
    stop("Gene ID to Gene Name mapping table is empty. Check input files and column names.")
}
message("Created Gene ID -> Gene Name mapping table with ", nrow(gene_map), " unique entries.")


# 2. Load and collapse TPM by GeneID
# ---------------------------------
tpm_list <- map2(names(tpm_paths), tpm_paths, function(sample_id, file_path) {
  df <- tryCatch({
      read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)
  }, error = function(e) {
      tryCatch({
          read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
      }, error = function(e2) {
          stop("Failed to read file: ", file_path, ". Error: ", e2$message)
      })
  })

  gid_col <- grep("gene[_ ]?id", names(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(gid_col)) stop("Could not find a 'Gene ID' column in file: ", file_path)
  df <- df %>% rename(GeneID = all_of(gid_col))
  if (!"TPM" %in% names(df)) stop("Could not find a 'TPM' column in file: ", file_path)

  df %>%
    group_by(GeneID) %>%
    summarise(!!sample_id := sum(as.numeric(TPM), na.rm = TRUE), .groups = 'drop') %>%
    filter(!is.na(GeneID) & GeneID != "")
})

# 3. Merge all samples and prepare matrix
# ---------------------------------------
tpm_df <- reduce(tpm_list, full_join, by = "GeneID")
tpm_df[is.na(tpm_df)] <- 0
tpm_df <- as.data.frame(tpm_df)
rownames(tpm_df) <- tpm_df$GeneID
tpm_mat <- as.matrix(tpm_df[, -1])
log2_tpm <- log2(tpm_mat + 1)

# 4. Define sample groups and design matrix for limma
# --------------------------------------------------
samples <- colnames(log2_tpm)
groups <- ifelse(startsWith(samples, "BL"), "Blood", "Saliva")
design <- model.matrix(~ 0 + factor(groups))
colnames(design) <- levels(factor(groups))

# 5. Fit limma model and define contrasts
# ---------------------------------------
tfit <- lmFit(log2_tpm, design)
cont <- makeContrasts(
  Blood_vs_Rest  = Blood - Saliva, # Contrast of interest for full DE table
  Saliva_vs_Rest = Saliva - Blood, # Needed for Saliva markers later
  levels = design
)
fit2 <- contrasts.fit(tfit, cont) %>% eBayes()

# 6. Extract and Save FULL Differential Expression Results (Blood vs Saliva)
# -------------------------------------------------------------------------
# Get results for ALL genes tested for the Blood vs Saliva contrast
# Use n=Inf to get all genes, sort by p-value by default
full_de_results <- topTable(fit2, coef = "Blood_vs_Rest", n = Inf, sort.by = "P", adjust.method = "BH") %>%
  rownames_to_column("GeneID") %>%
  left_join(gene_map, by = "GeneID") %>% # Add Gene Names
  select(GeneID, GeneName, logFC, AveExpr, P.Value, adj.P.Val) %>% # Select and reorder columns
  # Optionally filter here if you only want significant genes in this table:
  # filter(abs(logFC) > lfc_cut, adj.P.Val < fdr_cut) %>%
  arrange(adj.P.Val, P.Value) # Sort by adjusted p-value, then p-value

# Write the full DE table to a TSV file
full_de_filename <- file.path(output_dir, "limma_full_DE_results_Blood_vs_Saliva.tsv")
write.table(full_de_results,
            full_de_filename,
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
message("Full Differential Expression results (Blood vs Saliva) saved to: ", full_de_filename)


# --- Now proceed with Marker Identification using the defined thresholds ---

# 7. Extract DE results passing initial LFC/FDR thresholds for marker finding
# ---------------------------------------------------------------------------
# Note: We re-run topTable here, applying the lfc/fdr cutoffs,
#       or we could filter the 'full_de_results' table. Re-running is clear.
res_blood_marker_candidates  <- topTable(fit2, coef = "Blood_vs_Rest",  n = Inf, adjust.method = "BH") %>%
  rownames_to_column("GeneID") %>%
  select(GeneID, logFC, AveExpr, P.Value, adj.P.Val) %>%
  filter(logFC > lfc_cut, adj.P.Val < fdr_cut) # Filter for UP in Blood

res_saliva_marker_candidates <- topTable(fit2, coef = "Saliva_vs_Rest", n = Inf, adjust.method = "BH") %>%
  rownames_to_column("GeneID") %>%
  select(GeneID, logFC, AveExpr, P.Value, adj.P.Val) %>%
  filter(logFC > lfc_cut, adj.P.Val < fdr_cut) # Filter for UP in Saliva


# 8. TPM specificity filter (using Gene IDs from marker candidates)
# -----------------------------------------------------------------
mean_tpm_safe <- function(gene_id, target_group, tpm_matrix, group_vector) {
  if (!gene_id %in% rownames(tpm_matrix)) return(NA)
  group_samples <- group_vector == target_group
  mean(tpm_matrix[gene_id, group_samples], na.rm = TRUE)
}

ids_passing_tpm_blood <- res_blood_marker_candidates$GeneID %>%
  keep(~ {
    mean_b <- mean_tpm_safe(., "Blood",  tpm_mat, groups)
    mean_s <- mean_tpm_safe(., "Saliva", tpm_mat, groups)
    !is.na(mean_b) && !is.na(mean_s) && mean_b >= min_tpm_target && mean_s < max_tpm_other
  })

ids_passing_tpm_saliva <- res_saliva_marker_candidates$GeneID %>%
  keep(~ {
    mean_s <- mean_tpm_safe(., "Saliva", tpm_mat, groups)
    mean_b <- mean_tpm_safe(., "Blood",  tpm_mat, groups)
    !is.na(mean_s) && !is.na(mean_b) && mean_s >= min_tpm_target && mean_b < max_tpm_other
  })

# Subset the candidate results tables to keep only genes passing TPM filter
res_blood_filtered  <- res_blood_marker_candidates %>% filter(GeneID %in% ids_passing_tpm_blood)
res_saliva_filtered <- res_saliva_marker_candidates %>% filter(GeneID %in% ids_passing_tpm_saliva)

# 9. Remove mitochondrial genes and enforce exclusivity for markers
# -----------------------------------------------------------------
markers_blood_ids  <- res_blood_filtered$GeneID %>% setdiff(grep("^MT-", ., value = TRUE))
markers_saliva_ids <- res_saliva_filtered$GeneID %>% setdiff(grep("^MT-", ., value = TRUE))
common_markers_ids <- intersect(markers_blood_ids, markers_saliva_ids)
markers_blood_ids_final  <- setdiff(markers_blood_ids,  common_markers_ids)
markers_saliva_ids_final <- setdiff(markers_saliva_ids, common_markers_ids)

# 10. Create final marker tables with stats and names
# ---------------------------------------------------
final_markers_blood_table <- res_blood_filtered %>%
  filter(GeneID %in% markers_blood_ids_final) %>%
  left_join(gene_map, by = "GeneID") %>%
  select(GeneID, GeneName, logFC, adj.P.Val, P.Value, AveExpr) %>%
  arrange(desc(logFC))

final_markers_saliva_table <- res_saliva_filtered %>%
  filter(GeneID %in% markers_saliva_ids_final) %>%
  left_join(gene_map, by = "GeneID") %>%
  select(GeneID, GeneName, logFC, adj.P.Val, P.Value, AveExpr) %>%
  arrange(desc(logFC))

# 11. Write marker tables (TSV for Word) and simple lists
# -------------------------------------------------------
# Write detailed MARKER tables
write.table(final_markers_blood_table,
            file.path(output_dir, "blood_vs_rest_markers_DETAILED.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(final_markers_saliva_table,
            file.path(output_dir, "saliva_vs_rest_markers_DETAILED.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
message("Detailed MARKER tables (TSV) saved to: ", output_dir)

# Write simple MARKER Gene Name lists
write.table(final_markers_blood_table$GeneName,
            file.path(output_dir, "blood_vs_rest_marker_NAMES.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(final_markers_saliva_table$GeneName,
            file.path(output_dir, "saliva_vs_rest_marker_NAMES.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
message("Simple marker gene name lists saved to: ", output_dir)


# 12. Generate Heatmap with Top Markers (Gene Names)
# --------------------------------------------------
all_marker_ids_final <- c(markers_blood_ids_final, markers_saliva_ids_final)
combined_final_table <- bind_rows(
    final_markers_blood_table %>% mutate(Group = "Blood"),
    final_markers_saliva_table %>% mutate(Group = "Saliva")
) %>% arrange(desc(logFC))

hm_gene_ids <- head(combined_final_table$GeneID, heatmap_max_genes)

if (length(hm_gene_ids) >= 2) {
  hm_data <- log2_tpm[hm_gene_ids, , drop = FALSE]
  hm_gene_names <- gene_map %>%
                     filter(GeneID %in% hm_gene_ids) %>%
                     arrange(match(GeneID, hm_gene_ids)) %>%
                     pull(GeneName)

   # Fallback for missing names
   if(length(hm_gene_names) != length(hm_gene_ids)) {
       warning("Mismatch between number of gene IDs and mapped names for heatmap. Using IDs.")
       rownames(hm_data) <- hm_gene_ids
   } else if (any(is.na(hm_gene_names) | hm_gene_names == "")) {
       warning("Some gene names are missing/empty for heatmap rows. Using IDs for those.")
       missing_indices <- which(is.na(hm_gene_names) | hm_gene_names == "")
       hm_gene_names[missing_indices] <- hm_gene_ids[missing_indices]
       rownames(hm_data) <- hm_gene_names
   } else {
        rownames(hm_data) <- hm_gene_names
   }

  group_order <- order(groups)
  hm_data <- hm_data[, group_order, drop = FALSE]
  ordered_groups <- groups[group_order]
  ordered_samples <- samples[group_order]

  anno_col <- data.frame(Group = factor(ordered_groups, levels = c("Saliva", "Blood")),
                         row.names = ordered_samples)
  anno_colors <- list(Group = c(Saliva = "#4575b4", Blood = "#d73027"))

  heatmap_file <- file.path(output_dir, "one_vs_rest_TPM_limma_heatmap_GeneNames.png")
  png_height = max(1000, length(hm_gene_names) * 25)
  png(heatmap_file, width = 1200, height = png_height, res = 100)
  pheatmap(
    hm_data,
    scale             = "row",
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = 8,
    fontsize_col      = 10,
    annotation_col    = anno_col,
    annotation_colors = anno_colors,
    main              = paste("Top", length(hm_gene_names), "Specific Markers (One-vs-Rest TPM Limma)")
  )
  dev.off()
  message("Heatmap with Top Marker Gene Names saved to: ", heatmap_file)

} else {
  message("Fewer than 2 specific markers found after filtering. Heatmap generation skipped.")
}

message("Script finished successfully.")

