# Load required libraries
library(ggplot2)
library(dplyr)
merged_counts <- read.delim("/home/datavis/Desktop/BL4/filtered_stringtie_abundance_BL4.txt", header=TRUE)
# Ensure TPM is numeric
merged_counts$TPM <- as.numeric(merged_counts$TPM)

# Filter to keep rows that have a non-empty Gene.Name
filtered_all <- merged_counts[!is.na(merged_counts$Gene.Name) & merged_counts$Gene.Name != "", ]

# Sort the data by TPM in descending order
sorted_all <- filtered_all[order(filtered_all$TPM, decreasing = TRUE), ]

# Extract the top 20 genes
top20 <- head(sorted_all, 20)

# Convert 'Gene.Name' to a factor, ordered by TPM (lowest at bottom for a horizontal bar chart)
top20$Gene.Name <- factor(top20$Gene.Name, levels = unique(top20$Gene.Name[order(top20$TPM)]))


# Create a horizontal bar plot using ggplot2
ggplot(top20, aes(x = Gene.Name, y = TPM)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 Highly Expressed Genes (TPM)",
       x = "Gene",
       y = "TPM") +
  theme_minimal()
  
  ggsave("/home/datavis/Desktop/BL4/filtered_Top20Genes_BL4.png", width = 10, height = 8, dpi = 300)
  
  
  
  
  
  ###after filtering
  library(ggplot2)
library(dplyr)
library(scales)  # for comma formatting of numbers

# Read in the file and ensure TPM is numeric
merged_counts <- read.delim("/home/datavis/Desktop/SL3/filtered_stringtie_abundance_SL3.txt", header = TRUE)
# If necessary, clean the TPM column (e.g., remove commas) and convert it
# merged_counts$TPM <- as.numeric(gsub(",", "", merged_counts$TPM))
merged_counts$TPM <- as.numeric(merged_counts$TPM)

# Filter rows with a valid Gene.Name
filtered_all <- merged_counts %>%
  filter(!is.na(Gene.Name) & Gene.Name != "")

# Aggregate TPM values by gene to combine duplicated transcripts
aggregated <- filtered_all %>%
  group_by(Gene.Name) %>%
  summarise(TPM = sum(TPM, na.rm = TRUE)) %>%
  ungroup()

# Optionally, if mitochondrial genes are to be removed, uncomment the following:
# aggregated <- aggregated %>% filter(!grepl("^MT", Gene.Name, ignore.case = TRUE))

# For a blood sample it may be useful to keep hemoglobin genes so we leave them in.
# Sort and select the top 20 genes
top20 <- aggregated %>%
  arrange(desc(TPM)) %>%
  head(20)

# Order 'Gene.Name' factor for proper ordering in the horizontal bar plot
top20$Gene.Name <- factor(top20$Gene.Name, levels = unique(top20$Gene.Name[order(top20$TPM)]))

# Create a horizontal bar plot using a log scale and annotate the bars with TPM values
p <- ggplot(top20, aes(x = Gene.Name, y = TPM)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  scale_y_log10(labels = comma) +  # Log scale compresses large values and uses comma formatting
  coord_flip() +
  labs(x = "Gene Name", y = "TPM (log scale)", title = "Top 20 Expressed Genes in Saliva") +
  theme_minimal() +
  geom_text(aes(label = comma(round(TPM, 0))), 
            hjust = -0.1, 
            size = 3)

# Save the plot in the same directory
ggsave("/home/datavis/Desktop/SL3/filtered_Top20Genes_SL3.png", plot = p, width = 10, height = 8, dpi = 300)




#####heat map
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(plotly)

# Define file paths
blood_files <- c("/home/datavis/Desktop/BL1/filtered_stringtie_abundance_BL1.txt",
                 "/home/datavis/Desktop/BL2/filtered_stringtie_abundance_BL2.txt",
                 "/home/datavis/Desktop/BL3/filtered_stringtie_abundance_BL3.txt",
                 "/home/datavis/Desktop/BL4/filtered_stringtie_abundance_BL4.txt")

saliva_files <- c("/home/datavis/Desktop/SL1/filtered_stringtie_abundance_SL1.txt",
                  "/home/datavis/Desktop/SL2/filtered_stringtie_abundance_SL2.txt",
                  "/home/datavis/Desktop/SL3/filtered_stringtie_abundance_SL3.txt")

# Function to load and clean TPM data
load_tpm <- function(file) {
  df <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  df$TPM <- as.numeric(gsub(",", "", df$TPM))
  df %>%
    group_by(Gene.Name) %>%
    summarise(TPM = sum(TPM, na.rm = TRUE)) %>%
    ungroup()
}

# Load and rename blood data
blood_list <- lapply(seq_along(blood_files), function(i) {
  df <- load_tpm(blood_files[i])
  sample_name <- paste0("BL", i)
  df %>% rename(!!sample_name := TPM)
})
blood_merged <- Reduce(function(x, y) full_join(x, y, by = "Gene.Name"), blood_list)

# Load and rename saliva data
saliva_list <- lapply(seq_along(saliva_files), function(i) {
  df <- load_tpm(saliva_files[i])
  sample_name <- paste0("SL", i)
  df %>% rename(!!sample_name := TPM)
})
saliva_merged <- Reduce(function(x, y) full_join(x, y, by = "Gene.Name"), saliva_list)

# Combine all samples
combined_df <- full_join(blood_merged, saliva_merged, by = "Gene.Name")
combined_mat <- as.data.frame(combined_df)
rownames(combined_mat) <- combined_mat$Gene.Name
combined_mat <- combined_mat[, -1]
combined_mat[is.na(combined_mat)] <- 0

# Filter: Keep genes with mean TPM > 2, then top 100 by variance
filtered <- combined_mat[rowMeans(combined_mat) > 2, ]
filtered$variance <- apply(filtered, 1, var, na.rm = TRUE)
top_genes <- filtered %>%
  arrange(desc(variance)) %>%
  slice(1:200) %>%
  select(-variance)
log_mat <- log2(top_genes + 1)

# Create column annotation
sample_labels <- colnames(log_mat)
annotation_col <- data.frame(Group = ifelse(grepl("^BL", sample_labels), "Blood", "Saliva"))
rownames(annotation_col) <- sample_labels
annotation_colors <- list(Group = c(Blood = "#d73027", Saliva = "#4575b4"))

# Reorder columns manually and disable clustering
log_mat <- log_mat[, c("BL1", "BL2", "BL3", "BL4", "SL1", "SL2", "SL3")]
annotation_col <- annotation_col[c("BL1", "BL2", "BL3", "BL4", "SL1", "SL2", "SL3"), , drop = FALSE]

# Define color gradient breaks
breaksList <- seq(0, 15, by = 0.5)

# Save heatmap to PNG
png("heatmap_blood_vs_saliva_top10.png", width = 1200, height = 1000)
pheatmap(log_mat,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)),
         breaks = breaksList,
         main = "Top 200 Most Variable Genes: Blood vs Saliva",
         fontsize_row = 8)
dev.off()
# -- PCA Visualization --

# Perform PCA on the samples (transpose matrix so that samples are rows)
pca_res <- prcomp(t(log_mat), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$Sample <- rownames(pca_df)
# Assign group labels based on sample names: assume BL* are blood and SL* are saliva
pca_df$Group <- ifelse(grepl("^BL", pca_df$Sample), "Blood", "Saliva")

# Plot PCA, coloring samples by Group
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1) +
  labs(title = "PCA of Gene Expression: Blood vs Saliva",
       x = "PC1", y = "PC2") +
  theme_minimal()

