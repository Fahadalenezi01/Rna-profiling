import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np

# --- 1. Configuration ---

# Input file (the gene-level count matrix from your pipeline)
# This script assumes it's in the same directory.
GENE_COUNT_FILE = "./merged_gene_counts.csv"

# Output settings
OUTPUT_DIR = "./expression_analysis_results"
# New filename for the two-gene sex determination plot
SEX_PLOT_FILE = os.path.join(OUTPUT_DIR, "Sex_Determination_Plot_Final.png")

# Genes to check for sex determination
FEMALE_MARKER_GENE = 'XIST'
# --- UPDATED: Using a list of male markers for robustness ---
MALE_MARKER_GENES = ['RPS4Y1', 'DDX3Y', 'EIF1AY', 'KDM5D', 'UTY']


# --- 2. Main Analysis Script ---

def analyze_sex_markers(gene_count_file):
    """
    Analyzes the expression of XIST and a top Y-chromosome gene to
    infer biological sex and generates a grouped bar plot.
    """
    print("--- 1. Loading and Preparing Data ---")

    # Create output directory if it doesn't exist
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print(f"Created output directory: {OUTPUT_DIR}")

    # Load the gene count data
    try:
        counts_df = pd.read_csv(gene_count_file, index_col=0)
        print(f"Successfully loaded {gene_count_file}. Shape: {counts_df.shape}")
    except FileNotFoundError:
        print(f"FATAL ERROR: The file '{gene_count_file}' was not found.")
        return

    # --- 2. Clean Gene Names and Find Marker Genes ---
    print(f"\n--- 2. Searching for {FEMALE_MARKER_GENE} and male marker genes ---")
    
    def get_clean_gene_name(raw_name):
        if '|' in raw_name:
            return raw_name.split('|')[1]
        return raw_name
    
    counts_df.index = [get_clean_gene_name(name) for name in counts_df.index]
    counts_df = counts_df[~counts_df.index.duplicated(keep='first')]

    # --- UPDATED: Find the best available male marker ---
    available_male_markers = [gene for gene in MALE_MARKER_GENES if gene in counts_df.index]
    
    if not available_male_markers:
        print("WARNING: No male marker genes found. Plot will only show XIST.")
        selected_male_marker = None
    else:
        # Select the male marker with the highest average expression
        mean_male_expr = counts_df.loc[available_male_markers].mean(axis=1)
        selected_male_marker = mean_male_expr.idxmax()
        print(f"Found available male markers: {available_male_markers}")
        print(f"Selected '{selected_male_marker}' as the representative male marker.")

    # Check if XIST is present
    if FEMALE_MARKER_GENE not in counts_df.index:
        print(f"WARNING: Female marker gene '{FEMALE_MARKER_GENE}' not found.")
        genes_to_plot = [selected_male_marker] if selected_male_marker else []
    else:
        genes_to_plot = [FEMALE_MARKER_GENE] + ([selected_male_marker] if selected_male_marker else [])

    if not genes_to_plot:
        print("FATAL ERROR: No sex marker genes could be found in the data.")
        return

    # Extract the expression data for the final marker genes
    marker_data = counts_df.loc[genes_to_plot].T
    
    # Use log2(counts + 1) for better visualization
    log_marker_data = np.log2(marker_data + 1)
    
    # Prepare data for grouped bar plot
    plot_data = log_marker_data.reset_index().melt(id_vars='index', var_name='Gene', value_name='Expression')
    plot_data.rename(columns={'index': 'Sample'}, inplace=True)

    print("Found marker gene expression data.")

    # --- 3. Generate and Save Grouped Bar Plot ---
    print(f"\n--- 3. Generating grouped bar plot for sex marker expression ---")

    plt.figure(figsize=(14, 8))
    
    # Define a color palette
    palette = {FEMALE_MARKER_GENE: 'hotpink'}
    if selected_male_marker:
        palette[selected_male_marker] = 'deepskyblue'

    # Create the grouped bar plot
    sns.barplot(data=plot_data, x='Sample', y='Expression', hue='Gene', palette=palette)
    
    plt.title(f'Sex Marker Expression ({", ".join(genes_to_plot)})', fontsize=16)
    plt.xlabel('Sample Name', fontsize=12)
    plt.ylabel('log2(Counts + 1)', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='Gene')
    plt.tight_layout()
    
    # Save the figure with high resolution
    plt.savefig(SEX_PLOT_FILE, dpi=600)
    
    print(f"\n--- Analysis Complete ---")
    print(f"Sex determination plot saved to: {SEX_PLOT_FILE}")


if __name__ == "__main__":
    analyze_sex_markers(GENE_COUNT_FILE)
