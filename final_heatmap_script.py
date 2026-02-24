import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.patches as mpatches # Import for creating custom legends

# --- 1. Configuration ---

# Input file (the gene-level count matrix from your pipeline)
# This script assumes it's in the same directory.
GENE_COUNT_FILE = "./merged_gene_counts.csv"

# Output settings
OUTPUT_DIR = "./expression_analysis_results"
# UPDATED: New filenames for the separate heatmaps
HEATMAP_TOP50_PURE_FILE = os.path.join(OUTPUT_DIR, "Top_50_Genes_Heatmap_Pure_Samples.png")
HEATMAP_MARKERS_MIXTURE_FILE = os.path.join(OUTPUT_DIR, "Marker_Genes_Heatmap_Mixture_Samples.png")
FILTERED_COUNTS_FILE = os.path.join(OUTPUT_DIR, "filtered_gene_counts.csv")

# Number of top expressed genes to plot for the pure samples heatmap
NUM_TOP_GENES = 50

# --- 2. Define Gene Lists for Filtering and Analysis ---

# These lists contain common prefixes or symbols for genes to be removed.
RRNA_PREFIXES = ['RNA5S', 'RNA5-8S', 'RNA18S', 'RNA28S']
TRNA_SUBSTRINGS = ['TRN']
MITOCHONDRIAL_PREFIXES = ['MT-']
HOUSEKEEPING_GENES = [
    'ACTB', 'GAPDH', 'B2M', 'UBC', 'HPRT1', 'TBP', 'RPL13A', 'RPS18',
    'PPIA', 'GUSB', 'YWHAZ', 'PGK1', 'HMBS', 'SDHA'
]

# Refined marker gene lists based on data observations
BLOOD_PBMC_MARKERS = ['CD3D', 'CD8A', 'MS4A1', 'CD19', 'NCAM1', 'CD14']
SALIVA_MARKERS = ['AMY1A', 'AMY1B', 'MUC7', 'BPIFA1', 'STATH']
MARKER_GENES_TO_INCLUDE = BLOOD_PBMC_MARKERS + SALIVA_MARKERS


# --- 3. Main Analysis Script ---

def analyze_expression(gene_count_file):
    """
    Generates two separate heatmaps:
    1. Top 50 expressed genes in pure samples.
    2. Marker gene expression in mixture samples.
    """
    print("--- 1. Loading and Preparing Data ---")

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print(f"Created output directory: {OUTPUT_DIR}")

    try:
        counts_df = pd.read_csv(gene_count_file, index_col=0)
        print(f"Successfully loaded {gene_count_file}. Shape: {counts_df.shape}")
    except FileNotFoundError:
        print(f"FATAL ERROR: The file '{gene_count_file}' was not found.")
        return

    # --- 2. Filtering and Cleaning Gene Names ---
    print("\n--- 2. Filtering unwanted genes and cleaning names ---")
    
    def get_clean_gene_name(raw_name):
        if '|' in raw_name:
            return raw_name.split('|')[1]
        return raw_name
    
    counts_df.index = [get_clean_gene_name(name) for name in counts_df.index]
    counts_df = counts_df[~counts_df.index.duplicated(keep='first')]

    genes_to_remove = counts_df.index.isin(HOUSEKEEPING_GENES) | counts_df.index.str.startswith(tuple(RRNA_PREFIXES))
    filtered_counts_df = counts_df[~genes_to_remove]
    
    # --- 3. Define Final Sample Groups (Excluding BL1 Outlier) ---
    print("\n--- 3. Defining final sample groups and excluding BL1 outlier ---")
    
    users_saliva_samples = sorted(['SL3', 'SL3_rep', 'SL4'])
    commercial_saliva_samples = sorted(['SL1', 'SL2'])
    blood_pbmc_samples = sorted(['BL2', 'BL3', 'BL4'])
    mixture_samples = sorted([s for s in filtered_counts_df.columns if s not in users_saliva_samples + commercial_saliva_samples + blood_pbmc_samples + ['BL1']])
    
    pure_samples = users_saliva_samples + commercial_saliva_samples + blood_pbmc_samples
    
    pure_samples_df = filtered_counts_df[pure_samples]
    mixture_samples_df = filtered_counts_df[mixture_samples]

    # --- 4. Generate Heatmap 1: Top 50 Genes in PURE Samples ---
    print(f"\n--- 4. Generating heatmap for Top 50 genes in PURE samples ---")
    
    mean_expression_pure = pure_samples_df.mean(axis=1)
    top_50_genes_pure = mean_expression_pure.sort_values(ascending=False).head(NUM_TOP_GENES)
    top_genes_df = pure_samples_df.loc[top_50_genes_pure.index]

    z_scored_df = top_genes_df.apply(lambda x: (x - x.mean()) / x.std(), axis=1).fillna(0)

    pure_sample_types = []
    for s in pure_samples:
        if s in users_saliva_samples: pure_sample_types.append("User's Saliva")
        elif s in commercial_saliva_samples: pure_sample_types.append("Commercial Saliva")
        else: pure_sample_types.append("Blood (PBMC)")
        
    pure_color_map = {"User's Saliva": 'cyan', "Commercial Saliva": 'blue', "Blood (PBMC)": 'red'}
    pure_col_colors = pd.Series(pure_sample_types, index=pure_samples).map(pure_color_map)

    g_pure = sns.clustermap(
        z_scored_df, cmap="coolwarm", yticklabels=True, col_cluster=False, row_cluster=True,
        col_colors=pure_col_colors, figsize=(12, 14), center=0
    )
    g_pure.fig.suptitle('Top 50 Expressed Genes in Pure Samples', fontsize=16)
    plt.savefig(HEATMAP_TOP50_PURE_FILE, dpi=600, bbox_inches='tight')
    print(f"Pure samples heatmap saved to: {HEATMAP_TOP50_PURE_FILE}")
    plt.close()

    # --- 5. Generate Heatmap 2: Marker Genes in MIXTURE Samples ---
    print(f"\n--- 5. Generating heatmap for Marker genes in MIXTURE samples ---")
    
    available_markers = [marker for marker in MARKER_GENES_TO_INCLUDE if marker in mixture_samples_df.index]
    marker_genes_mixture_df = mixture_samples_df.loc[available_markers]
    z_scored_mixture_df = marker_genes_mixture_df.apply(lambda x: (x - x.mean()) / x.std(), axis=1).fillna(0)

    g_mixture = sns.clustermap(
        z_scored_mixture_df, cmap="coolwarm", yticklabels=True, col_cluster=False, row_cluster=True,
        figsize=(8, 8), center=0
    )
    g_mixture.fig.suptitle('Expression of Markers in Mixture Samples', fontsize=16)
    plt.savefig(HEATMAP_MARKERS_MIXTURE_FILE, dpi=600, bbox_inches='tight')
    print(f"Mixture samples heatmap saved to: {HEATMAP_MARKERS_MIXTURE_FILE}")
    plt.close()

    print(f"\n--- Analysis Complete ---")

if __name__ == "__main__":
    analyze_expression(GENE_COUNT_FILE)
