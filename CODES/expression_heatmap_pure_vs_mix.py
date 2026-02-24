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
# Define output filenames for the two separate heatmaps
HEATMAP_PURE_SAMPLES_FILE = os.path.join(OUTPUT_DIR, "Marker_Genes_Heatmap_Pure_Samples.png")
HEATMAP_MIXTURE_SAMPLES_FILE = os.path.join(OUTPUT_DIR, "Marker_Genes_Heatmap_Mixture_Samples.png")
FILTERED_COUNTS_FILE = os.path.join(OUTPUT_DIR, "filtered_gene_counts.csv")

# --- 2. Define Gene Lists for Filtering and Analysis ---

# These lists contain common prefixes or symbols for genes to be removed.
RRNA_PREFIXES = ['RNA5S', 'RNA5-8S', 'RNA18S', 'RNA28S']
TRNA_SUBSTRINGS = ['TRN']
MITOCHONDRIAL_PREFIXES = ['MT-']
HOUSEKEEPING_GENES = [
    'ACTB', 'GAPDH', 'B2M', 'UBC', 'HPRT1', 'TBP', 'RPL13A', 'RPS18',
    'PPIA', 'GUSB', 'YWHAZ', 'PGK1', 'HMBS', 'SDHA'
]

# Define known marker genes to ensure they are included in the heatmap
BLOOD_MARKERS = ['HBB', 'HBA1', 'HBA2', 'ALAS2']
SALIVA_MARKERS = ['AMY1A', 'AMY1B', 'MUC7', 'BPIFA1']
MARKER_GENES_TO_INCLUDE = BLOOD_MARKERS + SALIVA_MARKERS


# --- 3. Main Analysis Script ---

def analyze_expression(gene_count_file):
    """
    Filters a gene count matrix and generates two separate heatmaps:
    1. For pure Saliva vs. Blood samples.
    2. For Mixture samples only.
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

    # --- 2. Filtering Unwanted Genes ---
    print("\n--- 2. Filtering out rRNA, tRNA, MT-RNA, and housekeeping genes ---")
    initial_gene_count = len(counts_df)
    
    # Clean gene names before filtering
    def get_clean_gene_name(raw_name):
        if '|' in raw_name:
            return raw_name.split('|')[1]
        return raw_name
    
    counts_df.index = [get_clean_gene_name(name) for name in counts_df.index]
    # Remove duplicate gene names that might arise after cleaning, keeping the first one
    counts_df = counts_df[~counts_df.index.duplicated(keep='first')]

    is_rrna = counts_df.index.str.startswith(tuple(RRNA_PREFIXES))
    is_trna = pd.Series([any(sub in idx for sub in TRNA_SUBSTRINGS) for idx in counts_df.index], index=counts_df.index)
    is_mito = counts_df.index.str.startswith(tuple(MITOCHONDRIAL_PREFIXES))
    is_housekeeping = counts_df.index.isin(HOUSEKEEPING_GENES)
    genes_to_remove = is_rrna | is_trna | is_mito | is_housekeeping
    filtered_counts_df = counts_df[~genes_to_remove]
    
    # --- 3. Separate Pure and Mixture Samples ---
    print("\n--- 3. Separating Pure and Mixture Samples ---")
    all_samples = filtered_counts_df.columns.tolist()
    sl_samples = sorted([s for s in all_samples if s.startswith('SL')])
    bl_samples = sorted([s for s in all_samples if s.startswith('BL')])
    mixture_samples = sorted([s for s in all_samples if s not in sl_samples and s not in bl_samples])
    
    pure_samples_df = filtered_counts_df[sl_samples + bl_samples]
    mixture_samples_df = filtered_counts_df[mixture_samples]

    # --- 4. Select Marker Genes ---
    available_markers = [marker for marker in MARKER_GENES_TO_INCLUDE if marker in filtered_counts_df.index]
    if not available_markers:
        print("FATAL ERROR: None of the specified marker genes were found in the dataset.")
        return
    
    # --- 5. Generate Heatmap for PURE Samples ---
    print(f"\n--- 5. Generating heatmap for PURE samples (Saliva vs. Blood) ---")
    
    marker_genes_pure_df = pure_samples_df.loc[available_markers]
    z_scored_pure_df = marker_genes_pure_df.apply(lambda x: (x - x.mean()) / x.std(), axis=1).fillna(0)
    
    pure_sample_types = ['Saliva' if s.startswith('SL') else 'Blood' for s in z_scored_pure_df.columns]
    pure_color_map = {'Saliva': 'blue', 'Blood': 'red'}
    pure_col_colors = pd.Series(pure_sample_types, index=z_scored_pure_df.columns).map(pure_color_map)

    g_pure = sns.clustermap(
        z_scored_pure_df,
        cmap="coolwarm", yticklabels=True, col_cluster=False, row_cluster=True,
        col_colors=pure_col_colors, figsize=(10, 8), center=0
    )
    g_pure.fig.suptitle('Expression of Markers in Pure Samples', fontsize=16)
    plt.savefig(HEATMAP_PURE_SAMPLES_FILE, dpi=600, bbox_inches='tight')
    print(f"Pure samples heatmap saved to: {HEATMAP_PURE_SAMPLES_FILE}")
    plt.close() # Close the figure to prepare for the next one

    # --- 6. Generate Heatmap for MIXTURE Samples ---
    print(f"\n--- 6. Generating heatmap for MIXTURE samples ---")
    
    marker_genes_mixture_df = mixture_samples_df.loc[available_markers]
    z_scored_mixture_df = marker_genes_mixture_df.apply(lambda x: (x - x.mean()) / x.std(), axis=1).fillna(0)
    
    g_mixture = sns.clustermap(
        z_scored_mixture_df,
        cmap="coolwarm", yticklabels=True, col_cluster=False, row_cluster=True,
        figsize=(8, 8), center=0
    )
    g_mixture.fig.suptitle('Expression of Markers in Mixture Samples', fontsize=16)
    plt.savefig(HEATMAP_MIXTURE_SAMPLES_FILE, dpi=600, bbox_inches='tight')
    print(f"Mixture samples heatmap saved to: {HEATMAP_MIXTURE_SAMPLES_FILE}")
    plt.close()

    print(f"\n--- Analysis Complete ---")

if __name__ == "__main__":
    analyze_expression(GENE_COUNT_FILE)
