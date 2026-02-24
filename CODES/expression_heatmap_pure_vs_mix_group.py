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
# UPDATED: New filename for the final data-driven heatmap
HEATMAP_FILE = os.path.join(OUTPUT_DIR, "Data_Driven_Groups_Heatmap.png")
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

# --- UPDATED: Refined marker gene lists based on data observations ---
# Using immune cell markers for blood (PBMCs) instead of hemoglobin
BLOOD_PBMC_MARKERS = ['CD3D', 'CD8A', 'MS4A1', 'CD19', 'NCAM1', 'CD14']
SALIVA_MARKERS = ['AMY1A', 'AMY1B', 'MUC7', 'BPIFA1', 'STATH'] # Added Statherin
MARKER_GENES_TO_INCLUDE = BLOOD_PBMC_MARKERS + SALIVA_MARKERS


# --- 3. Main Analysis Script ---

def analyze_expression(gene_count_file):
    """
    Generates a heatmap based on data-driven sample groupings to clarify
    expression patterns of known marker genes.
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
    
    # --- 3. Define Data-Driven Sample Groups ---
    print("\n--- 3. Defining sample groups based on observed data patterns ---")
    
    # Re-assigning samples to groups based on our analysis
    users_saliva_samples = sorted(['SL3', 'SL3_rep', 'SL4', 'BL1'])
    commercial_saliva_samples = sorted(['SL1', 'SL2'])
    blood_pbmc_samples = sorted(['BL2', 'BL3', 'BL4'])
    mixture_samples = sorted([s for s in filtered_counts_df.columns if s not in users_saliva_samples + commercial_saliva_samples + blood_pbmc_samples])
    
    # Define the final order for the heatmap
    final_sample_order = users_saliva_samples + commercial_saliva_samples + blood_pbmc_samples + mixture_samples
    reordered_df = filtered_counts_df[final_sample_order]

    # --- 4. Select Marker Genes for the Heatmap ---
    print(f"\n--- 4. Selecting marker genes for the heatmap ---")
    available_markers = [marker for marker in MARKER_GENES_TO_INCLUDE if marker in reordered_df.index]
    marker_genes_df = reordered_df.loc[available_markers]
    
    print(f"Found {len(marker_genes_df)} marker genes to plot.")

    # --- 5. Generate Final Heatmap ---
    print(f"\n--- 5. Generating final heatmap with data-driven groups ---")
    
    z_scored_df = marker_genes_df.apply(lambda x: (x - x.mean()) / x.std(), axis=1).fillna(0)
    
    # Create color mapping for the new, data-driven groups
    sample_types = []
    for s in final_sample_order:
        if s in users_saliva_samples: sample_types.append("User's Saliva")
        elif s in commercial_saliva_samples: sample_types.append("Commercial Saliva")
        elif s in blood_pbmc_samples: sample_types.append("Blood (PBMC)")
        else: sample_types.append("Mixture")
        
    color_map = {"User's Saliva": 'cyan', "Commercial Saliva": 'blue', "Blood (PBMC)": 'red', "Mixture": 'grey'}
    col_colors = pd.Series(sample_types, index=final_sample_order).map(color_map)

    g = sns.clustermap(
        z_scored_df,
        cmap="coolwarm", yticklabels=True, col_cluster=False, row_cluster=True,
        col_colors=col_colors, figsize=(12, 10), center=0
    )
    
    g.fig.suptitle('Expression of Markers in Data-Driven Sample Groups', fontsize=16)
    plt.savefig(HEATMAP_FILE, dpi=600, bbox_inches='tight')
    print(f"Final heatmap saved to: {HEATMAP_FILE}")

    print(f"\n--- Analysis Complete ---")

if __name__ == "__main__":
    analyze_expression(GENE_COUNT_FILE)
