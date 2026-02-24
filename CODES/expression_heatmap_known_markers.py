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
# UPDATED: New filename for the marker-specific heatmap
HEATMAP_FILE = os.path.join(OUTPUT_DIR, "Marker_Genes_Heatmap.png")
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
    Filters a gene count matrix and generates a heatmap for a specific
    list of known marker genes.
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
    final_gene_count = len(filtered_counts_df)
    print(f"Removed {initial_gene_count - final_gene_count} genes.")
    print(f"Remaining genes for analysis: {final_gene_count}")
    # We can skip saving this file again if it's already created
    if not os.path.exists(FILTERED_COUNTS_FILE):
        filtered_counts_df.to_csv(FILTERED_COUNTS_FILE)
        print(f"Saved filtered count data to: {FILTERED_COUNTS_FILE}")

    # --- 3. Group and Reorder Samples ---
    print("\n--- 3. Grouping samples by type (SL, BL, Mixture) ---")
    all_samples = filtered_counts_df.columns.tolist()
    sl_samples = sorted([s for s in all_samples if s.startswith('SL')])
    bl_samples = sorted([s for s in all_samples if s.startswith('BL')])
    mixture_samples = sorted([s for s in all_samples if s not in sl_samples and s not in bl_samples])
    
    grouped_sample_order = sl_samples + bl_samples + mixture_samples
    grouped_counts_df = filtered_counts_df[grouped_sample_order]
    print("Samples have been reordered for the heatmap.")

    # --- 4. Select ONLY Known Marker Genes for the Heatmap ---
    print(f"\n--- 4. Selecting only the known marker genes for the heatmap ---")
    
    # Find which of our defined markers are actually present in the filtered data
    available_markers = [marker for marker in MARKER_GENES_TO_INCLUDE if marker in grouped_counts_df.index]
    
    if not available_markers:
        print("FATAL ERROR: None of the specified marker genes were found in the dataset.")
        return
        
    # Create a dataframe containing only the marker genes
    marker_genes_df = grouped_counts_df.loc[available_markers]
    
    print(f"Found {len(marker_genes_df)} marker genes to plot.")

    # --- 5. Generate and Save Heatmap ---
    print(f"\n--- 5. Generating heatmap for marker genes ---")
    # Z-score normalization for better color scaling (row-wise)
    z_scored_df = marker_genes_df.apply(lambda x: (x - x.mean()) / x.std(), axis=1).fillna(0)

    # Create color mapping for sample types
    sample_types = ['Saliva' if s.startswith('SL') else 'Blood' if s.startswith('BL') else 'Mixture' for s in grouped_sample_order]
    color_map = {'Saliva': 'blue', 'Blood': 'red', 'Mixture': 'grey'}
    col_colors = pd.Series(sample_types, index=grouped_sample_order).map(color_map)

    # Create the clustermap
    g = sns.clustermap(
        z_scored_df,
        cmap="coolwarm",
        yticklabels=True,
        col_cluster=False,
        row_cluster=True,
        col_colors=col_colors,
        figsize=(12, 8), # Adjusted size for fewer genes
        center=0
    )

    # Adjust plot aesthetics
    g.fig.suptitle('Expression of Known Blood and Saliva Markers', fontsize=16)
    g.ax_heatmap.set_xlabel('Samples', fontsize=12)
    g.ax_heatmap.set_ylabel('Marker Genes', fontsize=12)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')
    
    # Add a proper legend for the sample type colors
    handles = [mpatches.Patch(color=color, label=label) for label, color in color_map.items()]
    g.fig.legend(handles=handles, title='Sample Type', bbox_to_anchor=(1.02, 0.9), loc='upper left', borderaxespad=0.)

    # Save the figure with high resolution, adjusting bbox to prevent cropping
    plt.savefig(HEATMAP_FILE, dpi=600, bbox_inches='tight')
    
    print(f"\n--- Analysis Complete ---")
    print(f"Marker gene heatmap saved to: {HEATMAP_FILE}")

if __name__ == "__main__":
    analyze_expression(GENE_COUNT_FILE)
