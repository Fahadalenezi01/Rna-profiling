import pandas as pd
import numpy as np
import os
import umap
import matplotlib.pyplot as plt
import seaborn as sns

# --- 1. Configuration ---

# Input file (the gene-level count matrix from your pipeline)
# This script assumes it's in the same directory.
GENE_COUNT_FILE = "./merged_gene_counts.csv"

# Output settings
OUTPUT_DIR = "./expression_analysis_results"
# UPDATED: New filename for the static, labeled 2D UMAP plot
UMAP_PLOT_FILE = os.path.join(OUTPUT_DIR, "Static_2D_UMAP_Plot.png")

# UMAP parameters (can be tuned for different visualizations)
N_NEIGHBORS = 5  # Controls how UMAP balances local vs. global structure
MIN_DIST = 0.3   # Controls how tightly UMAP is allowed to pack points together

# --- 2. Main Analysis Script ---

def create_umap_plot(gene_count_file):
    """
    Performs UMAP dimensionality reduction and generates a static,
    labeled 2D scatter plot suitable for publication.
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

    # --- 2. Define Final Data-Driven Sample Groups ---
    print("\n--- 2. Defining sample groups based on previous analysis ---")
    
    # Transpose the dataframe so samples are rows and genes are columns
    data_for_umap = counts_df.T
    
    # Define sample groups
    users_saliva_samples = [s for s in data_for_umap.index if s in ['SL3', 'SL3_rep', 'SL4']]
    commercial_saliva_samples = [s for s in data_for_umap.index if s in ['SL1', 'SL2']]
    blood_pbmc_samples = [s for s in data_for_umap.index if s in ['BL2', 'BL3', 'BL4']]
    technical_outlier = ['BL1']
    mixture_samples = [s for s in data_for_umap.index if s.startswith(('30SL', '50SL', '70SL'))]

    # Assign a group to each sample
    group_map = {}
    for s in data_for_umap.index:
        if s in users_saliva_samples: group_map[s] = "User's Saliva"
        elif s in commercial_saliva_samples: group_map[s] = "Commercial Saliva"
        elif s in blood_pbmc_samples: group_map[s] = "Blood (PBMC)"
        elif s in technical_outlier: group_map[s] = "Technical Outlier (BL1)"
        elif s in mixture_samples: group_map[s] = "Mixture"
        else: group_map[s] = "Unknown"
        
    data_for_umap['SampleGroup'] = data_for_umap.index.map(group_map)

    # --- 3. Normalize Data and Run UMAP ---
    print("\n--- 3. Normalizing data and running UMAP ---")
    
    # Use a log2(x+1) transformation on the count data (features)
    numeric_data = data_for_umap.drop(columns=['SampleGroup'])
    log_transformed_data = np.log2(numeric_data + 1)

    # Initialize and run the UMAP algorithm for 2 dimensions
    reducer = umap.UMAP(n_neighbors=N_NEIGHBORS, min_dist=MIN_DIST, n_components=2, random_state=42)
    embedding = reducer.fit_transform(log_transformed_data)

    # Create a new dataframe with the UMAP results
    umap_df = pd.DataFrame(embedding, columns=['UMAP1', 'UMAP2'])
    umap_df['Sample'] = data_for_umap.index
    umap_df['Group'] = data_for_umap['SampleGroup']
    
    print("UMAP calculation complete.")

    # --- 4. Generate Static 2D Plot with Seaborn ---
    print("\n--- 4. Generating static 2D plot ---")

    plt.figure(figsize=(12, 10))
    
    # Create the scatter plot
    ax = sns.scatterplot(
        data=umap_df,
        x='UMAP1',
        y='UMAP2',
        hue='Group', # Color points by group
        s=100, # Marker size
        alpha=0.8
    )

    # Add labels to each point
    for i, row in umap_df.iterrows():
        plt.text(row['UMAP1'] + 0.1, row['UMAP2'], row['Sample'], fontsize=9)

    # Set labels, title, and legend
    plt.title('2D UMAP Projection of Samples', fontsize=16)
    plt.xlabel('UMAP Dimension 1', fontsize=12)
    plt.ylabel('UMAP Dimension 2', fontsize=12)
    plt.legend(title='Sample Groups', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True)
    plt.tight_layout()

    # Save the plot as a high-resolution PNG file
    plt.savefig(UMAP_PLOT_FILE, dpi=300, bbox_inches='tight')

    print(f"\n--- Analysis Complete ---")
    print(f"Static 2D UMAP plot saved to: {UMAP_PLOT_FILE}")


# --- Corrected main execution block ---
if __name__ == "__main__":
    create_umap_plot(GENE_COUNT_FILE)
