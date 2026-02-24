# --- User-Defined Variables (EDIT THESE) ---

# Input Files/Directories:
# Assumes you created these files inside /home/datavis/Desktop/lncRNA
SAMPLE_LIST_INITIAL_GTFS="/home/datavis/Desktop/lncRNA/sample_list_stringtie_gtfs.txt"
SAMPLE_LIST_BAMS="/home/datavis/Desktop/lncRNA/sample_list_bams.txt"

# Reference Files (REPLACE <placeholders> with actual paths):
REFERENCE_GTF="/home/datavis/Desktop/gencode.v47.primary_assembly.basic.annotation.gtf"
REFERENCE_FASTA="/home/datavis/Desktop/GRCh38.primary_assembly.genome.fa"
KNOWN_LNCRNA_GTF="/home/datavis/Desktop/gencode.v47.long_noncoding_RNAs.gtf"
KNOWN_CODING_GTF="/home/datavis/Desktop/gencode.v47.primary_assembly.basic.annotation.gtf"
# Output Files/Directories (Script will create these relative to current location, or specify full paths):
OUTPUT_DIR="/home/datavis/Desktop/lncRNA/lncRNA_pipeline_output" # Creates subdirectory in your current location
# Or use a full path like: OUTPUT_DIR="/home/datavis/Desktop/lncRNA/lncRNA_pipeline_output"
MERGED_GTF="${OUTPUT_DIR}/merged_transcripts.gtf"
REQUANT_DIR="${OUTPUT_DIR}/requantified_gtfs"
PREPDE_COUNTS_GENE="${OUTPUT_DIR}/merged_gene_counts.csv"
PREPDE_COUNTS_TX="${OUTPUT_DIR}/merged_transcript_counts.csv" # We'll use this one
LNCRNA_ID_DIR="${OUTPUT_DIR}/lncRNA_identification"
LNCRNA_IDS_FINAL="${LNCRNA_ID_DIR}/final_lncrna_transcript_ids.txt"

# Parameters:
MIN_LNCRNA_LEN=200 # Minimum length for a transcript to be considered lncRNA
CPU_THREADS=11     # Set based on your 'nproc' output

# Path to prepDE.py (or prepDE.py3) - Adjust if not in PATH (REPLACE <placeholder>):
PREPDE_SCRIPT="/home/datavis/miniconda3/pkgs/stringtie-2.1.7-h978d192_0/bin/prepDE.py"

# --- End of User Configuration ---

# --- Create Output Directories ---
mkdir -p ${OUTPUT_DIR}
mkdir -p ${REQUANT_DIR}
mkdir -p ${LNCRNA_ID_DIR}

# --- Set Current Directory (Optional but Recommended) ---
cd /home/datavis/Desktop/lncRNA # Uncomment this line if you want to run everything from this directory

echo "Configuration set. Output will be in: ${OUTPUT_DIR}"
echo "Make sure you have created:"
echo "  ${SAMPLE_LIST_INITIAL_GTFS}"
echo "  ${SAMPLE_LIST_BAMS}"
echo "And filled in all <placeholder> paths above."
