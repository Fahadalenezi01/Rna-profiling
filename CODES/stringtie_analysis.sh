#!/bin/bash
set -e # Exit immediately if any command fails

# ==============================================================================
#                           CONFIGURATION
# This script is designed to run StringTie on existing aligned BAM files.
# ==============================================================================

# 1. REFERENCE ANNOTATION (Required for StringTie)
#    The path to your reference annotation file in GTF or GFF3 format.
REFERENCE_ANNOTATION_GTF="/mnt/d/gencode.v47.primary_assembly.basic.annotation.gtf"

# 2. ALIGNED DATA DIRECTORY
#    The main folder where your sample-specific result folders are located.
ALIGNED_DATA_DIR="~/analysis_pipeline_results"

# 3. STRINGTIE OUTPUT DIRECTORY
#    A new folder where the StringTie results will be saved.
STRINGTIE_OUTPUT_DIR="~/stringtie_quantification_results"

# 4. SAMPLE NAMES
#    The list of your five samples.
SAMPLE_NAMES=(
    "SL3_rep"
    "SL4"
    "50SL50BL"
    "70SL30BL"
    "30SL70BL"
)

# 5. THREADS
#    Number of CPU threads to use for assembly.
THREADS=20


# ==============================================================================
#                        TOOL INSTALLATION (Run this section once)
# ==============================================================================
#
# If you haven't already, you only need to install stringtie.
#
# sudo apt-get update
# sudo apt-get install -y stringtie
#
# ==============================================================================
#                           QUANTIFICATION PIPELINE
# ==============================================================================

echo "--- Starting StringTie Quantification Pipeline ---"

# --- Create the main output directory ---
eval mkdir -p "$STRINGTIE_OUTPUT_DIR"
echo "[INFO] StringTie output directory is at: $STRINGTIE_OUTPUT_DIR"

# --- Verify Reference Annotation File Exists ---
if [ ! -f "$REFERENCE_ANNOTATION_GTF" ]; then
    echo "!!! CRITICAL ERROR: Reference Annotation GTF file not found at $REFERENCE_ANNOTATION_GTF" && exit 1
fi
echo "[INFO] Successfully found reference annotation file."


# --- Loop through each sample and run StringTie ---
echo ""
echo "========================================="
echo "      Starting Quantification Loop"
echo "========================================="

for sample_name in "${SAMPLE_NAMES[@]}"; do
    echo ""
    echo "--- Processing Sample: ${sample_name} ---"

    # --- Define input and output paths for this sample ---
    eval input_sorted_bam="${ALIGNED_DATA_DIR}/${sample_name}/${sample_name}.aligned.sorted.bam"
    output_gtf="${STRINGTIE_OUTPUT_DIR}/${sample_name}.transcripts.gtf"

    # --- Verify the input sorted BAM file exists ---
    if [ ! -f "$input_sorted_bam" ]; then
        echo "!!! WARNING: Input sorted BAM for sample ${sample_name} not found at ${input_sorted_bam}. Skipping."
        continue
    fi
    echo "[INFO] Found input file: ${input_sorted_bam}"

    # --- Run StringTie ---
    echo "[INFO] Assembling transcripts with StringTie..."
    # The -L flag is important for long reads like Nanopore.
    # The -G flag provides the reference annotation to guide assembly.
    stringtie "$input_sorted_bam" -L -p "$THREADS" -G "$REFERENCE_ANNOTATION_GTF" -o "$output_gtf"

    echo "--- Finished processing sample ${sample_name}. ---"
    echo "Final assembled transcripts saved to: ${output_gtf}"
done

echo ""
echo "========================================="
echo "      All Samples Processed!"
echo "========================================="
echo "Pipeline finished successfully."
