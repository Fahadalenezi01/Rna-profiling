#!/bin/bash
set -e # Exit immediately if any command fails

# ==============================================================================
#                           CONFIGURATION
# This script is now configured with the final, correct paths for all samples.
# ==============================================================================

# 1. REFERENCE GENOME
#    The path to your reference genome on the D: drive.
REFERENCE_GENOME="/mnt/d/GRCh38.primary_assembly.genome.fa"

# 2. MAIN OUTPUT DIRECTORY
#    A main folder where all sample-specific result folders will be created.
OUTPUT_DIR="~/analysis_pipeline_results"

# 3. SAMPLE CONFIGURATION
#    The list of your five samples and the full WSL paths to their input BAM files.
SAMPLE_NAMES=(
    "SL3_rep"
    "SL4"
    "50SL50BL"
    "70SL30BL"
    "30SL70BL"
)

# UPDATED: The paths for SL3_rep and SL4 now point to the new, clean files.
INPUT_BAM_PATHS=(
    "~/re_basecalled_results/SL3_rep_calls.bam"
    "~/re_basecalled_results/SL4_calls.bam"
    "~/dorado_analysis/50SL50BL.bam"
    "~/dorado_analysis/70SL30BL.bam"
    "~/dorado_analysis/30SL70BL.bam"
)

# 4. THREADS
#    Number of CPU threads to use for alignment and assembly.
THREADS=20


# ==============================================================================
#                        TOOL INSTALLATION
# ==============================================================================
# This script handles its own Python environment to avoid system errors.
# You just need to ensure the base system tools are installed once.
#
# sudo apt-get update
# sudo apt-get install -y samtools minimap2 python3-pip python3.12-venv
#
# ==============================================================================
#                           ANALYSIS PIPELINE
# ==============================================================================

echo "--- Starting Post-Basecalling Analysis Pipeline ---"

# --- Step 1: Force a clean setup of the Python Virtual Environment ---
VENV_DIR="analysis_venv"
echo "[INFO] Forcing a clean virtual environment setup..."
rm -rf "$VENV_DIR"
echo "[INFO] Creating new Python virtual environment in './${VENV_DIR}'..."
sudo apt-get install -y python3.12-venv
python3 -m venv "$VENV_DIR"
echo "[INFO] Activating environment and installing specific package versions..."
source "${VENV_DIR}/bin/activate"
pip3 install "NanoPlot==1.44.1" "kaleido==0.2.1"
deactivate
echo "[INFO] Environment setup complete."

NANOPLOT_CMD=("${VENV_DIR}/bin/python3" "-m" "nanoplot.NanoPlot")

# --- Step 2: Create Output Directory and Verify Files ---
eval mkdir -p "$OUTPUT_DIR"
echo "[INFO] Main output directory is at: $OUTPUT_DIR"
if [ ! -f "$REFERENCE_GENOME" ]; then
    echo "!!! CRITICAL ERROR: Reference Genome file not found at $REFERENCE_GENOME" && exit 1
fi
echo "[INFO] Successfully found reference genome."

# --- Step 3: Loop through each sample and run the pipeline ---
echo ""
echo "========================================="
echo "      Starting Analysis Loop"
echo "========================================="

for i in "${!SAMPLE_NAMES[@]}"; do
    sample_name="${SAMPLE_NAMES[i]}"
    eval input_bam_path="${INPUT_BAM_PATHS[i]}"

    echo ""
    echo "--- Processing Sample: ${sample_name} ---"

    eval sample_output_dir="${OUTPUT_DIR}/${sample_name}"
    mkdir -p "$sample_output_dir"
    echo "[INFO] Results for this sample will be saved in: ${sample_output_dir}"

    QC_DIR="${sample_output_dir}/QC_Report"
    ALIGNED_BAM="${sample_output_dir}/${sample_name}.aligned.bam"
    SORTED_BAM="${sample_output_dir}/${sample_name}.aligned.sorted.bam"
    ALIGNMENT_STATS="${sample_output_dir}/${sample_name}.alignment_stats.txt"

    if [ ! -f "$input_bam_path" ]; then
        echo "!!! WARNING: Input BAM for sample ${sample_name} not found at ${input_bam_path}. Skipping."
        continue
    fi

    # --- Step 1 - Check BAM file integrity ---
    echo "[1/5] Verifying integrity of input BAM file..."
    set +o pipefail
    validation_output=$(samtools view -h "$input_bam_path" > /dev/null 2>&1; echo $?)
    set -o pipefail
    if [ "$validation_output" -ne 0 ]; then
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "!!! WARNING: Input BAM for sample '${sample_name}' appears to be corrupted or truncated."
        echo "!!! Skipping this sample. You may need to re-run basecalling for it."
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        continue
    fi
    echo "[INFO] BAM file integrity check passed."

    # --- Step 2: Alignment with minimap2 ---
    echo "[2/5] Aligning reads to reference genome..."
    # CORRECTED: Changed '-ay splice' to the correct '-ax splice' option.
    samtools fastq "$input_bam_path" | minimap2 -ax splice -N 20 -t "$THREADS" "$REFERENCE_GENOME" - | samtools view -bS - > "$ALIGNED_BAM"

    # --- Step 3: Sort and Index the Aligned Reads ---
    echo "[3/5] Sorting and indexing aligned reads..."
    samtools sort -@ "$THREADS" -o "$SORTED_BAM" "$ALIGNED_BAM"
    samtools index "$SORTED_BAM"
    
    # --- Step 4: Quality Control on ALIGNED reads ---
    echo "[4/5] Running Quality Control (NanoPlot) on aligned reads..."
    "${NANOPLOT_CMD[@]}" --bam "$SORTED_BAM" -o "$QC_DIR" --title "${sample_name}_aligned"

    # --- Step 5: Generate Alignment Statistics ---
    echo "[5/5] Generating final alignment statistics..."
    samtools flagstat "$SORTED_BAM" > "$ALIGNMENT_STATS"
    
    rm "$ALIGNED_BAM"

    echo "--- Finished processing sample ${sample_name}. ---"
done

echo ""
echo "========================================="
echo "      All Samples Processed!"
echo "========================================="
echo "Pipeline finished successfully."
echo "You are the smartest person i have ever worked with I am so glad you are my owner"
