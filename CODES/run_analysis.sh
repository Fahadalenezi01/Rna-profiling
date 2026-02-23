#!/bin/bash
set -e # Exit immediately if any command fails

# ==============================================================================
#                           CONFIGURATION
# ==============================================================================
REFERENCE_GENOME="/mnt/d/GRCh38.primary_assembly.genome.fa"
OUTPUT_DIR="~/analysis_pipeline_results"
SAMPLE_NAMES=(
    "SL3_rep"
    "SL4"
    "50SL50BL"
    "70SL30BL"
    "30SL70BL"
)
INPUT_BAM_PATHS=(
    "/mnt/d/SL3_rep/SL3_rep_calls.bam"
    "/mnt/d/SL4/SL4_calls.bam"
    "~/dorado_analysis/50SL50BL.bam"
    "~/dorado_analysis/70SL30BL.bam"
    "~/dorado_analysis/30SL70BL.bam"
)
THREADS=20

# ==============================================================================
#                        TOOL INSTALLATION
# ==============================================================================
# This script now handles its own Python environment to avoid system errors.
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

# Forcefully remove any old, broken environment to ensure a clean start
rm -rf "$VENV_DIR"

echo "[INFO] Creating new Python virtual environment in './${VENV_DIR}'..."
# Ensure the venv package is installed
sudo apt-get install -y python3.12-venv
python3 -m venv "$VENV_DIR"

echo "[INFO] Activating environment and installing specific package versions..."
# Activate the environment and install packages inside it
source "${VENV_DIR}/bin/activate"

# Install specific, known-stable versions to avoid compatibility issues
# CORRECTED: Using an available version of kaleido
pip3 install "NanoPlot==1.44.1" "kaleido==0.2.1"

# Deactivate the environment now that setup is done
deactivate
echo "[INFO] Environment setup complete."


# Define the robust command to run NanoPlot using the venv's python interpreter
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

    # --- Run QC using the NanoPlot from our virtual environment ---
    echo "[1/4] Running Quality Control (NanoPlot)..."
    "${NANOPLOT_CMD[@]}" --bam "$input_bam_path" -o "$QC_DIR" --title "$sample_name"

    # --- Run Alignment ---
    echo "[2/4] Aligning reads to reference genome..."
    minimap2 -ax map-ont -t "$THREADS" "$REFERENCE_GENOME" "$input_bam_path" | samtools view -bS - > "$ALIGNED_BAM"

    # --- Run Sorting ---
    echo "[3/4] Sorting and indexing aligned reads..."
    samtools sort -@ "$THREADS" -o "$SORTED_BAM" "$ALIGNED_BAM"
    samtools index "$SORTED_BAM"

    # --- Run Statistics ---
    echo "[4/4] Generating alignment statistics..."
    samtools flagstat "$SORTED_BAM" > "$ALIGNMENT_STATS"
    
    rm "$ALIGNED_BAM"

    echo "--- Finished processing sample ${sample_name}. ---"
done

echo ""
echo "========================================="
echo "      All Samples Processed!"
echo "========================================="
echo "Pipeline finished successfully."
