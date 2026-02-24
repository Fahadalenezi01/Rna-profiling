#!/bin/bash
set -e # Exit immediately if any command fails

# ==============================================================================
#                           CONFIGURATION
# This script is configured to re-basecall only the specified samples.
# ==============================================================================

# 1. Base directory on your D: drive where the sample folders are located.
#    UPDATED: The path now correctly points to your "Fahad" folder.
BASE_INPUT_DIR="/mnt/d/Fahad"

# 2. List of the specific samples you want to re-basecall.
SAMPLES=(
    "SL3_rep"
    "SL4"
)

# 3. A new output directory to store the results of this specific run.
OUTPUT_DIR="~/re_basecalled_results"

# ==============================================================================
#                           SCRIPT
# You shouldn't need to edit anything below this line.
# ==============================================================================

echo "--- Starting Re-Basecalling Run ---"

# --- Define software versions and names ---
DORADO_VERSION="1.0.2" 
MODEL_NAME="dna_r10.4.1_e8.2_400bps_sup@v5.0.0"

# --- Create a dedicated directory for the analysis ---
# The 'eval' command correctly handles the '~' in the path
eval mkdir -p "$OUTPUT_DIR"
eval cd "$OUTPUT_DIR"
echo "[INFO] Created and changed to working directory: $(pwd)"

# --- Step 1: Find and Unpack Dorado ---
DORADO_TAR_FILE="dorado-${DORADO_VERSION}-linux-x64.tar.gz"
DORADO_EXECUTABLE_DIR="dorado-${DORADO_VERSION}-linux-x64"

if [ ! -d "$DORADO_EXECUTABLE_DIR" ]; then
    # Since we've run this before, we expect the file to be in the old directory.
    # We will copy it from there instead of re-downloading.
    OLD_DORADO_PATH="~/dorado_analysis/${DORADO_TAR_FILE}"
    eval cp "$OLD_DORADO_PATH" .
    
    if [ ! -f "$DORADO_TAR_FILE" ]; then
        echo "!!! CRITICAL ERROR: Dorado tarball not found in old location."
        echo "!!! Please manually copy ${DORADO_TAR_FILE} to $(pwd)"
        exit 1
    fi
    echo "[INFO] Dorado tarball found. Unpacking..."
    tar -zxvf "$DORADO_TAR_FILE"
fi

DORADO_BIN="./${DORADO_EXECUTABLE_DIR}/bin/dorado"

# --- Step 2: Get the Model ---
MODEL_PATH="./${MODEL_NAME}"
if [ ! -d "$MODEL_PATH" ]; then
    echo "[INFO] Model not found. Downloading: ${MODEL_NAME}..."
    "$DORADO_BIN" download --model "$MODEL_NAME"
else
    echo "[INFO] Model directory already exists. Skipping download."
fi

# --- Step 3: Run Basecalling Loop ---
echo ""
echo "========================================="
echo "      Starting Basecalling Loop"
echo "========================================="

for sample_name in "${SAMPLES[@]}"; do
    echo ""
    echo "--- Processing Sample: ${sample_name} ---"
    INPUT_DIR="${BASE_INPUT_DIR}/${sample_name}"
    OUTPUT_BAM="${sample_name}_calls.bam" # Naming it to match your original files

    if [ ! -d "$INPUT_DIR" ]; then
        echo "!!! WARNING: Directory for sample ${sample_name} not found at ${INPUT_DIR}. Skipping."
        continue
    fi

    echo "Model:  ${MODEL_PATH}"
    echo "Input:  ${INPUT_DIR} (recursive search)"
    echo "Output: ${OUTPUT_BAM}"
    
    "$DORADO_BIN" basecaller --recursive "$MODEL_PATH" "$INPUT_DIR" > "$OUTPUT_BAM"

    echo "--- Finished processing sample ${sample_name}. ---"
done

echo ""
echo "========================================="
echo "      All Samples Processed!"
echo "========================================="
echo "Pipeline finished successfully."
echo "New BAM files are located in: ${OUTPUT_DIR}"
