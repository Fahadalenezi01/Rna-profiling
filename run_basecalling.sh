#!/bin/bash
set -e # This makes the script stop immediately if any command fails.

# ==============================================================================
#                           CONFIGURATION
# ==============================================================================
BASE_INPUT_DIR="/mnt/d/Fahad"
SAMPLES=("50SL50BL" "70SL30BL" "30SL70BL")

# This matches the version of the file you downloaded.
DORADO_VERSION="1.0.2" 

# ==============================================================================
#                           SCRIPT
# ==============================================================================
echo "--- Starting Automated Dorado Basecalling for Multiple Samples ---"

WORKDIR="dorado_analysis"
mkdir -p "$WORKDIR"
cd "$WORKDIR"
echo "[INFO] Created and changed to working directory: $(pwd)"

# --- Step 1: Find and Unpack Dorado ---
DORADO_TAR_FILE="dorado-${DORADO_VERSION}-linux-x64.tar.gz"
DORADO_EXECUTABLE_DIR="dorado-${DORADO_VERSION}-linux-x64"

if [ ! -d "$DORADO_EXECUTABLE_DIR" ]; then
    if [ ! -f "$DORADO_TAR_FILE" ]; then
        echo "!!! CRITICAL ERROR: Dorado tarball not found: ${DORADO_TAR_FILE}"
        exit 1
    fi
    echo "[INFO] Dorado tarball found. Unpacking..."
    tar -zxvf "$DORADO_TAR_FILE"
fi

DORADO_BIN="./${DORADO_EXECUTABLE_DIR}/bin/dorado"

# --- Step 2: Get the Model ---
# FINAL CORRECTION: This is the correct model name for this version of Dorado.
MODEL_NAME="dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
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
    OUTPUT_BAM="${sample_name}.bam"

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
echo "Script finished successfully."
