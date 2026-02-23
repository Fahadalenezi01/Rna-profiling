#!/usr/bin/env bash
set -o errexit   # exit on any command failure
set -o nounset   # treat unset variables as errors
set -o pipefail  # pipelines fail on the first non-zero exit

# Uncomment to capture a full log:
# exec > >(tee "${OUTPUT_DIR}/pipeline.$(date +%Y%m%d_%H%M).log") 2>&1

# ==============================================================================
# lncRNA Identification and Quantification Pipeline
# ==============================================================================

# --- User Configuration ---
GTF_DIR_D_SAMPLES="/mnt/d"
GTF_DIR_D_LCRNA="/mnt/d/lncRNA"
GTF_DIR_HOME="/home/alajdal/stringtie_quantification_results"

BAM_DIR_D="/mnt/d/lncRNA"
BAM_DIR_HOME="/home/alajdal/analysis_pipeline_results"

OUTPUT_DIR="/mnt/d/lncRNA/lncRNA_pipeline_output_full"

REFERENCE_GTF="/mnt/d/gencode.v47.primary_assembly.basic.annotation.gtf"
REFERENCE_FASTA="/mnt/d/GRCh38.primary_assembly.genome.fa"
KNOWN_LNCRNA_GTF="/mnt/d/gencode.v47.long_noncoding_RNAs.gtf"
KNOWN_CODING_GTF_FILTERED="${OUTPUT_DIR}/reference_protein_coding.gtf"

CPU_THREADS=22

SAMPLES=(
    "SL1" "SL2" "SL3" "SL3_rep" "SL4"
    "BL1" "BL2" "BL3" "BL4"
    "50SL50BL" "70SL30BL" "30SL70BL"
)

# --- Prep directories ---
mkdir -p "${OUTPUT_DIR}/requantified_gtfs" "${OUTPUT_DIR}/lncRNA_identification" "${OUTPUT_DIR}/feelnc_output"

echo "--- Initializing Pipeline ---"

# --- STEP 0: Filter for Protein‑Coding Transcripts ---
echo "Filtering for protein‑coding transcripts..."
if [[ ! -f "${REFERENCE_GTF}" ]]; then
echo "ERROR: Reference GTF not found at ${REFERENCE_GTF}. Exiting." >&2
exit 1
fi
awk '$3=="transcript" && /protein_coding/' \
    "${REFERENCE_GTF}" > "${KNOWN_CODING_GTF_FILTERED}"
echo "Filtering for protein‑coding exons/CDS..."
if [[ ! -f "${REFERENCE_GTF}" ]]; then
    echo "ERROR: Reference GTF not found at ${REFERENCE_GTF}. Exiting." >&2
    exit 1
fi
# keep exon and CDS entries for protein_coding genes
awk '$3 ~ /exon|CDS/ && /protein_coding/' \
    "${REFERENCE_GTF}" > "${KNOWN_CODING_GTF_FILTERED}"
 echo "→ Filtered GTF: ${KNOWN_CODING_GTF_FILTERED}"

# ==============================================================================
# STEP 1: Merge Transcripts from all Samples
# ==============================================================================
echo -e "\n--- STEP 1: Merging transcripts ---"
MERGE_LIST_FILE="${OUTPUT_DIR}/gtf_merge_list.txt"
: > "${MERGE_LIST_FILE}"

for sample in "${SAMPLES[@]}"; do
    if [[ -f "${GTF_DIR_D_SAMPLES}/${sample}/stringtie_output_${sample}.gtf" ]]; then
        echo "${GTF_DIR_D_SAMPLES}/${sample}/stringtie_output_${sample}.gtf" >> "${MERGE_LIST_FILE}"
    elif [[ -f "${GTF_DIR_D_LCRNA}/stringtie_output_${sample}.gtf" ]]; then
        echo "${GTF_DIR_D_LCRNA}/stringtie_output_${sample}.gtf" >> "${MERGE_LIST_FILE}"
    elif [[ -f "${GTF_DIR_HOME}/${sample}.transcripts.gtf" ]]; then
        echo "${GTF_DIR_HOME}/${sample}.transcripts.gtf" >> "${MERGE_LIST_FILE}"
    else
        echo "WARNING: Missing GTF for ${sample}, skipping." >&2
    fi
done

if [[ ! -s "${MERGE_LIST_FILE}" ]]; then
    echo "ERROR: No GTF files to merge in ${MERGE_LIST_FILE}. Exiting." >&2
    exit 1
fi

echo "Merging the following GTFs:"
cat "${MERGE_LIST_FILE}"
stringtie --merge -p "${CPU_THREADS}" \
    -G "${REFERENCE_GTF}" \
    -o "${OUTPUT_DIR}/merged_transcripts.gtf" \
    "${MERGE_LIST_FILE}"
echo "→ Merged GTF: ${OUTPUT_DIR}/merged_transcripts.gtf"

# ==============================================================================
# STEP 2: Re‑quantify Expression for All Samples
# ==============================================================================
echo -e "\n--- STEP 2: Re‑quantifying expression ---"
for sample in "${SAMPLES[@]}"; do
    echo "Processing: ${sample}"
    BAM_FILE=""
    OUTPUT_GTF="${OUTPUT_DIR}/requantified_gtfs/${sample}.gtf"

    # potential BAM locations
    for path in \
        "${BAM_DIR_D}/filtered_aligned_${sample}.sorted.bam" \
        "${BAM_DIR_D}/aligned_sorted_${sample}.sorted.bam" \
        "${BAM_DIR_HOME}/${sample}/${sample}.aligned.sorted.bam" \
        "${BAM_DIR_HOME}/${sample}.bam"
    do
        if [[ -f "${path}" ]]; then
            BAM_FILE="${path}"
            break
        fi
    done

    if [[ -z "${BAM_FILE}" ]]; then
        echo "ERROR: No BAM for ${sample}, skipping." >&2
        continue
    fi

    echo "→ Found BAM: ${BAM_FILE}"
    stringtie -e -B -p "${CPU_THREADS}" \
        -G "${OUTPUT_DIR}/merged_transcripts.gtf" \
        -o "${OUTPUT_GTF}" \
        "${BAM_FILE}"
done
echo "Re‑quantification complete."

# ==============================================================================
# STEP 3: Generate Count Matrices
# ==============================================================================
echo -e "\n--- STEP 3: Generating count matrices ---"
PREPDE_INPUT_LIST="${OUTPUT_DIR}/prepDE_input_list.txt"
: > "${PREPDE_INPUT_LIST}"

for sample in "${SAMPLES[@]}"; do
    gtf="${OUTPUT_DIR}/requantified_gtfs/${sample}.gtf"
    if [[ -f "${gtf}" ]]; then
        echo "${sample} ${gtf}" >> "${PREPDE_INPUT_LIST}"
    fi
done

prepDE.py -i "${PREPDE_INPUT_LIST}" \
    -g "${OUTPUT_DIR}/merged_gene_counts.csv" \
    -t "${OUTPUT_DIR}/merged_transcript_counts.csv"
echo "Count matrices saved to ${OUTPUT_DIR}."

# ==============================================================================
# STEP 4: Identify Novel lncRNAs with FEELnc
# ==============================================================================
 # ==============================================================================
 # STEP 4: Identify Novel lncRNAs with FEELnc
 # ==============================================================================
echo -e "\n--- STEP 4: Identifying novel lncRNAs ---"
FEELNC_OUT="${OUTPUT_DIR}/feelnc_output"

echo "Running FEELnc_filter.pl..."
FEELnc_filter.pl \
    -i "${OUTPUT_DIR}/merged_transcripts.gtf" \
    -a "${KNOWN_CODING_GTF_FILTERED}" \
    > "${FEELNC_OUT}/candidate.gtf"
echo -e "\n--- STEP 4: Identifying novel lncRNAs ---"
FEELNC_OUT="${OUTPUT_DIR}/feelnc_output"

echo "Running FEELnc_filter.pl on exons/CDS‐filtered GTF..."
FEELnc_filter.pl \
    -i "${OUTPUT_DIR}/merged_transcripts.gtf" \
    -a "${KNOWN_CODING_GTF_FILTERED}" \
    > "${FEELNC_OUT}/candidate.gtf"

 echo "Running FEELnc_codpot.pl..."
 FEELnc_codpot.pl \
     -i "${FEELNC_OUT}/candidate.gtf" \
     -a "${KNOWN_CODING_GTF_FILTERED}" \
     -g "${REFERENCE_FASTA}" \
     -l "${KNOWN_LNCRNA_GTF}" \
     > "${FEELNC_OUT}/feelnc_codpot_output.txt"
 echo "FEELnc analysis complete."

# ==============================================================================
# STEP 5: Finalize lncRNA Transcript IDs
# ==============================================================================
echo -e "\n--- STEP 5: Finalizing lncRNA list ---"
KNOWN_IDS="${OUTPUT_DIR}/known_lncrna_ids.txt"
NOVEL_IDS="${OUTPUT_DIR}/novel_lncrna_ids.txt"
FINAL_IDS="${OUTPUT_DIR}/final_lncrna_transcript_ids.txt"

# extract known lncRNA transcript_id attributes
awk '$3=="transcript" {
    match($0, /transcript_id "([^"]+)"/, a);
    print a[1]
}' "${KNOWN_LNCRNA_GTF}" > "${KNOWN_IDS}"

# extract novel non‑coding IDs from FEELnc codpot output
awk -F $'\t' '$2=="noncoding" {print $1}' \
    "${FEELNC_OUT}/feelnc_codpot_output.txt" \
    > "${NOVEL_IDS}"

# combine, sort, uniq
cat "${KNOWN_IDS}" "${NOVEL_IDS}" | sort | uniq > "${FINAL_IDS}"

echo "Final lncRNA transcript IDs in ${FINAL_IDS}"
echo -e "\n--- Pipeline Finished Successfully ---"
