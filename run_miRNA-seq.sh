#!/bin/bash

# --- Stricter Shell Safety ---
# Exit immediately if a command exits with a non-zero status.
set -e
# Fail a pipeline if any command fails, not just the last one.
set -o pipefail
# Set a sane default for the Internal Field Separator.
IFS=$'\n\t'

# Author: ZHANG HE
#
# Description:
#   This script executes a complete miRNA-seq workflow with a single command.
#   It handles raw FASTQ processing, miRNA quantification with miRDeep2,
#   and differential expression with DESeq2.

# --- 1. Default Parameters & Help Function ---
THREADS=8
REF_LEVEL=""
# Set a default for the SIF path, can be overridden by -c
SIF_PATH="$(dirname "$0")/miRNA.sif"

HELP_MSG="Usage: $0 -s <sample_sheet.csv> -o <out_dir> -r <ref_dir> -c <container.sif> [OPTIONS]

Required:
  -s  Sample sheet (CSV). Format for single-end reads: sample,condition,fastq_path
  -o  Main output directory.
  -r  Reference data directory.
  -c  Path to the Singularity container (miRNA.sif).

Optional:
  -t  Number of threads to use (Default: 8).
  -L  Reference level for DESeq2 comparison (e.g., 'Control'). If unset, DESeq2 uses alphabetical order.
  -h  Display this help message.
"

# --- 2. Parse Command-line Arguments ---
while getopts "s:o:r:c:t:L:h" opt; do
  case ${opt} in
    s ) SAMPLE_SHEET=$(realpath "${OPTARG}") ;;
    o ) OUT_DIR=$(realpath "${OPTARG}") ;;
    r ) REF_DIR=$(realpath "${OPTARG}") ;;
    c ) SIF_PATH=$(realpath "${OPTARG}") ;;
    t ) THREADS=${OPTARG} ;;
    L ) REF_LEVEL=${OPTARG} ;;
    h ) echo "${HELP_MSG}"; exit 0 ;;
    \? ) echo "Invalid option: -${OPTARG}" >&2; echo "${HELP_MSG}"; exit 1 ;;
  esac
done

# Check for mandatory arguments
if [ -z "${SAMPLE_SHEET}" ] || [ -z "${OUT_DIR}" ] || [ -z "${REF_DIR}" ] || [ -z "${SIF_PATH}" ]; then
    echo "Error: Missing mandatory arguments." >&2; echo "${HELP_MSG}"; exit 1
fi

# --- 3. Setup Environment and Directories ---
if [ ! -f "${SIF_PATH}" ]; then
    echo "Error: Singularity container not found at: ${SIF_PATH}" >&2; exit 1
fi
if [ ! -d "${REF_DIR}" ]; then
    echo "Error: Reference directory not found at: ${REF_DIR}" >&2; exit 1
fi

mkdir -p "${OUT_DIR}"

# Define the main Singularity execution command template
SINGULARITY_BASE_CMD="singularity exec --cleanenv -B ${OUT_DIR}:/output -B ${REF_DIR}:/reference"

echo "================================================="
echo "====== miRNA-seq Pipeline Started ======="
echo "================================================="
echo "Sample Sheet: ${SAMPLE_SHEET}"
echo "Output Directory: ${OUT_DIR}"
echo "Reference Directory: ${REF_DIR}"
echo "Threads: ${THREADS}"
echo "Singularity Container: ${SIF_PATH}"
if [ -n "${REF_LEVEL}" ]; then echo "DESeq2 Ref Level: ${REF_LEVEL}"; fi
echo "================================================="

# --- 4. Loop Through and Process Each Sample ---
# Use process substitution and 'tr' to robustly handle CSV files with CRLF line endings
while IFS=, read -r SAMPLE_NAME CONDITION FQ1; do
    # Trim whitespace from all fields read from the sample sheet
    SAMPLE_NAME=$(echo "${SAMPLE_NAME}" | tr -d '[:space:]')
    CONDITION=$(echo "${CONDITION}" | tr -d '[:space:]')
    FQ1=$(echo "${FQ1}" | tr -d '[:space:]')

    echo -e "\n\n======== Starting to process sample: [${SAMPLE_NAME}] ========"

    if [ ! -f "${FQ1}" ]; then
        echo "Error: FASTQ file for sample ${SAMPLE_NAME} not found: ${FQ1}. Skipping." >&2
        continue
    fi

    # --- 4.1. Set up directories and dynamic binds ---
    FQ_DIR=$(dirname "${FQ1}")
    SAMPLE_OUT_DIR="/output/${SAMPLE_NAME}" # Path inside the container
    mkdir -p "${OUT_DIR}/${SAMPLE_NAME}"   # Create directory on the host

    SINGULARITY_CMD="singularity exec --cleanenv -B ${OUT_DIR}:/output -B ${REF_DIR}:/reference -B ${FQ_DIR}:/reads ${SIF_PATH}"

    # --- 4.2. Per-Sample Processing Pipeline ---
    echo "======== [${SAMPLE_NAME}] Step 1: Raw FastQC ========"
    eval "$SINGULARITY_CMD fastqc -t ${THREADS} /reads/$(basename ${FQ1}) -o ${SAMPLE_OUT_DIR}"

    echo "======== [${SAMPLE_NAME}] Step 2: Adapter Trimming (Trim Galore) ========"
    TRIMMED_FQ="${SAMPLE_OUT_DIR}/${SAMPLE_NAME}_trimmed.fq.gz"
    # Adapter sequence is now hardcoded as it is essential for this step
    eval "$SINGULARITY_CMD trim_galore \
        --fastqc \
        --length 18 --max_length 35 \
        -o ${SAMPLE_OUT_DIR} \
        /reads/$(basename ${FQ1})"
    # Find the trimmed file name as Trim Galore appends _trimmed.fq.gz
    TRIMMED_FQ_RAW=$(find "${OUT_DIR}/${SAMPLE_NAME}" -name "*_trimmed.fq.gz")
    mv "${TRIMMED_FQ_RAW}" "${OUT_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.fq.gz"
    TRIMMED_FQ_IN_CONTAINER="${SAMPLE_OUT_DIR}/${SAMPLE_NAME}.fq.gz"

    echo "======== [${SAMPLE_NAME}] Step 3: Convert to FASTA and Clean Headers ========"
    TRIMMED_FA="${SAMPLE_OUT_DIR}/${SAMPLE_NAME}.fa"
    CLEANED_FA="${SAMPLE_OUT_DIR}/${SAMPLE_NAME}_clean.fa"

    # Step 3a: Convert FASTQ to FASTA
    eval "$SINGULARITY_CMD seqkit fq2fa "${TRIMMED_FQ_IN_CONTAINER}" -o "${TRIMMED_FA}""

    # Step 3b: Clean FASTA headers using the definitive, correct method
    eval "$SINGULARITY_CMD seqkit replace -p '([^\s\/]+).*' -r '{1}' "${TRIMMED_FA}" -o "${CLEANED_FA}""

    echo "======== [${SAMPLE_NAME}] Step 4: Map reads to genome (miRDeep2 Mapper) ========"
    MAPPED_FA="${SAMPLE_OUT_DIR}/${SAMPLE_NAME}_mapped.fa"
    MAPPED_ARF="${SAMPLE_OUT_DIR}/${SAMPLE_NAME}.arf"
    
    # Step 4: Use the correctly created CLEANED_FA file for mapping
    # The redundant adapter parameter '-k' has been removed.
    eval "$SINGULARITY_CMD mapper.pl "${CLEANED_FA}" \
        -c -i -j \
        -l 18 -m \
        -p /reference/hsa \
        -s ${MAPPED_FA} \
        -t ${MAPPED_ARF} -v"

    echo "======== [${SAMPLE_NAME}] Step 5: Quantify known miRNAs (miRDeep2 Quantifier) ========"
    eval "$SINGULARITY_CMD quantifier.pl \
        -p /reference/hairpin.fa \
        -m /reference/mature.fa \
        -r ${MAPPED_FA} \
        -t hsa \
        -y ${SAMPLE_NAME} \
        -d "
    # Move quantifier output to the sample directory
    mv miRNAs_expressed*.csv "${OUT_DIR}/${SAMPLE_NAME}/miRNAs_expressed.csv"
    mv bowtie.log "${OUT_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_bowtie.log"
    rm expression*.html
    rm -r ./expression_analyses
    rm -r ./mapper_logs
    echo "======== [${SAMPLE_NAME}] Processing finished. ========"
done < <(tail -n +2 "${SAMPLE_SHEET}" | tr -d '\r')

echo -e "\n\n==============================================="
echo "======= All samples processed. Starting final analysis... ======="
echo "==============================================="

# --- 5. Merge miRNA counts and Run DESeq2 ---
echo "======== Step 6: Merging counts and Running Differential Expression Analysis (DESeq2) ========"
MERGED_COUNTS_FILE="${OUT_DIR}/merged_mirna_counts.csv"
DESEQ_R_SCRIPT="${OUT_DIR}/deseq_analysis_runner.R"
DESEQ_SAMPLESHEET="${OUT_DIR}/deseq_samplesheet.csv"
OUTPUT_PREFIX="/output/final_results"

# Create header for merged counts
echo -n "miRNA" > "${MERGED_COUNTS_FILE}"
while IFS=, read -r SAMPLE_NAME _; do
    SAMPLE_NAME=$(echo "${SAMPLE_NAME}" | tr -d '[:space:]')
    echo -n ",${SAMPLE_NAME}" >> "${MERGED_COUNTS_FILE}"
done < <(tail -n +2 "${SAMPLE_SHEET}" | tr -d '\r')
echo "" >> "${MERGED_COUNTS_FILE}"

# Use awk to parse and merge quantifier results
awk '
    FNR==1 { next }
    FNR==NR {
        counts[$1] = $5;
        mirnas[$1] = 1;
        next;
    }
    {
        for (mirna in mirnas) {
            if ($1 == mirna) {
                counts[mirna] = counts[mirna] "," $5;
            }
        }
    }
    END {
        for (mirna in mirnas) {
            print mirna "," counts[mirna];
        }
    }
' $(find ${OUT_DIR} -name "miRNAs_expressed.csv") >> "${MERGED_COUNTS_FILE}"


# Create the simplified sample sheet for DESeq2
echo "sample,condition" > "${DESEQ_SAMPLESHEET}"
while IFS=, read -r SAMPLE_NAME CONDITION _; do
    SAMPLE_NAME=$(echo "${SAMPLE_NAME}" | tr -d '[:space:]')
    CONDITION=$(echo "${CONDITION}" | tr -d '[:space:]')
    echo "${SAMPLE_NAME},${CONDITION}" >> "${DESEQ_SAMPLESHEET}"
done < <(tail -n +2 "${SAMPLE_SHEET}" | tr -d '\r')


# Create the R script for DESeq2
cat > "${DESEQ_R_SCRIPT}" << 'EOF'
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript deseq_runner.R <counts_file> <sample_sheet> <output_prefix> [ref_level]", call. = FALSE)
}

COUNTS_PATH <- args[1]
SAMPLE_SHEET_PATH <- args[2]
OUTPUT_PREFIX <- args[3]
REF_LEVEL <- if (length(args) >= 4) args[4] else ""


message("Reading merged counts from: ", COUNTS_PATH)
temp_df <- read.csv(COUNTS_PATH, header = TRUE, check.names = FALSE)

if (any(duplicated(temp_df[, 1]))) {
    message("Found duplicate miRNA names. Aggregating rows by summing counts...")
    grouping_col_name <- names(temp_df)[1]
    count_df <- aggregate(as.formula(paste(". ~", grouping_col_name)), data = temp_df, FUN = sum)
    rownames(count_df) <- count_df[, 1]
    count_df <- count_df[, -1]
} else {
    message("No duplicate miRNA names found. Proceeding directly.")
    count_df <- temp_df
    rownames(count_df) <- count_df[, 1]
    count_df <- count_df[, -1]
}

count_data <- as.matrix(round(count_df))

message("Reading sample sheet from: ", SAMPLE_SHEET_PATH)
sample_table <- read.csv(SAMPLE_SHEET_PATH, header = TRUE, row.names = 1)

if (!all(colnames(count_data) %in% rownames(sample_table))) {
  stop("Column names in counts file do not all match row names in sample sheet.")
}
count_data <- count_data[, rownames(sample_table)]

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_table,
                              design = ~ condition)

if (REF_LEVEL != "" && REF_LEVEL %in% dds$condition) {
    message("Setting reference level for condition to: ", REF_LEVEL)
    dds$condition <- relevel(dds$condition, ref = REF_LEVEL)
} else if (REF_LEVEL != "") {
    warning("Provided reference level '", REF_LEVEL, "' not found in conditions. Using alphabetical default.")
}

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

message("Running DESeq2 analysis...")
dds <- DESeq(dds)
res <- results(dds)
res_df <- as.data.frame(res)

res_df <- res_df[order(res_df$padj, na.last = TRUE),]
res_df <- res_df[grepl("-", rownames(res_df)), ]
significant_mirnas <- subset(res_df, padj < 0.05)
write.table(rownames(significant_mirnas),
            file = paste0(OUTPUT_PREFIX, "_significant_mirnas.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

message("Writing results to output files prefixed with: ", OUTPUT_PREFIX)
write.table(res_df,
            file = paste0(OUTPUT_PREFIX, "_DEG_results.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE) # Set row.names=TRUE to include miRNA names

message("Analysis complete.")
EOF

# Execute the DESeq2 R script
eval "$SINGULARITY_BASE_CMD ${SIF_PATH} Rscript /output/deseq_analysis_runner.R \
    /output/merged_mirna_counts.csv \
    /output/deseq_samplesheet.csv \
    ${OUTPUT_PREFIX} \
    ${REF_LEVEL}"

mv "${OUT_DIR}/final_results_DEG_results.txt" "${OUT_DIR}/deg_results.txt"
SIG_MIRNAS="${OUT_DIR}/final_results_significant_mirnas.txt"
SIG_MIRNAS_FA="${OUT_DIR}/final_results_significant_mirnas.fa"
echo "Differential expression results saved to: ${OUT_DIR}/deg_results.txt"
eval "$SINGULARITY_BASE_CMD ${SIF_PATH} seqkit grep -f ${SIG_MIRNAS} /reference/mature.fa -o ${SIG_MIRNAS_FA}"

# ==============================================================================
# --- 5b. Map Significant miRNAs to Target Genes ---
# ==============================================================================
echo "======== Mapping Significant miRNAs to Families and Target Genes ========"

# --- Define Paths ---
MIR_FAMILY_INFO="${REF_DIR}/miR_Family_Info.txt"
TARGET_INFO="${REF_DIR}/Predicted_Targets_Info.default_predictions.txt"
FINAL_TARGET_FILE="${OUT_DIR}/miRNA_targets_genes.txt"

# --- Check for necessary reference files ---
if [ ! -f "${SIG_MIRNAS}" ] || [ ! -s "${SIG_MIRNAS}" ]; then
    echo "Warning: Significant miRNA file is missing or empty. Skipping target mapping."
elif [ ! -f "${MIR_FAMILY_INFO}" ] || [ ! -f "${TARGET_INFO}" ]; then
    echo "Warning: Skipping miRNA target mapping." >&2
    echo "Reason: Missing reference files: miR_Family_Info.txt or Predicted_Targets_Info.default_predictions.txt in ${REF_DIR}" >&2
else
    # --- Create the final miRNA -> Family -> Gene mapping file ---
    # We use a series of joins to create the final three-column file.
    
    # 1. Prepare temporary sorted files for joining
    sort -u "${SIG_MIRNAS}" > "${OUT_DIR}/tmp.sig_mirnas.sorted"
    # File 2: MiRBase ID (col 4) and miR Family (col 1)
    awk -F'\t' 'NR > 1 {print $4, $1}' "${MIR_FAMILY_INFO}" | sort -k1,1 -u > "${OUT_DIR}/tmp.family_info.sorted"
    # File 3: miR Family (col 1) and Gene Symbol (col 3)
    awk -F'\t' 'NR > 1 {print $1, $3}' "${TARGET_INFO}" | sort -k1,1 -u > "${OUT_DIR}/tmp.target_info.sorted"

    # 2. Join significant miRNAs with their families
    # Output: MiRBase_ID miR_Family
    join -1 1 -2 1 "${OUT_DIR}/tmp.sig_mirnas.sorted" "${OUT_DIR}/tmp.family_info.sorted" \
        | sort -k2,2 > "${OUT_DIR}/tmp.mir_to_family.sorted" # Sort by miR_Family for the next join

    # 3. Join the result with target genes
    # Input 1: MiRBase_ID miR_Family (sorted by miR_Family)
    # Input 2: miR_Family Gene_Symbol (sorted by miR_Family)
    # Join on miR_Family. Output format: miR_Family MiRBase_ID Gene_Symbol
    join -1 2 -2 1 "${OUT_DIR}/tmp.mir_to_family.sorted" "${OUT_DIR}/tmp.target_info.sorted" \
        | awk 'BEGIN{OFS="\t"; print "MiRBase_ID\tmiR_Family\tGene_Symbol"} {print $2, $1, $3}' \
        > "${FINAL_TARGET_FILE}"

    echo "Successfully created miRNA to target gene mapping file."
    echo "Output available at: ${FINAL_TARGET_FILE}"

    # 4. Clean up temporary files
    rm -f "${OUT_DIR}/tmp.sig_mirnas.sorted" \
          "${OUT_DIR}/tmp.family_info.sorted" \
          "${OUT_DIR}/tmp.target_info.sorted" \
          "${OUT_DIR}/tmp.mir_to_family.sorted"
fi

# --- 6. Generate MultiQC Report ---
echo "======== Generating Final MultiQC Report ========"
mkdir -p "${OUT_DIR}/multiqc_report"
eval "$SINGULARITY_BASE_CMD ${SIF_PATH} multiqc /output -o /output/multiqc_report --force"
rm -f "${OUT_DIR}/deseq_samplesheet.csv"
rm -f "${OUT_DIR}/deseq_analysis_runner.R"
while IFS=, read -r SAMPLE_NAME CONDITION FQ1 FQ2; do
    SAMPLE_NAME=$(echo "${SAMPLE_NAME}" | tr -d '[:space:]')
    find "${OUT_DIR}/${SAMPLE_NAME}" -type f \
        ! -name "*mapped.fa" \
        ! -name "*.arf" \
        ! -name "*.csv" \
        -delete
done < <(tail -n +2 "${SAMPLE_SHEET}" | tr -d '\r')
echo -e "\n\n==============================================="
echo "====== miRNA-seq Pipeline Completed Successfully! ======"
echo "==============================================="
echo "All results are located in: ${OUT_DIR}"
echo "View the interactive QC report at: ${OUT_DIR}/multiqc_report/multiqc_report.html"
echo "==============================================="