#!/usr/bin/env bash
# ======================================================================
# READY-USE TXT FILE FOR UBUNTU
# GSE183437 miRNA-seq pipeline for all SRR accessions listed below
#
# HOW TO RUN
# 1) Put this TXT file in your project folder.
# 2) Put your reference files inside a folder named ref/ in the same folder.
# 3) Run directly as a bash script:
#      bash GSE183437_Ubuntu_Pipeline_READY_USE.txt
#
# Alternative, provide the project path:
#      bash GSE183437_Ubuntu_Pipeline_READY_USE.txt /media/rebaz/ECBEC142BEC105D4/YOUR_PROJECT_FOLDER
#
# Required reference files:
#      ref/hg38.fa
#      ref/hg38.1.ebwt
#      ref/hg38.2.ebwt
#      ref/hg38.3.ebwt
#      ref/hg38.4.ebwt
#      ref/hg38.rev.1.ebwt
#      ref/hg38.rev.2.ebwt
#      ref/mature_hsa.fa
#      ref/hairpin_hsa.fa
# ======================================================================
# Dataset: GSE183437
# Number of SRR runs: 30

set -euo pipefail

# ============================================================
# GSE183437 selected SRR miRNA-seq pipeline
# Prepared from the same Ubuntu/Linux command structure used
# in this project and applied to all GSE183437 accessions embedded in this TXT file
# ============================================================
#
# USAGE
#   1) Put this TXT file and the ref/ folder in one project folder.
#   2) The ref/ folder must contain:
#        hg38.fa
#        hg38.1.ebwt, hg38.2.ebwt, hg38.3.ebwt, hg38.4.ebwt,
#        hg38.rev.1.ebwt, hg38.rev.2.ebwt
#        mature_hsa.fa
#        hairpin_hsa.fa
#   3) Run:
#        chmod +x GSE183437_Ubuntu_Pipeline_READY_USE.txt
#        ./GSE183437_Ubuntu_Pipeline_READY_USE.txt
#
# Optional:
#   ./GSE183437_Ubuntu_Pipeline_READY_USE.txt /path/to/project_folder
#
# Notes:
#   - The step order, folder structure, and tools follow the uploaded
#     Ubuntu Linux Commands file.
#   - The script reads all SRR accessions from selected_srr_GSE183437.txt, which is created automatically if missing.
#   - It is resume-safe: existing outputs are skipped.
#   - It handles both outputs from fasterq-dump:
#       single-end: SRRxxxx.fastq.gz
#       split-files/read1: SRRxxxx_1.fastq.gz
#     If _1 exists, read 1 is used for the miRNA pipeline, as in the
#     original command structure.
# ============================================================

# ---------------------------
# USER SETTINGS
# ---------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${1:-${SCRIPT_DIR}}"
cd "${PROJECT_DIR}"

SRR_LIST="${PROJECT_DIR}/selected_srr_GSE183437.txt"

# ------------------------------------------------------------
# EMBEDDED SRR LIST
# This file is self-contained. If selected_srr_GSE183437.txt is not present
# in the project folder, it will be created automatically from the
# 30 GSE183437 SRR accessions supplied by the user.
# ------------------------------------------------------------
if [[ ! -f "${SRR_LIST}" ]]; then
  cat > "${SRR_LIST}" <<'SRR_LIST_EOF'
SRR15720387
SRR15720388
SRR15720389
SRR15720390
SRR15720391
SRR15720392
SRR15720393
SRR15720394
SRR15720395
SRR15720396
SRR15720397
SRR15720398
SRR15720399
SRR15720400
SRR15720401
SRR15720402
SRR15720403
SRR15720404
SRR15720405
SRR15720406
SRR15720407
SRR15720408
SRR15720409
SRR15720410
SRR15720411
SRR15720412
SRR15720413
SRR15720414
SRR15720415
SRR15720416
SRR_LIST_EOF
fi

# references
REF_DIR="${PROJECT_DIR}/ref"
GENOME_PREFIX="${REF_DIR}/hg38"
GENOME_FASTA="${REF_DIR}/hg38.fa"
MATURE_FASTA="${REF_DIR}/mature_hsa.fa"
HAIRPIN_FASTA="${REF_DIR}/hairpin_hsa.fa"

# threads
FASTQC_THREADS="${FASTQC_THREADS:-8}"
FASTP_THREADS="${FASTP_THREADS:-4}"
FASTERQ_THREADS="${FASTERQ_THREADS:-8}"
MAPPER_THREADS="${MAPPER_THREADS:-4}"

# trimming and mapper settings
FASTP_MIN_LEN="${FASTP_MIN_LEN:-15}"
MAPPER_MIN_LEN="${MAPPER_MIN_LEN:-18}"

# ---------------------------
# FOLDERS
# ---------------------------
mkdir -p \
00_logs \
00_versions \
SRR_files \
FASTQ_files \
FASTQ_files_extra \
01_fastqc_raw \
02_trimmed_fastp \
02_trimmed_fastp_unzipped \
03_fastqc_trimmed \
04_mapper_collapsed \
05_mapper_arf \
06_quantifier \
07_multiqc_raw \
08_multiqc_trimmed \
09_fastp_reports \
11_counts_collected \
mapper_logs

# ---------------------------
# OPTIONAL HELPER R SCRIPT
# This creates the R helper file automatically if it is absent.
# ---------------------------
if [[ ! -f "02_build_miRNA_count_matrix_from_quantifier.R" ]]; then
  cat > "02_build_miRNA_count_matrix_from_quantifier.R" <<'R_SCRIPT_EOF'
#!/usr/bin/env Rscript

# ============================================================
# Optional helper: combine miRDeep2 quantifier miRNA_expressed.csv
# files into one miRNA count matrix.
# Input folder: 11_counts_collected
# Output folder: 12_count_matrix
# ============================================================

suppressPackageStartupMessages({
  if (!requireNamespace("readr", quietly = TRUE)) stop("Please install R package: readr")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install R package: dplyr")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Please install R package: stringr")
})

library(readr)
library(dplyr)
library(stringr)

input_dir <- "11_counts_collected"
output_dir <- "12_count_matrix"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(input_dir, pattern = "__miRNA_expressed\\.csv$", full.names = TRUE)

if (length(files) == 0) {
  stop("No collected miRNA_expressed.csv files found in 11_counts_collected")
}

read_one <- function(f) {
  sample <- basename(f)
  sample <- sub("__miRNA_expressed\\.csv$", "", sample)

  x <- readr::read_csv(f, show_col_types = FALSE)

  # miRDeep2/quantifier output names can vary slightly between versions.
  # This helper attempts to identify the miRNA ID and read-count columns robustly.
  cn <- names(x)
  mir_col <- cn[grepl("miRNA|mature|name|precursor", cn, ignore.case = TRUE)][1]
  count_col <- cn[grepl("read|count|total", cn, ignore.case = TRUE)][1]

  if (is.na(mir_col) || is.na(count_col)) {
    stop(paste("Could not identify miRNA/count columns in:", f,
               "\nColumns were:", paste(cn, collapse = ", ")))
  }

  x %>%
    transmute(miRNA_ID = .data[[mir_col]],
              !!sample := suppressWarnings(as.numeric(.data[[count_col]]))) %>%
    filter(!is.na(miRNA_ID), miRNA_ID != "") %>%
    group_by(miRNA_ID) %>%
    summarise(!!sample := sum(.data[[sample]], na.rm = TRUE), .groups = "drop")
}

mat <- Reduce(function(a, b) full_join(a, b, by = "miRNA_ID"), lapply(files, read_one))
mat[is.na(mat)] <- 0

readr::write_csv(mat, file.path(output_dir, "all_samples_miRBase_counts.csv"))
readr::write_tsv(mat, file.path(output_dir, "all_samples_miRBase_counts.tsv"))

message("Saved:")
message(file.path(output_dir, "all_samples_miRBase_counts.csv"))
message(file.path(output_dir, "all_samples_miRBase_counts.tsv"))
R_SCRIPT_EOF
fi

# ---------------------------
# CAPTURE SOFTWARE VERSIONS
# ---------------------------
{
  echo "=== DATE ==="
  date
  echo

  echo "=== PROJECT DIRECTORY ==="
  echo "${PROJECT_DIR}"
  echo

  echo "=== SYSTEM ==="
  uname -a
  echo
  if command -v lsb_release >/dev/null 2>&1; then
    lsb_release -a 2>/dev/null || true
  fi
  echo

  echo "=== fastqc ==="
  fastqc --version 2>&1 || true
  echo

  echo "=== multiqc ==="
  multiqc --version 2>&1 || true
  echo

  echo "=== fastp ==="
  fastp --version 2>&1 || true
  echo

  echo "=== prefetch ==="
  prefetch --version 2>&1 || true
  echo

  echo "=== fasterq-dump ==="
  fasterq-dump --version 2>&1 || true
  echo

  echo "=== bowtie ==="
  bowtie --version 2>&1 || true
  echo

  echo "=== mapper.pl ==="
  mapper.pl 2>&1 | head -20 || true
  echo

  echo "=== quantifier.pl ==="
  quantifier.pl 2>&1 | head -20 || true
  echo
} > 00_versions/software_versions.txt

# ---------------------------
# CHECK INPUTS
# ---------------------------
if [[ ! -f "${SRR_LIST}" ]]; then
  echo "ERROR: ${SRR_LIST} not found"
  exit 1
fi

for f in \
  "${GENOME_FASTA}" \
  "${MATURE_FASTA}" \
  "${HAIRPIN_FASTA}" \
  "${GENOME_PREFIX}.1.ebwt" \
  "${GENOME_PREFIX}.2.ebwt" \
  "${GENOME_PREFIX}.3.ebwt" \
  "${GENOME_PREFIX}.4.ebwt" \
  "${GENOME_PREFIX}.rev.1.ebwt" \
  "${GENOME_PREFIX}.rev.2.ebwt"
do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing reference file: $f"
    exit 1
  fi
done

# Clean SRR list into a reproducible working copy
grep -v '^[[:space:]]*$' "${SRR_LIST}" | sed 's/\r$//' | sort -u > 00_logs/selected_srr.clean.txt

echo "Project directory: ${PROJECT_DIR}"
echo "SRR list: ${SRR_LIST}"
echo "Number of unique selected SRRs: $(wc -l < 00_logs/selected_srr.clean.txt)"

# ---------------------------
# DOWNLOAD SRA FILES
# ---------------------------
echo "=== Step 1: Downloading SRA files ==="

while read -r srr; do
  [[ -z "${srr}" ]] && continue

  existing_sra="$(find SRR_files -type f -name "${srr}.sra" | head -n 1 || true)"
  if [[ -n "${existing_sra}" && -f "${existing_sra}" ]]; then
    echo "Skipping ${srr}: already downloaded at ${existing_sra}"
    continue
  fi

  echo "Downloading ${srr}"
  prefetch "${srr}" -O SRR_files > "00_logs/${srr}_prefetch.log" 2>&1
done < 00_logs/selected_srr.clean.txt

# ---------------------------
# CONVERT TO FASTQ
# ---------------------------
echo "=== Step 2: Converting SRA to FASTQ ==="

while read -r srr; do
  [[ -z "${srr}" ]] && continue

  if ls "FASTQ_files/${srr}"*.fastq.gz >/dev/null 2>&1; then
    echo "Skipping ${srr}: FASTQ already exists"
    continue
  fi

  sra_file="$(find SRR_files -type f -name "${srr}.sra" | head -n 1 || true)"
  if [[ -z "${sra_file}" || ! -f "${sra_file}" ]]; then
    echo "WARNING: missing SRA file for ${srr}" | tee -a 00_logs/missing_sra.log
    continue
  fi

  echo "Converting ${srr}"
  fasterq-dump "${sra_file}" \
    -O FASTQ_files \
    --split-files \
    --threads "${FASTERQ_THREADS}" \
    > "00_logs/${srr}_fasterq.log" 2>&1

  gzip -f FASTQ_files/${srr}*.fastq
done < 00_logs/selected_srr.clean.txt

# Move single plain FASTQ away only when split read-1 also exists
while read -r srr; do
  [[ -z "${srr}" ]] && continue
  if [[ -f "FASTQ_files/${srr}.fastq.gz" && -f "FASTQ_files/${srr}_1.fastq.gz" ]]; then
    mv "FASTQ_files/${srr}.fastq.gz" FASTQ_files_extra/ 2>/dev/null || true
  fi
done < 00_logs/selected_srr.clean.txt

# Build the FASTQ input list used by downstream steps
# If SRR_1.fastq.gz exists, use read 1; otherwise use SRR.fastq.gz.
: > all_fastq.txt
: > 00_logs/fastq_input_manifest.tsv
echo -e "SRR\tSample_ID\tInput_FASTQ" > 00_logs/fastq_input_manifest.tsv

while read -r srr; do
  [[ -z "${srr}" ]] && continue

  if [[ -f "FASTQ_files/${srr}_1.fastq.gz" ]]; then
    sample="${srr}_1"
    input_fastq="FASTQ_files/${srr}_1.fastq.gz"
  elif [[ -f "FASTQ_files/${srr}.fastq.gz" ]]; then
    sample="${srr}"
    input_fastq="FASTQ_files/${srr}.fastq.gz"
  else
    echo "WARNING: no usable FASTQ found for ${srr}" | tee -a 00_logs/missing_fastq.log
    continue
  fi

  echo "${sample}" >> all_fastq.txt
  echo -e "${srr}\t${sample}\t${input_fastq}" >> 00_logs/fastq_input_manifest.tsv
done < 00_logs/selected_srr.clean.txt

sort -u all_fastq.txt -o all_fastq.txt

# ---------------------------
# RAW FASTQC: ONLY MISSING
# ---------------------------
echo "=== Step 3: Raw FastQC on missing files ==="

ls 01_fastqc_raw/*_fastqc.html 2>/dev/null \
| sed 's#.*/##' \
| sed 's/_fastqc.html//' \
| sort > done_fastqc.txt || true

comm -23 all_fastq.txt done_fastqc.txt > fastqc_missing.txt || true

while read -r sample; do
  [[ -z "${sample}" ]] && continue
  fastqc "FASTQ_files/${sample}.fastq.gz" \
    -o 01_fastqc_raw \
    -t "${FASTQC_THREADS}" \
    >> 00_logs/fastqc_raw_missing.log 2>&1
done < fastqc_missing.txt

multiqc 01_fastqc_raw -o 07_multiqc_raw \
  >> 00_logs/multiqc_raw.log 2>&1 || true

# ---------------------------
# TRIMMING WITH FASTP
# ---------------------------
echo "=== Step 4: fastp trimming ==="

ls 02_trimmed_fastp/*.trimmed.fastq.gz 2>/dev/null \
| sed 's#.*/##' \
| sed 's/.trimmed.fastq.gz//' \
| sort > done_trimmed.txt || true

comm -23 all_fastq.txt done_trimmed.txt > trim_missing.txt || true

while read -r sample; do
  [[ -z "${sample}" ]] && continue

  fastp \
    -i "FASTQ_files/${sample}.fastq.gz" \
    -o "02_trimmed_fastp/${sample}.trimmed.fastq.gz" \
    --thread "${FASTP_THREADS}" \
    --length_required "${FASTP_MIN_LEN}" \
    --cut_front \
    --cut_tail \
    --cut_window_size 4 \
    --cut_mean_quality 20 \
    --json "09_fastp_reports/${sample}.fastp.json" \
    --html "09_fastp_reports/${sample}.fastp.html" \
    > "00_logs/${sample}_fastp.log" 2>&1
done < trim_missing.txt

# ---------------------------
# TRIMMED FASTQC: ONLY MISSING
# ---------------------------
echo "=== Step 5: Trimmed FastQC on missing files ==="

ls 03_fastqc_trimmed/*_fastqc.html 2>/dev/null \
| sed 's#.*/##' \
| sed 's/_fastqc.html//' \
| sort > done_fastqc_trimmed.txt || true

ls 02_trimmed_fastp/*.trimmed.fastq.gz 2>/dev/null \
| sed 's#.*/##' \
| sed 's/.trimmed.fastq.gz//' \
| sort > all_trimmed.txt || true

comm -23 all_trimmed.txt done_fastqc_trimmed.txt > fastqc_trimmed_missing.txt || true

while read -r sample; do
  [[ -z "${sample}" ]] && continue
  fastqc "02_trimmed_fastp/${sample}.trimmed.fastq.gz" \
    -o 03_fastqc_trimmed \
    -t "${FASTQC_THREADS}" \
    >> 00_logs/fastqc_trimmed_missing.log 2>&1
done < fastqc_trimmed_missing.txt

multiqc 03_fastqc_trimmed -o 08_multiqc_trimmed \
  >> 00_logs/multiqc_trimmed.log 2>&1 || true

# ---------------------------
# MIRDEEP2 MAPPER
# unzip one sample at a time to avoid disk-space issues
# ---------------------------
echo "=== Step 6: mapper.pl ==="

rm -rf 02_trimmed_fastp_unzipped
mkdir -p 02_trimmed_fastp_unzipped mapper_logs 04_mapper_collapsed 05_mapper_arf

ls 02_trimmed_fastp/*.trimmed.fastq.gz 2>/dev/null \
| sed 's#.*/##' \
| sed 's/.trimmed.fastq.gz//' \
| sort > all_trimmed.txt || true

while read -r sample; do
  [[ -z "${sample}" ]] && continue

  if [[ -f "04_mapper_collapsed/${sample}_collapsed.fa" && -f "05_mapper_arf/${sample}_vs_hg38.arf" ]]; then
    echo "Skipping ${sample}: mapper outputs already exist"
    continue
  fi

  echo "Processing ${sample} ..."

  gunzip -c "02_trimmed_fastp/${sample}.trimmed.fastq.gz" > "02_trimmed_fastp_unzipped/${sample}.trimmed.fastq" || {
    echo "Unzip failed for ${sample}" >> 00_logs/mapper_failed_samples.log
    rm -f "02_trimmed_fastp_unzipped/${sample}.trimmed.fastq"
    continue
  }

  mapper.pl "02_trimmed_fastp_unzipped/${sample}.trimmed.fastq" \
    -e -h -i -j -m -l "${MAPPER_MIN_LEN}" \
    -p ref/hg38 \
    -s "04_mapper_collapsed/${sample}_collapsed.fa" \
    -t "05_mapper_arf/${sample}_vs_hg38.arf" \
    -v -o "${MAPPER_THREADS}" \
    > "00_logs/${sample}_mapper.log" 2>&1 || {
      echo "mapper.pl failed for ${sample}" >> 00_logs/mapper_failed_samples.log
      rm -f "02_trimmed_fastp_unzipped/${sample}.trimmed.fastq"
      continue
    }

  rm -f "02_trimmed_fastp_unzipped/${sample}.trimmed.fastq"
done < all_trimmed.txt

# ---------------------------
# MIRDEEP2 QUANTIFIER
# ---------------------------
echo "=== Step 7: quantifier.pl ==="

find 06_quantifier -mindepth 1 -maxdepth 1 -type d 2>/dev/null \
| sed 's#.*/##' \
| sort > done_quantifier.txt || true

ls 04_mapper_collapsed/*_collapsed.fa 2>/dev/null \
| sed 's#.*/##' \
| sed 's/_collapsed.fa//' \
| sort > all_mapper_samples.txt || true

comm -23 all_mapper_samples.txt done_quantifier.txt > quantifier_missing.txt || true

while read -r sample; do
  [[ -z "${sample}" ]] && continue

  mkdir -p "06_quantifier/${sample}"

  quantifier.pl \
    -p ref/hairpin_hsa.fa \
    -m ref/mature_hsa.fa \
    -r "04_mapper_collapsed/${sample}_collapsed.fa" \
    -t hsa \
    -y "${sample}" \
    > "00_logs/${sample}_quantifier.log" 2>&1 || {
      echo "quantifier.pl failed for ${sample}" >> 00_logs/quantifier_failed_samples.log
      continue
    }

  mv expression_analyses "06_quantifier/${sample}/" 2>/dev/null || true
  mv miRBase.mrd "06_quantifier/${sample}/" 2>/dev/null || true
  mv *_mapped.arf "06_quantifier/${sample}/" 2>/dev/null || true
  mv *_not_mapped.fa "06_quantifier/${sample}/" 2>/dev/null || true
  mv *_converted.fa "06_quantifier/${sample}/" 2>/dev/null || true
done < quantifier_missing.txt

# ---------------------------
# COLLECT COUNT FILES
# ---------------------------
echo "=== Step 8: Collect miRNA_expressed.csv files ==="

mkdir -p 11_counts_collected

find 06_quantifier -type f -name "miRNA_expressed.csv" | while read -r f; do
  sample=$(echo "$f" | cut -d/ -f2)
  cp "$f" "11_counts_collected/${sample}__miRNA_expressed.csv"
done

# ---------------------------
# OPTIONAL: BUILD A SIMPLE COUNT MATRIX FROM COLLECTED FILES
# ---------------------------
# This is kept separate from the main pipeline logic. It will run only if Rscript is available
# and the helper file exists in the project folder.
if command -v Rscript >/dev/null 2>&1 && [[ -f "02_build_miRNA_count_matrix_from_quantifier.R" ]]; then
  echo "=== Optional Step: Build count matrix from collected miRNA_expressed.csv files ==="
  Rscript 02_build_miRNA_count_matrix_from_quantifier.R \
    > 00_logs/build_count_matrix.log 2>&1 || true
fi

# ---------------------------
# FINAL CHECKS
# ---------------------------
{
  echo "=== FINAL SUMMARY ==="
  echo "Selected SRRs:"
  wc -l < 00_logs/selected_srr.clean.txt
  echo "SRA files:"
  find SRR_files -type f -name "*.sra" | wc -l
  echo "FASTQ input files used:"
  wc -l < all_fastq.txt
  echo "Trimmed FASTQ files:"
  ls 02_trimmed_fastp/*.trimmed.fastq.gz 2>/dev/null | wc -l
  echo "Mapper collapsed files:"
  ls 04_mapper_collapsed/*_collapsed.fa 2>/dev/null | wc -l
  echo "Mapper ARF files:"
  ls 05_mapper_arf/*.arf 2>/dev/null | wc -l
  echo "Quantifier sample dirs:"
  find 06_quantifier -mindepth 1 -maxdepth 1 -type d | wc -l
  echo "Collected miRNA_expressed.csv files:"
  ls 11_counts_collected/*.csv 2>/dev/null | wc -l
  echo
  echo "Failed/missing logs, if present:"
  ls 00_logs/*missing* 00_logs/*failed* 2>/dev/null || true
  echo
  echo "Disk usage:"
  df -h .
} > 00_logs/final_summary.txt

echo "Pipeline finished."
echo "Version log: ${PROJECT_DIR}/00_versions/software_versions.txt"
echo "FASTQ manifest: ${PROJECT_DIR}/00_logs/fastq_input_manifest.tsv"
echo "Final summary: ${PROJECT_DIR}/00_logs/final_summary.txt"
