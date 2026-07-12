#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: 01_prepare_ready_vcfs.sh --input-dir DIR --output-dir DIR --newname-file FILE

Prepare per-sample CSP-ready VCFs:
1) keep PASS biallelic SNPs
2) rename contigs
3) normalize sample IDs (underscores -> hyphens)
4) index final VCFs
EOF
}

INPUT_DIR=""
OUTPUT_DIR=""
NEWNAME_FILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir) INPUT_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --newname-file) NEWNAME_FILE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$INPUT_DIR" ]] || { echo "--input-dir is required" >&2; exit 1; }
[[ -n "$OUTPUT_DIR" ]] || { echo "--output-dir is required" >&2; exit 1; }
[[ -n "$NEWNAME_FILE" ]] || { echo "--newname-file is required" >&2; exit 1; }

command -v bcftools >/dev/null 2>&1 || { echo "bcftools not found in PATH" >&2; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "tabix not found in PATH" >&2; exit 1; }

[[ -d "$INPUT_DIR" ]] || { echo "Input directory not found: $INPUT_DIR" >&2; exit 1; }
[[ -f "$NEWNAME_FILE" ]] || { echo "Rename file not found: $NEWNAME_FILE" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR"
log_file="${OUTPUT_DIR}/prepare_ready_vcfs.log"

shopt -s nullglob
vcf_files=("${INPUT_DIR}"/*.vcf.gz)
shopt -u nullglob

if [[ ${#vcf_files[@]} -eq 0 ]]; then
  echo "No .vcf.gz files found in ${INPUT_DIR}" | tee -a "$log_file"
  exit 1
fi

processed=0
for vcf in "${vcf_files[@]}"; do
  base="$(basename "$vcf" .vcf.gz)"
  echo "[$(date '+%F %T')] Processing ${base}" | tee -a "$log_file"

  snp_pass_vcf="${OUTPUT_DIR}/${base}.snp_pass.vcf.gz"
  renamed_vcf="${OUTPUT_DIR}/${base}.renamed.vcf.gz"
  ready_vcf="${OUTPUT_DIR}/${base}.ready.vcf.gz"
  sample_map="${OUTPUT_DIR}/${base}.sample_map.tsv"

  bcftools view \
    -m2 -M2 \
    -v snps \
    -f PASS \
    -Oz \
    -o "$snp_pass_vcf" \
    "$vcf"

  bcftools annotate \
    --rename-chrs "$NEWNAME_FILE" \
    -Oz \
    -o "$renamed_vcf" \
    "$snp_pass_vcf"

  paste \
    <(bcftools query -l "$renamed_vcf") \
    <(bcftools query -l "$renamed_vcf" | sed 's/_/-/g') \
    > "$sample_map"

  bcftools reheader \
    -s "$sample_map" \
    -o "$ready_vcf" \
    "$renamed_vcf"

  tabix -f -p vcf "$ready_vcf"

  processed=$((processed + 1))
done

echo "[$(date '+%F %T')] Completed 01_prepare_ready_vcfs.sh. Processed ${processed} file(s)." | tee -a "$log_file"
