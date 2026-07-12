#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: 02_merge_per_sample_vcfs.sh --input-dir DIR --output-file FILE --metadata-file FILE [--metadata-sheet SHEET]

Build one per-sample VCF by concatenating gene-level VCFs for each barcode,
rename sample to metadata sample_id, then merge all samples into one VCF.
EOF
}

INPUT_DIR=""
OUTPUT_FILE=""
METADATA_FILE=""
METADATA_SHEET="Sheet1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir) INPUT_DIR="$2"; shift 2 ;;
    --output-file) OUTPUT_FILE="$2"; shift 2 ;;
    --metadata-file) METADATA_FILE="$2"; shift 2 ;;
    --metadata-sheet) METADATA_SHEET="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$INPUT_DIR" ]] || { echo "--input-dir is required" >&2; exit 1; }
[[ -n "$OUTPUT_FILE" ]] || { echo "--output-file is required" >&2; exit 1; }
[[ -n "$METADATA_FILE" ]] || { echo "--metadata-file is required" >&2; exit 1; }

command -v bcftools >/dev/null 2>&1 || { echo "bcftools not found in PATH" >&2; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "tabix not found in PATH" >&2; exit 1; }
command -v Rscript >/dev/null 2>&1 || { echo "Rscript not found in PATH" >&2; exit 1; }

[[ -d "$INPUT_DIR" ]] || { echo "Input directory not found: $INPUT_DIR" >&2; exit 1; }
[[ -f "$METADATA_FILE" ]] || { echo "Metadata file not found: $METADATA_FILE" >&2; exit 1; }
mkdir -p "$(dirname "$OUTPUT_FILE")"

shopt -s nullglob
ready_vcfs=("${INPUT_DIR}"/*.ready.vcf.gz)
shopt -u nullglob

if [[ ${#ready_vcfs[@]} -eq 0 ]]; then
  echo "No *.ready.vcf.gz files found in: $INPUT_DIR" >&2
  exit 1
fi

tmp_root="$(dirname "$OUTPUT_FILE")/tmp_merge"
barcode_lists_dir="${tmp_root}/barcode_lists"
per_sample_dir="${tmp_root}/per_sample"
rm -rf "$tmp_root"
mkdir -p "$barcode_lists_dir" "$per_sample_dir"

barcode_sample_map="${tmp_root}/barcode_to_sample.tsv"
Rscript - "$METADATA_FILE" "$METADATA_SHEET" "$barcode_sample_map" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
metadata_file <- args[[1]]
metadata_sheet <- args[[2]]
out_file <- args[[3]]

suppressPackageStartupMessages({
  library(readxl)
})

read_meta <- function(path, sheet) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    as.data.frame(readxl::read_xlsx(path, sheet = sheet))
  }
}

meta <- read_meta(metadata_file, metadata_sheet)
cn <- tolower(gsub("[^a-z0-9]", "", names(meta)))
sample_idx <- which(cn == "sampleid")
barcode_idx <- which(cn == "ontbarcode")
status_idx <- which(cn == "status")

if (length(sample_idx) == 0 || length(barcode_idx) == 0) {
  stop("Metadata must include columns equivalent to sample_id and ont_barcode.")
}

control_ids <- c(
  "KH2", "NC", "PC",
  "Control_HB3", "Control_Dd2", "Control_KH2",
  "control_HB3", "control_Dd2", "control_KH2"
)

out <- data.frame(
  barcode = as.character(meta[[barcode_idx[1]]]),
  sample_id = as.character(meta[[sample_idx[1]]]),
  status = if (length(status_idx) > 0) as.character(meta[[status_idx[1]]]) else NA_character_,
  stringsAsFactors = FALSE
)

out <- out[nzchar(out$barcode) & nzchar(out$sample_id), , drop = FALSE]
out <- out[
  !(out$sample_id %in% control_ids) &
    !(toupper(trimws(out$status)) %in% c("PC", "NC", "CONTROL", "POSITIVE_CONTROL", "NEGATIVE_CONTROL")),
  c("barcode", "sample_id"),
  drop = FALSE
]
out <- unique(out)
write.table(out, out_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
RS

for vcf in "${ready_vcfs[@]}"; do
  base="$(basename "$vcf" .ready.vcf.gz)"
  if [[ "$base" =~ (barcode[0-9]+) ]]; then
    barcode="${BASH_REMATCH[1]}"
    echo "$vcf" >> "${barcode_lists_dir}/${barcode}.list"
  else
    echo "Could not extract barcode from filename: $base" >&2
    exit 1
  fi
done

shopt -s nullglob
barcode_list_files=("${barcode_lists_dir}"/*.list)
shopt -u nullglob
if [[ ${#barcode_list_files[@]} -eq 0 ]]; then
  echo "No barcode groups created from input VCFs." >&2
  exit 1
fi

per_sample_vcfs=()
for list_file in "${barcode_list_files[@]}"; do
  barcode="$(basename "$list_file" .list)"
  sample_id="$(awk -F '\t' -v b="$barcode" '$1 == b { print $2; exit }' "$barcode_sample_map")"
  if [[ -z "$sample_id" ]]; then
    echo "Skipping $barcode because it has no non-control sample_id mapping in metadata." >&2
    continue
  fi

  concat_vcf="${per_sample_dir}/${sample_id}.concat.vcf.gz"
  reheader_vcf="${per_sample_dir}/${sample_id}.vcf.gz"
  sample_map_file="${per_sample_dir}/${sample_id}.sample_map.tsv"

  mapfile_paths=()
  while IFS= read -r path; do
    [[ -n "$path" ]] && mapfile_paths+=("$path")
  done < "$list_file"
  [[ ${#mapfile_paths[@]} -gt 0 ]] || { echo "Empty list for $barcode" >&2; exit 1; }

  bcftools concat -a "${mapfile_paths[@]}" -Oz -o "$concat_vcf"
  original_sample="$(bcftools query -l "$concat_vcf" | head -n 1)"
  [[ -n "$original_sample" ]] || { echo "No sample found in $concat_vcf" >&2; exit 1; }
  printf '%s\t%s\n' "$original_sample" "$sample_id" > "$sample_map_file"
  bcftools reheader -s "$sample_map_file" -o "$reheader_vcf" "$concat_vcf"
  tabix -f -p vcf "$reheader_vcf"
  per_sample_vcfs+=("$reheader_vcf")
done

if [[ ${#per_sample_vcfs[@]} -eq 0 ]]; then
  echo "No per-sample VCFs were created." >&2
  exit 1
fi

bcftools merge "${per_sample_vcfs[@]}" -Oz -o "$OUTPUT_FILE"
tabix -f -p vcf "$OUTPUT_FILE"

echo "Merged ${#per_sample_vcfs[@]} per-sample VCF files into: $OUTPUT_FILE"
