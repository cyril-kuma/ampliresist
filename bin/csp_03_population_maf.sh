#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: 03_compute_population_maf.sh \
  --merged-vcf FILE \
  --metadata-file FILE \
  --metadata-sheet SHEET \
  --output-dir DIR \
  --coastal-regions "Region1|Region2" \
  --middlebelt-regions "Region1|Region2" \
  --savannah-regions "Region1|Region2"

Creates population-specific VCFs and MAF summaries.
EOF
}

MERGED_VCF=""
METADATA_FILE=""
METADATA_SHEET="Sheet1"
OUTPUT_DIR=""
COASTAL_REGIONS=""
MIDDLEBELT_REGIONS=""
SAVANNAH_REGIONS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --merged-vcf) MERGED_VCF="$2"; shift 2 ;;
    --metadata-file) METADATA_FILE="$2"; shift 2 ;;
    --metadata-sheet) METADATA_SHEET="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --coastal-regions) COASTAL_REGIONS="$2"; shift 2 ;;
    --middlebelt-regions) MIDDLEBELT_REGIONS="$2"; shift 2 ;;
    --savannah-regions) SAVANNAH_REGIONS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$MERGED_VCF" ]] || { echo "--merged-vcf is required" >&2; exit 1; }
[[ -n "$METADATA_FILE" ]] || { echo "--metadata-file is required" >&2; exit 1; }
[[ -n "$OUTPUT_DIR" ]] || { echo "--output-dir is required" >&2; exit 1; }

for cmd in bcftools vcftools bgzip tabix plink Rscript; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "Missing dependency in PATH: $cmd" >&2; exit 1; }
done

[[ -f "$MERGED_VCF" ]] || { echo "Merged VCF not found: $MERGED_VCF" >&2; exit 1; }
[[ -f "$METADATA_FILE" ]] || { echo "Metadata file not found: $METADATA_FILE" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR"

VCF_SAMPLES_FILE="${OUTPUT_DIR}/vcf_samples.txt"
bcftools query -l "$MERGED_VCF" > "$VCF_SAMPLES_FILE"

Rscript - "$METADATA_FILE" "$METADATA_SHEET" "$VCF_SAMPLES_FILE" "$OUTPUT_DIR" "$COASTAL_REGIONS" "$MIDDLEBELT_REGIONS" "$SAVANNAH_REGIONS" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
metadata_file <- args[[1]]
metadata_sheet <- args[[2]]
vcf_samples_file <- args[[3]]
out_dir <- args[[4]]
coastal_regions <- args[[5]]
middlebelt_regions <- args[[6]]
savannah_regions <- args[[7]]

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
})

read_meta <- function(path, sheet) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, stringsAsFactors = FALSE)
  } else {
    as.data.frame(readxl::read_xlsx(path, sheet = sheet))
  }
}

parse_regions <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  trimws(unlist(strsplit(x, "\\|", fixed = FALSE)))
}

meta <- read_meta(metadata_file, metadata_sheet)
need <- c("sample_id", "collection_region")
missing_cols <- setdiff(need, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Metadata missing required columns: ", paste(missing_cols, collapse = ", "))
}

vcf_samples <- readLines(vcf_samples_file, warn = FALSE)
meta2 <- meta %>%
  mutate(sample_id = as.character(sample_id),
         collection_region = as.character(collection_region)) %>%
  filter(!is.na(sample_id), sample_id != "", sample_id %in% vcf_samples)

coastal <- parse_regions(coastal_regions)
middlebelt <- parse_regions(middlebelt_regions)
savannah <- parse_regions(savannah_regions)

meta3 <- meta2 %>%
  mutate(group = case_when(
    collection_region %in% coastal ~ "Coastal",
    collection_region %in% middlebelt ~ "Middlebelt",
    collection_region %in% savannah ~ "Savannah",
    TRUE ~ NA_character_
  ))

write.table(meta3, file.path(out_dir, "filtered_samples.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

write_group <- function(df, group_name, out_name) {
  x <- df %>% filter(group == group_name) %>% distinct(sample_id) %>% pull(sample_id)
  writeLines(x, con = file.path(out_dir, out_name))
}

write_group(meta3, "Savannah", "savannah_samples.txt")
write_group(meta3, "Coastal", "coastal_samples.txt")
write_group(meta3, "Middlebelt", "middlebelt_samples.txt")
RS

# A cohort can legitimately contain no samples from one ecological zone (a zone
# not sampled that season, or a subset run). Skip the empty group with a warning
# rather than aborting the whole analysis — but fail if NONE of them has samples,
# because that means the region names in --*-regions don't match the sample sheet.
GROUPS_WITH_SAMPLES=()
for group in savannah coastal middlebelt; do
  sample_file="${OUTPUT_DIR}/${group}_samples.txt"
  if [[ ! -s "$sample_file" ]]; then
    echo "WARN: no samples for group '$group' - skipping it." >&2
    continue
  fi
  GROUPS_WITH_SAMPLES+=("$group")
done

if [[ ${#GROUPS_WITH_SAMPLES[@]} -eq 0 ]]; then
  echo "ERROR: no samples matched ANY population group." >&2
  echo "       The --coastal/--middlebelt/--savannah-regions values must match the" >&2
  echo "       sample sheet's collection_region column." >&2
  exit 1
fi
echo "Population groups with samples: ${GROUPS_WITH_SAMPLES[*]}"

for group in "${GROUPS_WITH_SAMPLES[@]}"; do
  vcftools --gzvcf "$MERGED_VCF" \
    --keep "${OUTPUT_DIR}/${group}_samples.txt" \
    --recode \
    --out "${OUTPUT_DIR}/${group}_filtered"

  recoded_vcf="${OUTPUT_DIR}/${group}_filtered.recode.vcf"
  [[ -f "$recoded_vcf" ]] || { echo "Expected file missing: $recoded_vcf" >&2; exit 1; }

  bgzip -f "$recoded_vcf"
  tabix -f -p vcf "${recoded_vcf}.gz"

  plink --vcf "${recoded_vcf}.gz" --freq --out "${OUTPUT_DIR}/${group}_maf" --allow-extra-chr

done

combined_maf="${OUTPUT_DIR}/combined_maf.tsv"
printf 'SNP\tMAF\tGroup\n' > "$combined_maf"
for group in "${GROUPS_WITH_SAMPLES[@]}"; do
  # Title-case the group for the report (savannah -> Savannah).
  label="$(tr '[:lower:]' '[:upper:]' <<< "${group:0:1}")${group:1}"
  awk -v g="$label" 'NR > 1 {print $2"\t"$5"\t"g}' "${OUTPUT_DIR}/${group}_maf.frq" >> "$combined_maf"
done

echo "MAF workflow complete. Outputs in: $OUTPUT_DIR"
