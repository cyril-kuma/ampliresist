#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: 04_compute_population_geneflow_pca.sh --merged-vcf FILE --maf-dir DIR --output-dir DIR

Computes:
- nucleotide diversity (pi)
- Tajima's D
- pairwise Fst
- group PCA and combined PCA
EOF
}

MERGED_VCF=""
MAF_DIR=""
OUTPUT_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --merged-vcf) MERGED_VCF="$2"; shift 2 ;;
    --maf-dir) MAF_DIR="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$MERGED_VCF" ]] || { echo "--merged-vcf is required" >&2; exit 1; }
[[ -n "$MAF_DIR" ]] || { echo "--maf-dir is required" >&2; exit 1; }
[[ -n "$OUTPUT_DIR" ]] || { echo "--output-dir is required" >&2; exit 1; }

for cmd in vcftools plink; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "Missing dependency in PATH: $cmd" >&2; exit 1; }
done

[[ -f "$MERGED_VCF" ]] || { echo "Merged VCF not found: $MERGED_VCF" >&2; exit 1; }
[[ -d "$MAF_DIR" ]] || { echo "MAF directory not found: $MAF_DIR" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR"

# Only the groups that stage 19 actually produced a VCF for: a cohort can
# legitimately contain no samples from one ecological zone.
GROUPS_WITH_SAMPLES=()
for group in savannah coastal middlebelt; do
  if [[ -f "${MAF_DIR}/${group}_filtered.recode.vcf.gz" && -s "${MAF_DIR}/${group}_samples.txt" ]]; then
    GROUPS_WITH_SAMPLES+=("$group")
  else
    echo "WARN: group '$group' has no data - skipping it." >&2
  fi
done
[[ ${#GROUPS_WITH_SAMPLES[@]} -gt 0 ]] || { echo "ERROR: no population group has data." >&2; exit 1; }

group_n () { wc -l < "${MAF_DIR}/${1}_samples.txt" | tr -d ' '; }

for group in "${GROUPS_WITH_SAMPLES[@]}"; do
  in_vcf="${MAF_DIR}/${group}_filtered.recode.vcf.gz"

  # Diversity statistics are defined for a single sample; PCA is not.
  vcftools --gzvcf "$in_vcf" --window-pi 100 --out "${OUTPUT_DIR}/${group}_pi"
  vcftools --gzvcf "$in_vcf" --TajimaD 100 --out "${OUTPUT_DIR}/${group}_tajimaD"

  # plink requires >= 2 individuals for a pairwise/PCA computation. A zone with a
  # single sample is a real possibility, and it should not abort the analysis.
  n=$(group_n "$group")
  if [[ "$n" -lt 2 ]]; then
    echo "WARN: group '$group' has ${n} sample(s) - skipping its PCA (needs >= 2)." >&2
    continue
  fi

  plink --vcf "$in_vcf" --make-bed --out "${OUTPUT_DIR}/${group}_plink" --allow-extra-chr
  plink --bfile "${OUTPUT_DIR}/${group}_plink" --pca 5 --out "${OUTPUT_DIR}/${group}_pca" --allow-extra-chr

done

fst_pair () {  # Fst is only defined for a pair of populations that both have data.
  local a="$1" b="$2"
  if [[ ! -s "${MAF_DIR}/${a}_samples.txt" || ! -s "${MAF_DIR}/${b}_samples.txt" ]]; then
    echo "WARN: skipping Fst ${a} vs ${b} - one of the groups has no samples." >&2
    return 0
  fi
  vcftools --gzvcf "$MERGED_VCF" \
    --weir-fst-pop "${MAF_DIR}/${a}_samples.txt" \
    --weir-fst-pop "${MAF_DIR}/${b}_samples.txt" \
    --out "${OUTPUT_DIR}/${a}_vs_${b}_fst"
}

fst_pair savannah coastal
fst_pair savannah middlebelt
fst_pair coastal middlebelt

n_total=$(bcftools query -l "$MERGED_VCF" | wc -l | tr -d ' ')
if [[ "$n_total" -lt 2 ]]; then
  echo "WARN: cohort has ${n_total} sample(s) - skipping the combined PCA (needs >= 2)." >&2
else
  plink --vcf "$MERGED_VCF" --make-bed --out "${OUTPUT_DIR}/combined_pca_input" --allow-extra-chr
  plink --bfile "${OUTPUT_DIR}/combined_pca_input" --pca 5 --out "${OUTPUT_DIR}/combined_pca" --allow-extra-chr
fi

echo "Gene flow and PCA workflow complete. Outputs in: $OUTPUT_DIR"
