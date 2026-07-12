#!/usr/bin/env bash
# Build the --input run samplesheet by scanning an ONT data directory.
#
# Each run directory is expected to contain exactly one flowcell subdirectory,
# which in turn holds fastq_pass/ and a sequencing_summary_*.txt.
#
#   bin/make_run_samplesheet.sh <data_dir> [run_name ...] > runs.csv
#
#   # all runs found under 02_data
#   bin/make_run_samplesheet.sh ../02_data > runs.csv
#
#   # just two of them
#   bin/make_run_samplesheet.sh ../02_data LUC_DRAG1_IMPAVESRUN23_20250920 \
#                                          LUC_DRAG1_IMPAVESRUN24_20250921 > runs.csv
set -euo pipefail

DATA_DIR="${1:?Usage: make_run_samplesheet.sh <data_dir> [run_name ...]}"
shift || true
[[ -d "$DATA_DIR" ]] || { echo "Not a directory: $DATA_DIR" >&2; exit 1; }
DATA_DIR="$(cd "$DATA_DIR" && pwd)"

if [[ $# -gt 0 ]]; then
    run_dirs=()
    for r in "$@"; do
        [[ -d "$DATA_DIR/$r" ]] || { echo "No such run: $DATA_DIR/$r" >&2; exit 1; }
        run_dirs+=("$DATA_DIR/$r")
    done
else
    mapfile -t run_dirs < <(find "$DATA_DIR" -mindepth 1 -maxdepth 1 -type d | sort)
fi

echo "run_name,sequencing_dir,sequencing_summary"

for run_dir in "${run_dirs[@]}"; do
    run_name="$(basename "$run_dir")"

    flowcell_dir="$(find "$run_dir" -mindepth 1 -maxdepth 1 -type d | head -n1)"
    if [[ -z "$flowcell_dir" ]]; then
        echo "WARN: no flowcell subdirectory in $run_dir - skipping" >&2
        continue
    fi

    summary="$(find "$flowcell_dir" -maxdepth 1 -type f -name 'sequencing_summary_*.txt' | head -n1)"
    if [[ -z "$summary" ]]; then
        echo "WARN: no sequencing_summary_*.txt in $flowcell_dir - skipping" >&2
        continue
    fi

    if [[ ! -d "$flowcell_dir/fastq_pass" ]]; then
        echo "WARN: no fastq_pass/ in $flowcell_dir - skipping" >&2
        continue
    fi

    printf '%s,%s,%s\n' "$run_name" "$flowcell_dir" "$summary"
done
