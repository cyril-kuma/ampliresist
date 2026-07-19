#!/usr/bin/env bash
# ============================================================================
# make_run_samplesheet.sh — build the --input run manifest (runs.csv) from an
# ONT data directory, validating the CANONICAL sequencing-run layout.
#
# CANONICAL LAYOUT (exactly one flowcell per run):
#
#   <data_dir>/
#     <run_name>/                       # == ont_multiplex_group in the workbook
#       <flowcell_dir>/                 # e.g. 20250920_1723_X1_FAZ03468_7d2e8f91
#         fastq_pass/                   # barcodeNN/*.fastq.gz
#         sequencing_summary_*.txt      # exactly one
#
# Usage:
#   bin/make_run_samplesheet.sh <data_dir> [run_name ...] > runs.csv
#
#   bin/make_run_samplesheet.sh ../02_data > runs.csv            # all runs
#   bin/make_run_samplesheet.sh ../02_data RUN_A RUN_B > runs.csv # a subset
#
# Discovery is DETERMINISTIC (runs sorted; exactly one flowcell required) and
# validated: a run with zero flowcells, multiple ambiguous flowcells, a missing
# fastq_pass/, or a missing/duplicated sequencing_summary_*.txt is a hard error
# that names the run and shows the expected structure. It does NOT recurse to
# arbitrary depth.
#
# DEPRECATED: a historical double-nested layout (<run>/<run>/<flowcell>/...) is
# still tolerated with a warning — correct it to the canonical layout above.
# ============================================================================
set -euo pipefail

err()  { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARN:  $*" >&2; }

usage() { echo "Usage: make_run_samplesheet.sh <data_dir> [run_name ...] > runs.csv" >&2; exit 2; }

EXPECTED_TREE=$'expected:\n  <data_dir>/<run_name>/<flowcell_dir>/fastq_pass/\n  <data_dir>/<run_name>/<flowcell_dir>/sequencing_summary_*.txt'

[[ $# -ge 1 ]] || usage
DATA_DIR="$1"; shift || true
[[ -d "$DATA_DIR" ]] || err "not a directory: $DATA_DIR"
[[ -r "$DATA_DIR" ]] || err "not readable: $DATA_DIR"
DATA_DIR="$(cd "$DATA_DIR" && pwd)"

# Select run directories: named subset (validated) or every immediate subdir.
declare -a run_dirs=()
if [[ $# -gt 0 ]]; then
    for r in "$@"; do
        [[ -d "$DATA_DIR/$r" ]] || err "no such run directory: $DATA_DIR/$r"
        run_dirs+=("$DATA_DIR/$r")
    done
else
    while IFS= read -r d; do run_dirs+=("$d"); done \
        < <(find "$DATA_DIR" -mindepth 1 -maxdepth 1 -type d | sort)
    [[ ${#run_dirs[@]} -gt 0 ]] || err "no run subdirectories in $DATA_DIR"
fi

# Resolve the single flowcell directory (one containing fastq_pass/) for a run.
resolve_flowcell() {
    local run_dir="$1" run_name="$2"
    local -a subs=() flow=() nested=()
    while IFS= read -r d; do subs+=("$d"); done \
        < <(find "$run_dir" -mindepth 1 -maxdepth 1 -type d | sort)
    [[ ${#subs[@]} -gt 0 ]] || err $''"$run_name"$': no flowcell directory under '"$run_dir"$'\n'"$EXPECTED_TREE"

    # Canonical: immediate subdir(s) that themselves contain fastq_pass/.
    local d
    for d in "${subs[@]}"; do [[ -d "$d/fastq_pass" ]] && flow+=("$d"); done
    if [[ ${#flow[@]} -eq 1 ]]; then printf '%s\n' "${flow[0]}"; return 0; fi
    [[ ${#flow[@]} -le 1 ]] || err "$run_name: AMBIGUOUS — ${#flow[@]} flowcell directories contain fastq_pass/: ${flow[*]}. Keep exactly one flowcell per run."

    # No canonical flowcell: check exactly one deprecated double-nested level.
    for d in "${subs[@]}"; do
        local n
        while IFS= read -r n; do [[ -d "$n/fastq_pass" ]] && nested+=("$n"); done \
            < <(find "$d" -mindepth 1 -maxdepth 1 -type d)
    done
    if [[ ${#nested[@]} -eq 1 ]]; then
        warn "$run_name: DEPRECATED double-nested layout ($run_name/*/flowcell). Correct it to <run>/<flowcell>/. Tolerated for now."
        printf '%s\n' "${nested[0]}"; return 0
    fi
    [[ ${#nested[@]} -eq 0 ]] || err "$run_name: AMBIGUOUS double-nested flowcells: ${nested[*]}"
    err $''"$run_name"$': found no flowcell directory containing fastq_pass/ (checked canonical depth and one legacy level).\n'"$EXPECTED_TREE"
}

echo "run_name,sequencing_dir,sequencing_summary"

declare -A seen_runs=()
declare -a rows=()
for run_dir in "${run_dirs[@]}"; do
    run_name="$(basename "$run_dir")"
    [[ -z "${seen_runs[$run_name]:-}" ]] || err "duplicate run_name '$run_name' — run names must be unique (they key the metadata join)."
    seen_runs[$run_name]=1
    [[ -r "$run_dir" ]] || err "$run_name: run directory not readable: $run_dir"

    flow="$(resolve_flowcell "$run_dir" "$run_name")"

    # Exactly one sequencing summary.
    mapfile -t summ < <(find "$flow" -maxdepth 1 -type f -name 'sequencing_summary_*.txt' | sort)
    [[ ${#summ[@]} -ge 1 ]] || err $''"$run_name"$': no sequencing_summary_*.txt in '"$flow"$'\n'"$EXPECTED_TREE"
    [[ ${#summ[@]} -eq 1 ]] || err "$run_name: AMBIGUOUS — ${#summ[@]} sequencing_summary_*.txt files in $flow: ${summ[*]}"
    [[ -r "${summ[0]}" ]] || err "$run_name: sequencing summary not readable: ${summ[0]}"

    # fastq_pass must hold at least one barcode directory.
    nbc="$(find "$flow/fastq_pass" -mindepth 1 -maxdepth 1 -type d -name 'barcode*' | wc -l)"
    [[ "$nbc" -gt 0 ]] || warn "$run_name: fastq_pass/ contains no barcode* directories."

    rows+=("$run_name,$flow,${summ[0]}")
done

printf '%s\n' "${rows[@]}" | sort
