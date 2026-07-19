#!/usr/bin/env bash
# Unit tests for bin/make_run_samplesheet.sh — canonical discovery, validation
# errors, and the deprecated double-nested fallback. Fully self-contained: builds
# synthetic fixtures, needs no study data, no containers.
#
#   tests/test_run_discovery.sh
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
SCRIPT="$HERE/bin/make_run_samplesheet.sh"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
E="$TMP/stderr"
pass=0; fail=0
ok(){ echo "  PASS: $1"; pass=$((pass+1)); }
no(){ echo "  FAIL: $1"; fail=$((fail+1)); }

# Build a canonical flowcell: <run_dir>/<flow>/fastq_pass/barcode01 + summary.
mkflow(){ local d="$1/$2"; mkdir -p "$d/fastq_pass/barcode01"
          : > "$d/fastq_pass/barcode01/x.fastq.gz"; : > "$d/sequencing_summary_${2}_abc.txt"; }

# 1. canonical layout: one row, run discovered
D="$TMP/canon"; mkflow "$D/RUN_A" flowA
if out=$("$SCRIPT" "$D" 2>"$E"); then
  [ "$(echo "$out" | tail -n +2 | grep -c .)" -eq 1 ] && grep -q RUN_A <<<"$out" \
    && ok "canonical run discovered" || no "canonical row count/content"
else no "canonical run wrongly rejected"; fi

# 2. no flowcell directory -> hard error
D="$TMP/noflow"; mkdir -p "$D/RUN_B"
"$SCRIPT" "$D" >/dev/null 2>"$E" && no "missing flowcell not rejected" \
  || { grep -qi 'no flowcell' "$E" && ok "missing flowcell errors" || no "wrong msg (missing flowcell)"; }

# 3. missing sequencing summary -> hard error
D="$TMP/nosum"; mkdir -p "$D/RUN_C/flowC/fastq_pass/barcode01"
"$SCRIPT" "$D" >/dev/null 2>"$E" && no "missing summary not rejected" \
  || { grep -qi 'sequencing_summary' "$E" && ok "missing summary errors" || no "wrong msg (missing summary)"; }

# 4. two flowcells with fastq_pass -> ambiguous error
D="$TMP/ambig"; mkflow "$D/RUN_D" flow1; mkflow "$D/RUN_D" flow2
"$SCRIPT" "$D" >/dev/null 2>"$E" && no "ambiguous flowcells not rejected" \
  || { grep -qi 'AMBIGUOUS' "$E" && ok "ambiguous flowcells error" || no "wrong msg (ambiguous)"; }

# 5. deprecated double-nesting (<run>/<run>/<flow>) -> warns but resolves
D="$TMP/nested"; mkflow "$D/RUN_E/RUN_E" flowE
if out=$("$SCRIPT" "$D" 2>"$E"); then
  grep -qi 'DEPRECATED' "$E" && grep -q RUN_E <<<"$out" \
    && ok "double-nesting warns + resolves" || no "double-nesting handling"
else no "double-nesting wrongly rejected"; fi

# 6. duplicate run_name (explicit args) -> hard error
"$SCRIPT" "$TMP/canon" RUN_A RUN_A >/dev/null 2>"$E" && no "duplicate run_name not rejected" \
  || { grep -qi 'duplicate run_name' "$E" && ok "duplicate run_name error" || no "wrong msg (duplicate)"; }

echo "run-discovery: $pass passed, $fail failed"
[ "$fail" -eq 0 ]
