#!/usr/bin/env bash
# =============================================================================
# 05_call_evidence_audit.sh
#
# Read-level evidence for a position the pipeline reported as NOT CALLABLE.
#
# The artefact rule (bin/02_genotype_calls.R) masks a position on a *statistical*
# signature: called heterozygous in ~every specimen and never homozygous. That is
# enough to refuse to report a frequency, but it is not enough to write in a
# thesis. This script produces the read-level evidence that says WHY, by
# discriminating the four mechanisms that could produce a spurious call:
#
#   1. insufficient depth        -> DP at the position
#   2. mismapping / paralogy     -> MAPQ and strand balance, PER ALLELE
#   3. basecaller context error  -> the same position in a CLONAL control, which
#                                   is a single genome and cannot be heterozygous
#   4. genuine mixed infection   -> allele fraction should VARY across specimens
#
# A position that survives 1-3 with clean reads but shows an invariant allele
# fraction across unrelated infections is a fixed-stoichiometry co-amplification,
# not biology. That is the finding for mdr1 Y184F (see docs/CALLABILITY.md).
#
# Usage:
#   bin/05_call_evidence_audit.sh <gene> <pos> <bam_dir> <control_barcode> [n_field]
# Example:
#   bin/05_call_evidence_audit.sh mdr1 551 work/ barcode95 8
# =============================================================================
set -euo pipefail

GENE=${1:?gene, e.g. mdr1}
POS=${2:?1-based position in the gene CDS, e.g. 551}
BAMDIR=${3:?directory to search for <run>_<barcode>_<gene>.sorted.bam}
CTRL_BC=${4:-barcode95}
N_FIELD=${5:-8}

# One read pileup at POS, broken down by the base the read carries.
# Reports MAPQ, strand and soft-clipping SEPARATELY for the reference-carrying
# and alternate-carrying reads: if the alt allele were mismapped noise, its reads
# would be MAPQ 0 / strand-skewed while the ref reads stayed clean.
pileup_one () {
    local label=$1 bam=$2
    [ -f "$bam" ] || { printf '  %-34s BAM not found\n' "$label"; return; }
    samtools view -F 0x904 "$bam" | awk -v P="$POS" -v L="$label" '
    {
        rp = $4; qp = 1; base = ""; n = ""
        for (i = 1; i <= length($6); i++) {
            c = substr($6, i, 1)
            if (c ~ /[0-9]/) { n = n c; continue }
            len = n + 0; n = ""
            if (c == "S" || c == "I")                       { qp += len }
            else if (c == "D" || c == "N")                  { if (P >= rp && P < rp + len) { base = "-"; break } rp += len }
            else if (c == "M" || c == "=" || c == "X")      { if (P >= rp && P < rp + len) { base = toupper(substr($10, qp + (P - rp), 1)); break } rp += len; qp += len }
        }
        if (base == "" || base == "-") next
        cnt[base]++; mq[base] += $5; if ($5 == 0) mq0[base]++
        if (and($2, 16)) minus[base]++
        tot++
    }
    END {
        if (tot == 0) { printf "  %-34s no reads spanning position\n", L; exit }
        # order bases by count, print the top two
        printf "  %-34s depth=%d\n", L, tot
        for (b in cnt) if (cnt[b] > tot * 0.02)
            printf "      base %s  n=%-6d (%5.1f%%)  meanMAPQ=%4.1f  MAPQ0=%-4d  minus-strand=%4.1f%%\n", \
                   b, cnt[b], 100*cnt[b]/tot, mq[b]/cnt[b], mq0[b]+0, 100*minus[b]/cnt[b]
    }'
}

echo "==========================================================================="
echo " Call-evidence audit: ${GENE} position ${POS}"
echo "==========================================================================="
echo
echo "[1] CLONAL CONTROL (${CTRL_BC}) -- a single genome. It CANNOT be heterozygous."
echo "    If the control is also ~50/50 here, the cause is technical (basecaller"
echo "    or alignment). If the control is clean, the chemistry works at this site."
echo
for bam in $(find "$BAMDIR" -name "*_${CTRL_BC}_${GENE}.sorted.bam" | head -2); do
    pileup_one "control $(basename "$bam" | cut -c1-30)" "$(realpath "$bam")"
done

echo
echo "[2] FIELD SPECIMENS -- independent infections."
echo "    A real polymorphism gives homozygotes and a VARIABLE allele fraction."
echo "    An invariant fraction across unrelated infections does not."
echo
find "$BAMDIR" -name "*_${GENE}.sorted.bam" ! -name "*_${CTRL_BC}_*" \
    | head -"$N_FIELD" | while read -r bam; do
    pileup_one "$(basename "$bam" | sed 's/.*_\(barcode[0-9]*\)_.*/field \1/')" "$(realpath "$bam")"
done

echo
echo "Interpretation: depth rules out (1); per-allele MAPQ and strand balance rule"
echo "out (2); a clean clonal control rules out (3); an invariant allele fraction"
echo "across field specimens rules out (4) -- leaving co-amplification of a second"
echo "template as the only mechanism consistent with all four observations."
