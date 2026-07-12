// Refuse to pool the same flowcell twice.
//
// This pipeline exists to merge several flowcell batches into ONE cohort, so the
// most damaging input error is a duplicated run: every sample in the repeated
// flowcell is then counted more than once under different sample_ids, and the
// resulting allele frequencies are pseudoreplicated. Nothing downstream can
// detect this -- the run completes cleanly and the numbers look plausible (they
// scale by a suspiciously clean integer factor).
//
// SORT_FASTQS fingerprints each run by hashing a sample of its read IDs (ONT
// UUIDs, unique per sequencing run). Identical fingerprints mean identical reads.

process VERIFY_RUNS_DISTINCT {
    label 'process_single'

    input:
    path fingerprints

    output:
    path "run_fingerprints.tsv", emit: report

    script:
    """
    printf 'run\\tread_id_fingerprint\\n' > run_fingerprints.tsv
    for f in *.fingerprint; do
        printf '%s\\t%s\\n' "\${f%.fingerprint}" "\$(cat "\$f")" >> run_fingerprints.tsv
    done

    dupes=\$(tail -n +2 run_fingerprints.tsv | cut -f2 | sort | uniq -d)
    if [ -n "\$dupes" ]; then
        echo "ERROR: two or more runs in --input contain the SAME reads." >&2
        echo "" >&2
        for d in \$dupes; do
            echo "  identical runs: \$(awk -v h="\$d" '\$2==h {printf "%s ", \$1}' run_fingerprints.tsv)" >&2
        done
        echo "" >&2
        echo "Pooling a flowcell twice pseudoreplicates every sample it carries and" >&2
        echo "invalidates the cohort allele frequencies. Fix --input and re-run." >&2
        exit 1
    fi

    echo "All \$(tail -n +2 run_fingerprints.tsv | wc -l) run(s) carry distinct reads."
    """

    stub:
    """
    printf 'run\\tread_id_fingerprint\\n' > run_fingerprints.tsv
    """
}
