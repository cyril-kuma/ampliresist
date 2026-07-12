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
//
// --allow_duplicate_runs downgrades the abort to a loud warning. That exists for
// pipeline development, where the same flowcell is deliberately copied under
// several run names to exercise the multi-batch path before real batches exist.
// It is not a switch to reach for on real data.

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
        echo "Two or more runs in --input contain the SAME reads:" >&2
        for d in \$dupes; do
            echo "  identical: \$(awk -v h="\$d" '\$2==h {printf "%s ", \$1}' run_fingerprints.tsv)" >&2
        done
        echo "" >&2
        echo "Pooling a flowcell twice pseudoreplicates every sample it carries, so the" >&2
        echo "cohort allele frequencies and population structure are not interpretable." >&2

        if [ "${params.allow_duplicate_runs}" = "true" ]; then
            echo "" >&2
            echo "WARNING: continuing anyway because --allow_duplicate_runs was set." >&2
            echo "         The results are structurally valid but NOT scientifically usable." >&2
        else
            echo "" >&2
            echo "Fix --input, or pass --allow_duplicate_runs true if you are deliberately" >&2
            echo "re-using a flowcell to develop the pipeline." >&2
            exit 1
        fi
    else
        echo "All \$(tail -n +2 run_fingerprints.tsv | wc -l) run(s) carry distinct reads."
    fi
    """

    stub:
    """
    printf 'run\\tread_id_fingerprint\\n' > run_fingerprints.tsv
    """
}
