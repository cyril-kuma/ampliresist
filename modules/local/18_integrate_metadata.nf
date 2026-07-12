// Stage 18 - join Objective 2 (qPCR parasite density) and Objective 1 (bloodmeal
// host) to the genotypes, and run the callability-bias test.
//
// The test is the point: dhfr is callable in only ~28% of specimens, and if
// callability tracked parasite density those frequencies would be computed on a
// biased high-parasitaemia subsample. Cq (Objective 2) is the instrument that
// settles it.
//
// Optional: skipped unless --qpcr_metadata and --host_calls are given.

process INTEGRATE_METADATA {
    tag "${params.cohort_name}"
    label 'process_low'
    container "drag1-downstream:1.0"

    input:
    path samplesheet
    path qpcr_metadata    // Obj2 mosquito_level.csv
    path host_calls       // Obj1 host_call_table.tsv
    path coverage_qc      // 01_coverage_qc/ - per-specimen callability
    path coi              // complexity_of_infection.tsv

    output:
    path "06_integrated", emit: block

    script:
    """
    mkdir -p 06_integrated

    # the per-specimen coverage table, not the per-gene run summary
    cov=\$(find ${coverage_qc} -name '*coverage_by_run_sample*.csv' | head -1)
    if [ -z "\$cov" ]; then
        echo "ERROR: no per-specimen coverage table in ${coverage_qc}" >&2
        exit 1
    fi

    07_integrate_metadata.R \\
        "${samplesheet}" \\
        "${qpcr_metadata}" \\
        "${host_calls}" \\
        "\$cov" \\
        "${coi}" \\
        06_integrated/specimen_master.tsv \\
        | tee 06_integrated/callability_bias_test.txt
    """

    stub:
    """
    mkdir -p 06_integrated
    touch 06_integrated/specimen_master.tsv 06_integrated/callability_bias_test.txt
    """
}
