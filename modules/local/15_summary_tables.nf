// Stage 15 - merged summary tables (the Genetic Report Card view) for the cohort.
//
// Runs once over ALL flowcells pooled: resistance-allele frequencies are a
// property of the sample set, not of a flowcell. The R stage keys samples on
// (ont_multiplex_group, ont_barcode), since barcode01 exists on every flowcell.

process SUMMARY_TABLES {
    tag "${params.cohort_name}"
    label 'process_low'
    container "${params.container_registry}/ampliresist-r:1.0.0"

    input:
    path coverage_qc
    path genotype_calls
    path frequencies
    val  run_names
    path resources
    path samplesheet

    output:
    path "04_summary", emit: block

    script:
    """
    mkdir -p run
    cp -rL ${coverage_qc} ${genotype_calls} ${frequencies} run/

    export HOME="\$PWD"
    export NANORAVE_RUN_NAME='${params.cohort_name}'
    export NANORAVE_RUNS='${run_names.join("|")}'
    export NANORAVE_ANALYSIS_DATE='${params.analysis_date}'
    export NANORAVE_MIN_COV='${params.min_cov}'
    export NANORAVE_RESOURCE_DIR="\$PWD/${resources}"
    export NANORAVE_METADATA_FILE="\$PWD/${samplesheet}"
    export NANORAVE_MULTIPLEX_SHEET='${params.multiplex_sheet}'
    export NANORAVE_RUN_OUTPUT_DIR="\$PWD/run"

    04_summary_tables.R

    mv run/04_summary .
    """

    stub:
    """
    mkdir -p 04_summary && touch 04_summary/stub.csv
    """
}
