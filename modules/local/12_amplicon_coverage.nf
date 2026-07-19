// Stage 12 - per-amplicon coverage QC from the bedGraphs of every flowcell.
//
// Runs once over ALL flowcells pooled: resistance-allele frequencies are a
// property of the sample set, not of a flowcell. The R stage keys samples on
// (ont_multiplex_group, ont_barcode), since barcode01 exists on every flowcell.

process AMPLICON_COVERAGE {
    tag "${params.cohort_name}"
    label 'process_low'
    container "${params.container_registry}/ampliresist-r:1.0.0"

    input:
    path bedgraphs
    val  run_names
    path resources
    path samplesheet

    output:
    path "01_coverage_qc", emit: block

    script:
    """
    mkdir -p input/genome_coverage run
    cp -L *.bedGraph input/genome_coverage/

    export HOME="\$PWD"
    export NANORAVE_RUN_NAME='${params.cohort_name}'
    export NANORAVE_RUNS='${run_names.join("|")}'
    export NANORAVE_ANALYSIS_DATE='${params.analysis_date}'
    export NANORAVE_MIN_COV='${params.min_cov}'
    export NANORAVE_RESOURCE_DIR="\$PWD/${resources}"
    export NANORAVE_METADATA_FILE="\$PWD/${samplesheet}"
    export NANORAVE_MULTIPLEX_SHEET='${params.multiplex_sheet}'
    export NANORAVE_RUN_OUTPUT_DIR="\$PWD/run"
    export NANORAVE_INPUT_DIR="\$PWD/input"

    01_amplicon_coverage.R

    mv run/01_coverage_qc .
    """

    stub:
    """
    mkdir -p 01_coverage_qc && touch 01_coverage_qc/stub.csv
    """
}
