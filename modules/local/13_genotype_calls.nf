// Stage 13 - drug-resistance genotype calls across the pooled cohort.
//
// Runs once over ALL flowcells pooled: resistance-allele frequencies are a
// property of the sample set, not of a flowcell. The R stage keys samples on
// (ont_multiplex_group, ont_barcode), since barcode01 exists on every flowcell.

process GENOTYPE_CALLS {
    tag "${params.cohort_name}"
    label 'process_medium'
    container "drag1-downstream:1.0"

    input:
    path vcfs
    val  run_names
    path resources
    path samplesheet

    output:
    path "02_genotype_calls", emit: block

    script:
    """
    mkdir -p input/variant_calling_unzip run
    cp -L *.vcf input/variant_calling_unzip/

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
    export NANORAVE_MIN_DP='${params.min_cov}'

    02_genotype_calls.R

    mv run/02_genotype_calls .
    """

    stub:
    """
    mkdir -p 02_genotype_calls && touch 02_genotype_calls/stub.csv
    """
}
