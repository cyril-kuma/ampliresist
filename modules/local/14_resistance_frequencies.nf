// Stage 14 - resistance and haplotype frequency tables + plots for the cohort.
//
// Runs once over ALL flowcells pooled: resistance-allele frequencies are a
// property of the sample set, not of a flowcell. The R stage keys samples on
// (ont_multiplex_group, ont_barcode), since barcode01 exists on every flowcell.

process RESISTANCE_FREQUENCIES {
    tag "${params.cohort_name}"
    label 'process_low'
    container "${params.container_registry}/ampliresist-r:1.0.0"

    input:
    path coverage_qc
    path genotype_calls
    val  run_names
    path resources
    path samplesheet

    // min_cov is a DECLARED INPUT, not a params reference inside the script.
    // Nextflow's cache key does not track params interpolated into a script body,
    // so `-profile sensitivity_10x -resume` silently replayed the 50x tasks and
    // published them as 10x results. Declared inputs do participate in the hash.
    val  min_cov

    output:
    path "03_frequencies", emit: block

    script:
    """
    mkdir -p run
    cp -rL ${coverage_qc} ${genotype_calls} run/

    export HOME="\$PWD"
    export NANORAVE_RUN_NAME='${params.cohort_name}'
    export NANORAVE_RUNS='${run_names.join("|")}'
    export NANORAVE_ANALYSIS_DATE='${params.analysis_date}'
    export NANORAVE_MIN_COV='${min_cov}'
    export NANORAVE_RESOURCE_DIR="\$PWD/${resources}"
    export NANORAVE_METADATA_FILE="\$PWD/${samplesheet}"
    export NANORAVE_MULTIPLEX_SHEET='${params.multiplex_sheet}'
    export NANORAVE_RUN_OUTPUT_DIR="\$PWD/run"

    03_resistance_frequencies.R

    mv run/03_frequencies .
    """

    stub:
    """
    mkdir -p 03_frequencies && touch 03_frequencies/stub.csv
    """
}
