// Stage 18 — CSP population genetics: build the population VCF.
//
// Concatenates the gene-level VCFs for each barcode into one per-sample VCF,
// renames each sample to its metadata sample_id, then merges every sample into a
// single cohort VCF.

process CSP_MERGE_VCFS {
    tag "${params.cohort_name}"
    label 'process_medium'
    container "drag1-popgen:1.0"

    input:
    path ready_vcfs
    path samplesheet

    output:
    path "population.vcf.gz*", emit: merged

    script:
    """
    csp_02_merge_vcfs.sh \\
        --input-dir ${ready_vcfs} \\
        --output-file population.vcf.gz \\
        --metadata-file ${samplesheet} \\
        --metadata-sheet '${params.multiplex_sheet}'
    """

    stub:
    """
    touch population.vcf.gz population.vcf.gz.tbi
    """
}
