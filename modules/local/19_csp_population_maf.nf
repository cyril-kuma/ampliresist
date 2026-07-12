// Stage 19 — CSP population genetics: per-population MAF.
//
// Splits the cohort VCF by ecological zone (the grouping is a --params setting,
// not hard-coded) and computes minor-allele-frequency summaries per group.

process CSP_POPULATION_MAF {
    tag "${params.cohort_name}"
    label 'process_low'
    container "drag1-popgen:1.0"

    input:
    path merged
    path samplesheet

    output:
    path "maf", emit: maf

    script:
    """
    mkdir -p maf
    csp_03_population_maf.sh \\
        --merged-vcf population.vcf.gz \\
        --metadata-file ${samplesheet} \\
        --metadata-sheet '${params.multiplex_sheet}' \\
        --output-dir maf \\
        --coastal-regions '${params.csp_coastal_regions}' \\
        --middlebelt-regions '${params.csp_middlebelt_regions}' \\
        --savannah-regions '${params.csp_savannah_regions}'
    """

    stub:
    """
    mkdir -p maf && touch maf/stub.frq
    """
}
