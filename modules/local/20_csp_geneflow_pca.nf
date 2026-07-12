// Stage 20 — CSP population genetics: diversity, differentiation and structure.
//
// Nucleotide diversity (pi), Tajima's D, pairwise Fst between ecological zones,
// and PCA of population structure.

process CSP_GENEFLOW_PCA {
    tag "${params.cohort_name}"
    label 'process_low'
    container "drag1-popgen:1.0"

    input:
    path merged
    path maf

    output:
    path "geneflow", emit: geneflow

    script:
    """
    mkdir -p geneflow
    csp_04_geneflow_pca.sh \\
        --merged-vcf population.vcf.gz \\
        --maf-dir ${maf} \\
        --output-dir geneflow
    """

    stub:
    """
    mkdir -p geneflow && touch geneflow/stub.eigenvec
    """
}
