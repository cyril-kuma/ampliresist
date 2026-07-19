// Per-base coverage across the amplicon, consumed by downstream stage 1.

process BEDTOOLS_GENOMECOV {
    tag "${meta.id}"
    label 'process_low'
    container "quay.io/biocontainers/bedtools@sha256:9199479a4142cb21a7dbd9f215a4f29811d95927cf57b49f7163d3c2d1d70abb"  // bedtools 2.29.2

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    path "versions.yml",                 emit: versions

    script:
    """
    bedtools genomecov -split -ibam ${bam} -bga | bedtools sort > ${meta.id}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/^bedtools v//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.bedGraph versions.yml
    """
}
