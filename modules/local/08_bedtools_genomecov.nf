// Per-base coverage across the amplicon, consumed by downstream stage 1.

process BEDTOOLS_GENOMECOV {
    tag "${meta.id}"
    label 'process_low'
    container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"

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
