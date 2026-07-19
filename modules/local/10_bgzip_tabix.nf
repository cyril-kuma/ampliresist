// Compress and index each VCF.

process BGZIP_TABIX {
    tag "${meta.id}"
    label 'process_single'
    container "quay.io/biocontainers/tabix@sha256:106e72ca3c7ca98c12b3971ba3d2699f4ec63673976f6037a38ebf1d46727515"  // htslib/tabix 1.11

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                     emit: versions

    script:
    """
    bgzip -c ${vcf} > ${meta.id}.vcf.gz
    tabix ${meta.id}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(tabix --version | head -n1 | sed 's/^tabix (htslib) //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.vcf.gz ${meta.id}.vcf.gz.tbi versions.yml
    """
}
