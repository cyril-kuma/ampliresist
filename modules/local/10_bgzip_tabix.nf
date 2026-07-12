// Compress and index each VCF.

process BGZIP_TABIX {
    tag "${meta.id}"
    label 'process_single'
    container "quay.io/biocontainers/tabix:1.11--hdfd78af_0"

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
