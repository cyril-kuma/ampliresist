// Medaka diploid variant calling (alternative to Clair3; not the validated path
// for this project — kept because nano-rave supports it).

process MEDAKA {
    tag "${meta.id}"
    label 'process_medium'
    container "quay.io/biocontainers/medaka:1.4.4--py38h130def0_0"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.vcf"), emit: vcf
    path "versions.yml",                     emit: versions

    script:
    """
    medaka_variant -f ${fasta} -i ${bam}

    awk '\$1 ~ /^#/ { print; next } { print | "sort -k1,1 -k2,2n" }' \\
        medaka_variant/round_1.vcf > ${meta.id}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$(medaka --version | sed 's/^medaka //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.vcf versions.yml
    """
}
