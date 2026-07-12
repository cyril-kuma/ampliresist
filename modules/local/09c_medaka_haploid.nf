// Medaka haploid variant calling (calls straight from reads, not the BAM).

process MEDAKA_HAPLOID {
    tag "${meta.id}"
    label 'process_medium'
    container "quay.io/biocontainers/medaka:1.4.4--py38h130def0_0"

    input:
    tuple val(meta), path(fastq), path(fasta), path(fai), path(mmi)

    output:
    tuple val(meta), path("${meta.id}.vcf"), emit: vcf
    path "versions.yml",                     emit: versions

    script:
    """
    medaka_haploid_variant -r ${fasta} -i ${fastq}

    awk '\$1 ~ /^#/ { print; next } { print | "sort -k1,1 -k2,2n" }' \\
        medaka/medaka.annotated.vcf > ${meta.id}.vcf

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
