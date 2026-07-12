// FreeBayes variant calling (alternative caller; not the validated path here).

process FREEBAYES {
    tag "${meta.id}"
    label 'process_medium'
    container "docker.io/gfanz/freebayes@sha256:d32bbce0216754bfc7e01ad6af18e74df3950fb900de69253107dc7bcf4e1351"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai)

    output:
    tuple val(meta), path("${meta.id}.vcf"), emit: vcf
    path "versions.yml",                     emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    freebayes -f ${fasta} ${args} ${bam} > raw.vcf

    awk '\$1 ~ /^#/ { print; next } { print | "sort -k1,1 -k2,2n" }' raw.vcf > ${meta.id}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freebayes: \$(freebayes --version | sed 's/^version:  *//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.vcf versions.yml
    """
}
