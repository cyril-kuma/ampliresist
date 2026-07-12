// Index the reference FASTA (Clair3 needs the .fai).

process SAMTOOLS_FAIDX {
    tag "${meta.id}"
    label 'process_single'
    container "quay.io/biocontainers/samtools:1.15.1--h1170115_0"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path(fasta), path("*.fai"), emit: fasta
    path "versions.yml",                         emit: versions

    script:
    """
    samtools faidx ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta}.fai versions.yml
    """
}
