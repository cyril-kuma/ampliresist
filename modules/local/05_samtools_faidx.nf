// Index the reference FASTA (Clair3 needs the .fai).

process SAMTOOLS_FAIDX {
    tag "${meta.id}"
    label 'process_single'
    container "quay.io/biocontainers/samtools@sha256:c3e0ba2add590177a2e6ea33ae9074dc1f82b99e6913338d3d1c3a70dc78b518"  // samtools 1.15.1

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
