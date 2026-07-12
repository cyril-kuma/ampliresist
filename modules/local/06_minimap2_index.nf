// Build the minimap2 index once per reference, reused across all barcodes.

process MINIMAP2_INDEX {
    tag "${meta.id}"
    label 'process_low'
    container "quay.io/biocontainers/minimap2:2.24--h7132678_1"

    input:
    tuple val(meta), path(fasta), path(fai)

    output:
    tuple val(meta), path(fasta), path(fai), path("*.mmi"), emit: index
    path "versions.yml",                                    emit: versions

    script:
    """
    minimap2 -x map-ont -t ${task.cpus} -d ${fasta.baseName}.mmi ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta.baseName}.mmi versions.yml
    """
}
