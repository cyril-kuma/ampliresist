// Build the minimap2 index once per reference, reused across all barcodes.

process MINIMAP2_INDEX {
    tag "${meta.id}"
    label 'process_low'
    container "quay.io/biocontainers/minimap2@sha256:1f23d5cfbefb25ef4f9a0ee5b4f78d3b6cb0b3c955028d80e1d8b00bc97e299a"  // minimap2 2.24

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
