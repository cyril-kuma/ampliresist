// Per-barcode read QC (length/quality distributions).

process NANOPLOT {
    tag "${meta.id}"
    label 'process_low'
    container "quay.io/biocontainers/nanoplot@sha256:57cffeb19b67d88126a8251069771190a65f73ec694fb2e496e0b4f1a98ef4ab"  // nanoplot 1.38.0

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}_nanoplot/*"), emit: report
    path "versions.yml",                            emit: versions

    script:
    """
    NanoPlot -t ${task.cpus} --fastq ${fastq} -o ${meta.id}_nanoplot

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(NanoPlot --version | sed 's/^NanoPlot //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}_nanoplot && touch ${meta.id}_nanoplot/report.html
    touch versions.yml
    """
}
