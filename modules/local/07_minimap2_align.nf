// Align one barcode against one reference, and sort/index in the same task.
//
// The upstream pipeline ran minimap2 -> samtools view -> samtools sort as three
// separate processes, writing a full SAM and an unsorted BAM to disk for every
// barcode x gene pair (665 of each for a single 96-barcode run). Piping straight
// into `samtools sort` removes both intermediates and cuts three tasks to one.
//
// The reference travels WITH the reads in a single tuple, so a BAM can never be
// paired with the wrong reference — the previous version emitted the reference on
// a separate channel and relied on positional ordering to re-pair them.

process MINIMAP2_ALIGN {
    tag "${meta.id}"
    label 'process_medium'
    container "quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0"

    input:
    tuple val(meta), path(fastq), path(fasta), path(fai), path(mmi)

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: bam
    path "versions.yml",                                             emit: versions

    script:
    def args  = task.ext.args  ?: '-ax map-ont --MD'
    def args2 = task.ext.args2 ?: ''
    """
    minimap2 ${args} -t ${task.cpus} ${mmi} ${fastq} \\
        | samtools sort ${args2} -@ ${task.cpus} -o ${meta.id}.sorted.bam -

    samtools index ${meta.id}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam ${meta.id}.sorted.bam.bai versions.yml
    """
}
