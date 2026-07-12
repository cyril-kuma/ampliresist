//
// 03 — Align every barcode against every reference amplicon, and derive coverage.
//
// This is the fan-out point of the pipeline: 95 barcodes x 7 genes = 665
// alignments per sequencing run.
//
// The reference travels in the SAME tuple as the reads, so a BAM is inseparable
// from the reference it was aligned to. (The original pipeline emitted the
// reference on a separate channel and re-paired it by channel position, which
// left Clair3's --ref_fn bound to the BAM only by ordering.)
//

include { MINIMAP2_ALIGN     } from '../../modules/local/07_minimap2_align'
include { BEDTOOLS_GENOMECOV } from '../../modules/local/08_bedtools_genomecov'

workflow ALIGN_AND_COVERAGE {
    take:
    reads        // [ meta, fastq ]
    references   // [ meta, fasta, fai, mmi ]

    main:
    ch_align_in = reads
        .combine(references)
        .map { meta, fastq, ref_meta, fasta, fai, mmi ->
            [ meta + [ id: "${meta.id}_${ref_meta.id}", gene: ref_meta.id ], fastq, fasta, fai, mmi ]
        }

    MINIMAP2_ALIGN(ch_align_in)
    BEDTOOLS_GENOMECOV(MINIMAP2_ALIGN.out.bam)

    // Re-attach the reference to each BAM by join() on the meta map — never by
    // channel position.
    ch_bam_with_ref = MINIMAP2_ALIGN.out.bam
        .join(ch_align_in.map { meta, fastq, fasta, fai, mmi -> [ meta, fasta, fai ] })

    // Flat channel of every bedGraph across every run: the resistance analysis
    // is a single pooled cohort, not one analysis per flowcell.
    ch_bedgraphs = BEDTOOLS_GENOMECOV.out.bedgraph.map { meta, bedgraph -> bedgraph }

    ch_versions = Channel.empty()
        .mix(MINIMAP2_ALIGN.out.versions.first())
        .mix(BEDTOOLS_GENOMECOV.out.versions.first())

    emit:
    bam             = ch_bam_with_ref       // [ meta, bam, bai, fasta, fai ]
    reads_with_ref  = ch_align_in           // [ meta, fastq, fasta, fai, mmi ]
    bedgraphs       = ch_bedgraphs          // bedGraph (flat, all runs)
    versions        = ch_versions
}
