//
// 02 — References: normalise, index, and build a minimap2 index per amplicon.
//
// Emits one fully-prepared reference per gene: [ meta, fasta, fai, mmi ].
// Everything downstream carries the reference in the same tuple as the reads, so
// a BAM can never drift away from the reference it was aligned to.
//

include { NORMALISE_FASTA } from '../../modules/local/04_normalise_fasta'
include { SAMTOOLS_FAIDX  } from '../../modules/local/05_samtools_faidx'
include { MINIMAP2_INDEX  } from '../../modules/local/06_minimap2_index'

workflow PREPARE_REFERENCES {
    take:
    reference_manifest   // path: CSV with reference_id,reference_path

    main:
    ch_manifest = Channel.fromPath(reference_manifest)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            if (!row.reference_id || !row.reference_path) {
                error("--reference_manifest needs columns 'reference_id,reference_path'; got: ${row}")
            }
            // Relative paths resolve against the pipeline directory, so the
            // manifest works no matter where the run is launched from.
            def p = row.reference_path.trim()
            def fasta = p.startsWith('/') ? file(p, checkIfExists: true)
                                          : file("${projectDir}/${p}", checkIfExists: true)
            [ [ id: row.reference_id.trim() ], fasta ]
        }

    NORMALISE_FASTA(ch_manifest)
    SAMTOOLS_FAIDX(NORMALISE_FASTA.out.fasta)
    MINIMAP2_INDEX(SAMTOOLS_FAIDX.out.fasta)

    ch_versions = Channel.empty()
        .mix(NORMALISE_FASTA.out.versions.first())
        .mix(SAMTOOLS_FAIDX.out.versions.first())
        .mix(MINIMAP2_INDEX.out.versions.first())

    emit:
    references = MINIMAP2_INDEX.out.index   // [ meta, fasta, fai, mmi ]
    versions   = ch_versions
}
