//
// 01 — Reads: concatenate each barcode's chunks, attach a meta map, and QC.
//
// Emits one [ meta, fastq ] per barcode, where meta carries the run name and
// barcode so every downstream artefact can be traced back to its sample.
//

include { SORT_FASTQS          } from '../../modules/local/01_sort_fastqs'
include { VERIFY_RUNS_DISTINCT } from '../../modules/local/01b_verify_runs_distinct'
include { NANOPLOT    } from '../../modules/local/02_nanoplot'
include { PYCOQC      } from '../../modules/local/03_pycoqc'

workflow PREPARE_READS {
    take:
    runs   // [ run_name, sequencing_dir, sequencing_summary ]

    main:
    SORT_FASTQS(runs)

    // Refuse to pool the same flowcell twice: that would pseudoreplicate every
    // sample it carries and silently invalidate the cohort frequencies.
    VERIFY_RUNS_DISTINCT(SORT_FASTQS.out.fingerprint.collect())

    ch_reads = SORT_FASTQS.out.fastqs
        .transpose()
        .map { run_name, fastq ->
            def barcode = (fastq.name - "${run_name}_" - '.fastq.gz')
            [ [ id: "${run_name}_${barcode}", run_name: run_name, barcode: barcode ], fastq ]
        }

    NANOPLOT(ch_reads)
    PYCOQC(runs.map { run_name, seq_dir, summary -> [ [ run_name: run_name ], summary ] })

    ch_versions = Channel.empty()
        .mix(SORT_FASTQS.out.versions.first())
        .mix(NANOPLOT.out.versions.first())
        .mix(PYCOQC.out.versions.first())

    emit:
    reads    = ch_reads      // [ meta, fastq ]
    versions = ch_versions
}
