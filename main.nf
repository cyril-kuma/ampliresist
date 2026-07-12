#!/usr/bin/env nextflow
/*
 * nano-rave-drag1
 *
 * ONT amplicon variant calling and P. falciparum drug-resistance genotyping,
 * end to end in a single pipeline: raw flowcell output -> resistance-allele
 * frequencies and summary tables.
 *
 * Adapted from sanger-pathogens/nano-rave (C) 2022,2023 Genome Research Ltd.
 */

include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

include { PREPARE_READS      } from './subworkflows/local/01_prepare_reads'
include { PREPARE_REFERENCES } from './subworkflows/local/02_prepare_references'
include { ALIGN_AND_COVERAGE } from './subworkflows/local/03_align_and_coverage'
include { CALL_VARIANTS      } from './subworkflows/local/04_call_variants'
include { DRUG_RESISTANCE    } from './subworkflows/local/05_drug_resistance'

include { DUMP_VERSIONS      } from './modules/local/16_dump_versions'

workflow {

    if (params.validate_params) {
        validateParameters()
    }
    log.info(paramsSummaryLog(workflow))

    // Cross-parameter checks the JSON schema cannot express on its own.
    if (params.variant_caller == 'clair3' && !params.clair3_model) {
        error("--clair3_model is required when --variant_caller is 'clair3'. It must match the basecaller chemistry (this data is R10.4.1, dna_r10.4.1_e8.2_400bps_sup@v4.3.0); the Clair3 image only ships R9.4.1 models.")
    }
    if (params.clair3_model && params.variant_caller != 'clair3') {
        log.warn("--clair3_model was given but --variant_caller is '${params.variant_caller}'; it will be ignored.")
    }

    // Make a development run impossible to mistake for a real one later.
    if (params.allow_duplicate_runs) {
        log.warn("""
        ============================================================
         --allow_duplicate_runs is set: DEVELOPMENT MODE
         Duplicate flowcells will NOT abort the run.
         Any cohort allele frequencies produced are
         pseudoreplicated and NOT scientifically usable.
        ============================================================
        """.stripIndent())
    }

    //
    // Input: one row per sequencing run. Any number of runs may be given; the
    // drug-resistance analysis is grouped per run automatically.
    //
    ch_runs = Channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: ',', strip: true)
        .map { row ->
            if (!row.run_name || !row.sequencing_dir || !row.sequencing_summary) {
                error("--input needs columns 'run_name,sequencing_dir,sequencing_summary'; got: ${row}")
            }
            // Relative paths resolve against the pipeline directory, so the
            // bundled test samplesheet works from any launch directory.
            def dir = row.sequencing_dir.startsWith('/')
                ? row.sequencing_dir
                : "${projectDir}/${row.sequencing_dir}"
            def summ = row.sequencing_summary.startsWith('/')
                ? row.sequencing_summary
                : "${projectDir}/${row.sequencing_summary}"

            [ row.run_name,
              file(dir, checkIfExists: true),
              file(summ, checkIfExists: true) ]
        }

    ch_versions = Channel.empty()

    // 01 - reads + QC
    PREPARE_READS(ch_runs)
    ch_versions = ch_versions.mix(PREPARE_READS.out.versions)

    // 02 - references
    PREPARE_REFERENCES(file(params.reference_manifest, checkIfExists: true))
    ch_versions = ch_versions.mix(PREPARE_REFERENCES.out.versions)

    // 03 - alignment + coverage
    ALIGN_AND_COVERAGE(PREPARE_READS.out.reads, PREPARE_REFERENCES.out.references)
    ch_versions = ch_versions.mix(ALIGN_AND_COVERAGE.out.versions)

    // 04 - variant calling
    ch_clair3_model = params.clair3_model
        ? Channel.value(file(params.clair3_model, checkIfExists: true))
        : Channel.value([])

    CALL_VARIANTS(
        ALIGN_AND_COVERAGE.out.bam,
        ALIGN_AND_COVERAGE.out.reads_with_ref,
        ch_clair3_model
    )
    ch_versions = ch_versions.mix(CALL_VARIANTS.out.versions)

    // 05 - drug-resistance analysis, over ALL flowcells pooled into one cohort.
    //
    // Resistance-allele frequencies are a property of the sample set, not of a
    // flowcell, so the denominators are only correct when the runs are analysed
    // together. Per-flowcell QC is still published separately (01_qc, 02_coverage).
    if (params.drug_resistance) {
        ch_run_names = ch_runs.map { run_name, seq_dir, summary -> run_name }.collect()

        DRUG_RESISTANCE(
            ALIGN_AND_COVERAGE.out.bedgraphs.collect(),
            CALL_VARIANTS.out.vcfs_plain.collect(),
            ch_run_names,
            file(params.resources, checkIfExists: true),
            file(params.samplesheet, checkIfExists: true)
        )
    }


    DUMP_VERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
}
