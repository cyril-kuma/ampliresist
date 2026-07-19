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

    //
    // Preflight — fail fast, before any compute, with actionable messages.
    // Covers what the JSON schema cannot: asset COMPLETENESS (not just existence
    // of a top path), model files, geo files, workbook type, and run uniqueness.
    //
    def preflight = []

    // Reference panel: every FASTA named in the manifest must exist.
    file(params.reference_manifest).readLines().drop(1).each { row ->
        if (!row?.trim()) return
        def cols = row.split(',')
        def rp = cols.size() > 1 ? cols[1].trim() : ''
        def rf = rp.startsWith('/') ? file(rp) : file("${projectDir}/${rp}")
        if (!rf.exists()) preflight << "reference FASTA missing for '${cols[0]}': ${rf}"
    }

    // Clair3 model: all four TensorFlow files must be present.
    if (params.variant_caller == 'clair3' && params.clair3_model) {
        ['pileup.index', 'pileup.data-00000-of-00001',
         'full_alignment.index', 'full_alignment.data-00000-of-00001'].each { f ->
            if (!file("${params.clair3_model}/${f}").exists())
                preflight << "Clair3 model file missing: ${params.clair3_model}/${f} (see assets/references/clair3_models/PROVENANCE.md)"
        }
    }

    if (params.drug_resistance) {
        // Geographic assets required by stage 18.
        ['ghana_ADM0.geojson', 'ghana_ADM1.geojson'].each { g ->
            if (!file("${params.plot_geo_dir}/${g}").exists())
                preflight << "geo asset missing: ${params.plot_geo_dir}/${g} (default is assets/geo; see assets/geo/PROVENANCE.md)"
        }
        // The metadata source of truth must be a provided .xlsx workbook.
        if (!params.samplesheet)
            preflight << "--samplesheet is required for the drug-resistance analysis: an .xlsx workbook with a '${params.multiplex_sheet}' sheet."
        else if (!(params.samplesheet ==~ /(?i).*\.xlsx$/))
            preflight << "--samplesheet must be an .xlsx workbook (read via read_xlsx): ${params.samplesheet}"
        else if (!file(params.samplesheet).exists())
            preflight << "--samplesheet not found: ${params.samplesheet}"
    }

    // Run manifest: run_name keys the metadata join and must be unique.
    def runNames = file(params.input).readLines().drop(1)
        .findAll { it?.trim() }.collect { it.split(',')[0].trim() }
    def counts = [:]
    runNames.each { counts[it] = (counts[it] ?: 0) + 1 }
    def dupRuns = counts.findAll { e -> e.value > 1 }.collect { e -> e.key }
    if (dupRuns) preflight << "duplicate run_name(s) in --input: ${dupRuns.join(', ')}"

    if (preflight) {
        error("Preflight validation failed:\n  - " + preflight.join("\n  - "))
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
            CALL_VARIANTS.out.per_call.map { _run, tsv -> tsv }.collect(),
            ch_run_names,
            file(params.resources, checkIfExists: true),
            file(params.samplesheet, checkIfExists: true),
            file(params.plot_geo_dir, checkIfExists: true)
        )
    }


    DUMP_VERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
}
