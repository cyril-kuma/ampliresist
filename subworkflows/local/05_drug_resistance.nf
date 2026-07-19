//
// 05 — Drug-resistance analysis, run ONCE over the pooled cohort.
//
//   bedGraphs ──> 12 amplicon coverage ──┬──> 14 frequencies ──> 15 summary tables
//                                        │
//   VCFs ───────> 13 genotype calls ─────┴──> 17 complexity of infection
//                        │
//   per_call.tsv ────────┘
//
// Every flowcell's output is pooled into a single analysis: resistance-allele
// frequencies are a property of the sample set, not of a flowcell, so splitting
// them per run would compute the wrong denominators. Per-flowcell QC is still
// emitted separately (01_qc/, 02_coverage/).
//
// The R stages key samples on (ont_multiplex_group, ont_barcode) — barcode01
// exists on every flowcell, so a barcode-only join would mis-assign samples once
// the runs are pooled.
//
// 12 -> 13 -> 14 -> 15. Stage 13 consumes stage 12's callability table: a
// specimen below the coverage threshold must be NOT CALLABLE, not wild type. They are separate processes because the stages are file-coupled —
// each reads the previous blocks off disk — so editing, say, the frequency tables
// re-runs only stage 14 under -resume.
//

include { AMPLICON_COVERAGE      } from '../../modules/local/12_amplicon_coverage'
include { GENOTYPE_CALLS         } from '../../modules/local/13_genotype_calls'
include { RESISTANCE_FREQUENCIES } from '../../modules/local/14_resistance_frequencies'
include { SUMMARY_TABLES         } from '../../modules/local/15_summary_tables'
include { COMPLEXITY_OF_INFECTION } from '../../modules/local/17_complexity_of_infection'
include { CLEAN_METADATA          } from '../../modules/local/18a_clean_metadata'
include { PUBLICATION_PLOTS       } from '../../modules/local/18_publication_plots'

workflow DRUG_RESISTANCE {
    take:
    bedgraphs     // [ bedGraph, ... ]  every run, collected
    vcfs          // [ vcf, ... ]       every run, collected
    per_calls     // [ per_call.tsv, ... ] every run, collected
    run_names     // val: list of run ids in the cohort
    resources     // path: assets/resources
    samplesheet   // path: sample metadata workbook (.xlsx) — authoritative
    geo_dir       // path: Ghana boundary GeoJSON directory (assets/geo)

    main:
    AMPLICON_COVERAGE(bedgraphs, run_names, resources, samplesheet)

    // 13 depends on 12: a genotype cannot be called honestly without knowing
    // whether the amplicon was callable in that specimen. Without it, a
    // low-coverage specimen with no variant record is indistinguishable from a
    // wild-type one, and gets reported as wild type.
    GENOTYPE_CALLS(AMPLICON_COVERAGE.out.block, vcfs, run_names, resources, samplesheet)

    RESISTANCE_FREQUENCIES(
        AMPLICON_COVERAGE.out.block,
        GENOTYPE_CALLS.out.block,
        run_names, resources, samplesheet
    )

    SUMMARY_TABLES(
        AMPLICON_COVERAGE.out.block,
        GENOTYPE_CALLS.out.block,
        RESISTANCE_FREQUENCIES.out.block,
        run_names, resources, samplesheet
    )

    // 17 - complexity of infection, from within-sample heterozygosity. Depends on
    // the artefact catalogue emitted by stage 13: without it every specimen looks
    // polygenomic.
    COMPLEXITY_OF_INFECTION(per_calls, GENOTYPE_CALLS.out.block)

    // Cleaned ecological metadata is DERIVED from the authoritative workbook
    // inside the workflow (default) and propagated to stage 18 by channel.
    // --plot_metadata overrides it only for standalone/independent plotting.
    if (params.plot_metadata) {
        ch_plot_metadata = channel.value(file(params.plot_metadata, checkIfExists: true))
    } else {
        CLEAN_METADATA(SUMMARY_TABLES.out.block, samplesheet)
        ch_plot_metadata = CLEAN_METADATA.out.clean_metadata
    }

    // Roadmap-driven figures integrate resistance, callability, complexity and
    // cleaned ecological metadata, and run once per pooled cohort.
    PUBLICATION_PLOTS(
        SUMMARY_TABLES.out.block,
        GENOTYPE_CALLS.out.block,
        COMPLEXITY_OF_INFECTION.out.block,
        ch_plot_metadata,
        geo_dir
    )

    emit:
    coverage_qc    = AMPLICON_COVERAGE.out.block
    genotype_calls = GENOTYPE_CALLS.out.block
    frequencies    = RESISTANCE_FREQUENCIES.out.block
    summary        = SUMMARY_TABLES.out.block
    complexity     = COMPLEXITY_OF_INFECTION.out.block
    plots          = PUBLICATION_PLOTS.out.plots
}
