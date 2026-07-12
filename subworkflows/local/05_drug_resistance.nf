//
// 05 — Drug-resistance analysis, run ONCE over the pooled cohort.
//
//   bedGraphs ──> 12 amplicon coverage ──┬──> 14 frequencies ──> 15 summary tables
//                                        │
//   VCFs ───────> 13 genotype calls ─────┴──> 17 complexity of infection
//                        │                            │
//   per_call.tsv ────────┘                            └──> 18 integrate metadata
//                                                            (+ callability-bias test)
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
// Stages 12 and 13 are independent and run in parallel; 14 needs both; 15 needs
// all three. They are separate processes because the stages are file-coupled —
// each reads the previous blocks off disk — so editing, say, the frequency tables
// re-runs only stage 14 under -resume.
//

include { AMPLICON_COVERAGE      } from '../../modules/local/12_amplicon_coverage'
include { GENOTYPE_CALLS         } from '../../modules/local/13_genotype_calls'
include { RESISTANCE_FREQUENCIES } from '../../modules/local/14_resistance_frequencies'
include { SUMMARY_TABLES         } from '../../modules/local/15_summary_tables'
include { COMPLEXITY_OF_INFECTION } from '../../modules/local/17_complexity_of_infection'
include { INTEGRATE_METADATA      } from '../../modules/local/18_integrate_metadata'

workflow DRUG_RESISTANCE {
    take:
    bedgraphs     // [ bedGraph, ... ]  every run, collected
    vcfs          // [ vcf, ... ]       every run, collected
    per_calls     // [ per_call.tsv, ... ] every run, collected
    run_names     // val: list of run ids in the cohort
    resources     // path: assets/resources
    samplesheet   // path: sample metadata .xlsx

    main:
    AMPLICON_COVERAGE(bedgraphs, run_names, resources, samplesheet)
    GENOTYPE_CALLS(vcfs, run_names, resources, samplesheet)

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

    // 18 - Objective 1/2 metadata + the callability-bias test. Optional: only
    // runs when both external tables are supplied.
    ch_integrated = Channel.empty()
    if (params.qpcr_metadata && params.host_calls) {
        INTEGRATE_METADATA(
            samplesheet,
            file(params.qpcr_metadata, checkIfExists: true),
            file(params.host_calls,    checkIfExists: true),
            AMPLICON_COVERAGE.out.block,
            COMPLEXITY_OF_INFECTION.out.coi
        )
        ch_integrated = INTEGRATE_METADATA.out.block
    }

    emit:
    coverage_qc    = AMPLICON_COVERAGE.out.block
    genotype_calls = GENOTYPE_CALLS.out.block
    frequencies    = RESISTANCE_FREQUENCIES.out.block
    summary        = SUMMARY_TABLES.out.block
    complexity     = COMPLEXITY_OF_INFECTION.out.block
    integrated     = ch_integrated
}
