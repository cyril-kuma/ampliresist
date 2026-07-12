// Stage 21 — CSP population genetics: report.
//
// Haplotype diversity and frequencies, the TH2R/TH3R haplotype network, the
// variant/site catalogue, and the MAF + PCA structure plots, rendered to HTML.

process CSP_REPORT {
    tag "${params.cohort_name}"
    label 'process_low'
    container "drag1-popgen:1.0"

    input:
    path merged
    path maf
    path geneflow
    path samplesheet

    output:
    path "report", emit: report

    script:
    """
    mkdir -p report

    # bin/ is on PATH; the renderer needs the directory itself to find the .Rmd.
    script_dir=\$(dirname \$(command -v csp_06_render_report.R))

    csp_05_report_data.R \\
        --merged-vcf population.vcf.gz \\
        --maf-dir ${maf} \\
        --geneflow-dir ${geneflow} \\
        --metadata-file ${samplesheet} \\
        --metadata-sheet '${params.multiplex_sheet}' \\
        --output-report-dir report \\
        --run-name '${params.cohort_name}'

    csp_06_render_report.R \\
        --output-report-dir report \\
        --run-name '${params.cohort_name}' \\
        --script-dir "\$script_dir"
    """

    stub:
    """
    mkdir -p report && touch report/report.html
    """
}
