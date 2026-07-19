// Stage 18 - consolidated, cohort-level publication figures.
//
// The process consumes only finalized analysis blocks plus cleaned ecological metadata;
// it never reaches into publishDir, so it remains cacheable and portable.

process PUBLICATION_PLOTS {
    tag "${params.cohort_name}"
    label 'process_low'
    container "${params.container_registry}/ampliresist-figures:1.0.0"

    input:
    path summary
    path genotype_calls
    path complexity
    path plot_metadata
    path geo_dir

    output:
    path "05_plots", emit: plots

    script:
    """
    mkdir -p resistance_input 05_plots
    export NANORAVE_PLOT_SUITE_VERSION=4
    ln -s "\$PWD/${summary}" resistance_input/04_summary
    ln -s "\$PWD/${genotype_calls}" resistance_input/02_genotype_calls
    ln -s "\$PWD/${complexity}" resistance_input/05_complexity

    python3 \$(command -v 07_build_roadmap_figures.py) \\
        --result-dir resistance_input \\
        --metadata "\$PWD/${plot_metadata}" \\
        --geo-dir "\$PWD/${geo_dir}" \\
        --outdir 05_plots
    """

    stub:
    """
    mkdir -p 05_plots
    touch 05_plots/plots.stub
    """
}
