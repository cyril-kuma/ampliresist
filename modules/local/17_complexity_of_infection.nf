// Stage 17 - complexity of infection (polygenomic vs mono-compatible).
//
// P. falciparum is haploid in the host, so a heterozygous call means more than
// one clone. Complexity is therefore read from within-sample heterozygosity
// across the amplicon panel, EXCLUDING the systematic-artefact positions --
// those are het in ~100% of specimens and would classify the whole cohort as
// polyclonal (see docs/CALLABILITY.md).
//
// msp1, the canonical COI marker, is unavailable: it carries zero reads in every
// specimen. The amplicon was never generated.

process COMPLEXITY_OF_INFECTION {
    tag "${params.cohort_name}"
    label 'process_low'
    container "drag1-downstream:1.0"

    input:
    path per_calls        // per_call.tsv, one per run, collected
    path genotype_calls   // 02_genotype_calls/ - carries the artefact catalogue

    output:
    path "05_complexity", emit: block
    path "05_complexity/complexity_of_infection.tsv", emit: coi

    script:
    """
    mkdir -p 05_complexity
    per_call_list=\$(ls *_per_call.tsv | paste -sd,)

    # -L: Nextflow stages a directory input as a SYMLINK, and find will not
    # descend into one without it. Without -L this silently finds nothing.
    catalogue=\$(find -L ${genotype_calls} -name '*_artefact_catalogue.csv' | head -1)
    if [ -z "\$catalogue" ]; then
        echo "ERROR: no artefact catalogue found in ${genotype_calls}." >&2
        echo "Without it the artefact positions cannot be excluded, and every" >&2
        echo "specimen would be miscalled polygenomic. Refusing to guess." >&2
        exit 1
    fi

    06_complexity_of_infection.R \\
        "\$per_call_list" \\
        "\$catalogue" \\
        "${params.min_cov}" \\
        05_complexity/complexity_of_infection.tsv
    """

    stub:
    """
    mkdir -p 05_complexity
    printf 'run_name\\tbarcode\\tn_het_real\\tinfection_class\\n' > 05_complexity/complexity_of_infection.tsv
    """
}
