// Stage 18a - generate the cleaned ecological metadata for the publication
// figures FROM the authoritative workbook (samplesheet.xlsx), inside the
// workflow. The cleaned CSV is a derived product propagated to stage 18 by
// channel; it is not a separately maintained input parameter.
//
// Depends on stage 15's summary block for the genotype<->metadata map, so the
// sample identifiers match the analysed cohort exactly (controls disambiguated).

process CLEAN_METADATA {
    tag "${params.cohort_name}"
    label 'process_single'
    container "${params.container_registry}/ampliresist-r:1.0.0"

    input:
    path summary_block          // stage 15 '04_summary' block
    path samplesheet            // authoritative metadata workbook (.xlsx)

    output:
    path "samplesheet_clean_metadata.csv", emit: clean_metadata
    path "versions.yml",                   emit: versions

    script:
    """
    08_clean_metadata.R ${summary_block} ${samplesheet} ${params.multiplex_sheet} samplesheet_clean_metadata.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | sed -n '1s/.*version \\([0-9.]*\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    """
    printf 'sample_id,bioclimatic_zone,sibling_species,collection_site,latitude,longitude,host_feeding_type,pf_ct,sample_status\\n' > samplesheet_clean_metadata.csv
    touch versions.yml
    """
}
