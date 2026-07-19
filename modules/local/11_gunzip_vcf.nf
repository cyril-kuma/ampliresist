// Expand the bgzipped VCFs back to plain text, per run.
//
// Downstream stage 2 (vcfR) reads plain VCFs. Deriving them from the bgzipped
// output — rather than tapping each caller's raw .vcf — keeps the downstream
// variant-caller agnostic, and replaces the manual `cp -r … && gunzip` step that
// used to be run by hand after the pipeline finished.

process GUNZIP_VCF {
    tag "${run_name}"
    label 'process_low'
    container "quay.io/biocontainers/tabix@sha256:106e72ca3c7ca98c12b3971ba3d2699f4ec63673976f6037a38ebf1d46727515"  // htslib/tabix 1.11

    input:
    tuple val(run_name), path(vcf_gz)

    output:
    tuple val(run_name), path("*.vcf"), emit: vcf

    script:
    """
    for f in *.vcf.gz; do
        gunzip -c "\$f" > "\${f%.gz}"
    done
    """

    stub:
    """
    touch ${run_name}_stub.vcf
    """
}
