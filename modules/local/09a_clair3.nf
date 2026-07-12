// Clair3 variant calling.
//
// The pinned hkubal/clair3 image only ships R9.4.1 models, but this data is
// R10.4.1 (dna_r10.4.1_e8.2_400bps_sup@v4.3.0). The previous version solved that
// by bind-mounting the model into the container via `docker.runOptions` in a
// side config — which hard-coded a host path and broke portability (and only
// worked under Docker). Here the model is a normal staged input channel, so it
// works identically under Docker, Singularity and any executor.

process CLAIR3 {
    tag "${meta.id}"
    label 'process_medium'
    container "docker.io/hkubal/clair3@sha256:57cf5d20f2ee39c1b91493ad1fb5c1b9fa838691efce818c3139caa5e6c6b974"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai)
    path model

    output:
    tuple val(meta), path("${meta.id}.vcf"), emit: vcf
    path "versions.yml",                     emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    run_clair3.sh \\
        --bam_fn=${bam} \\
        --ref_fn=${fasta} \\
        --model_path=\$PWD/${model} \\
        --threads=${task.cpus} \\
        --platform=ont \\
        --output=. \\
        ${args}

    # Amplicons are small: when the full-alignment stage finds nothing to merge,
    # Clair3 leaves merge_output empty. The pileup calls are still valid.
    if [[ ! -s merge_output.vcf.gz ]]; then
        cp pileup.vcf.gz merge_output.vcf.gz
    fi

    # Sort by coordinate so tabix can index it.
    zcat merge_output.vcf.gz \\
        | awk '\$1 ~ /^#/ { print; next } { print | "sort -k1,1 -k2,2n" }' \\
        > ${meta.id}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh --version 2>&1 | sed 's/^Clair3 //' | head -n1)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.vcf versions.yml
    """
}
