// Stage 17 — CSP population genetics: prepare per-sample VCFs.
//
// Keeps PASS biallelic SNPs, renames contigs to chromosome-style names (plink
// needs that), normalises sample IDs and indexes the result.

process CSP_PREPARE_VCFS {
    tag "${params.cohort_name}"
    label 'process_low'
    container "drag1-popgen:1.0"

    input:
    path vcf_gz
    path rename_map

    output:
    path "ready_vcfs", emit: vcfs

    script:
    """
    mkdir -p in ready_vcfs
    cp -L *.vcf.gz in/ 2>/dev/null || true
    cp -L *.vcf.gz.tbi in/ 2>/dev/null || true

    csp_01_prepare_vcfs.sh \\
        --input-dir in \\
        --output-dir ready_vcfs \\
        --newname-file ${rename_map}
    """

    stub:
    """
    mkdir -p ready_vcfs && touch ready_vcfs/stub.vcf.gz
    """
}
