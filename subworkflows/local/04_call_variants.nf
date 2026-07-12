//
// 04 — Variant calling. One caller is selected by --variant_caller; all emit the
// same [ meta, vcf ] shape, so everything downstream is caller-agnostic.
//
// Clair3 is the validated path for this assay (P. falciparum is haploid, and the
// R10.4.1 model is supplied via --clair3_model). medaka / medaka_haploid /
// freebayes are carried over from upstream nano-rave and are NOT validated here.
//

include { CLAIR3         } from '../../modules/local/09a_clair3'
include { MEDAKA         } from '../../modules/local/09b_medaka'
include { MEDAKA_HAPLOID } from '../../modules/local/09c_medaka_haploid'
include { FREEBAYES      } from '../../modules/local/09d_freebayes'
include { BGZIP_TABIX    } from '../../modules/local/10_bgzip_tabix'
include { GUNZIP_VCF     } from '../../modules/local/11_gunzip_vcf'

workflow CALL_VARIANTS {
    take:
    bam            // [ meta, bam, bai, fasta, fai ]
    reads_with_ref // [ meta, fastq, fasta, fai, mmi ] - medaka_haploid calls from reads
    clair3_model   // value channel: path to the Clair3 model dir

    main:
    ch_raw_vcf  = Channel.empty()
    ch_versions = Channel.empty()

    if (params.variant_caller == 'clair3') {
        CLAIR3(bam, clair3_model)
        ch_raw_vcf  = CLAIR3.out.vcf
        ch_versions = ch_versions.mix(CLAIR3.out.versions.first())
    }
    else if (params.variant_caller == 'medaka') {
        MEDAKA(bam)
        ch_raw_vcf  = MEDAKA.out.vcf
        ch_versions = ch_versions.mix(MEDAKA.out.versions.first())
    }
    else if (params.variant_caller == 'medaka_haploid') {
        MEDAKA_HAPLOID(reads_with_ref)
        ch_raw_vcf  = MEDAKA_HAPLOID.out.vcf
        ch_versions = ch_versions.mix(MEDAKA_HAPLOID.out.versions.first())
    }
    else if (params.variant_caller == 'freebayes') {
        FREEBAYES(bam)
        ch_raw_vcf  = FREEBAYES.out.vcf
        ch_versions = ch_versions.mix(FREEBAYES.out.versions.first())
    }

    BGZIP_TABIX(ch_raw_vcf)
    ch_versions = ch_versions.mix(BGZIP_TABIX.out.versions.first())

    // Plain VCFs for the analysis stages, grouped per sequencing run. Deriving
    // them from the bgzipped output keeps this caller-agnostic.
    ch_vcf_by_run = BGZIP_TABIX.out.vcf
        .map { meta, vcf, tbi -> [ meta.run_name, vcf ] }
        .groupTuple()

    GUNZIP_VCF(ch_vcf_by_run)

    // GUNZIP_VCF batches per run (a sensible unit of work), but the resistance
    // analysis pools every run, so flatten the plain VCFs into one channel.
    ch_vcfs_plain = GUNZIP_VCF.out.vcf
        .map { run_name, vcfs -> vcfs }
        .flatten()

    emit:
    vcf        = BGZIP_TABIX.out.vcf   // [ meta, vcf.gz, tbi ]
    vcfs_plain = ch_vcfs_plain         // vcf (flat, all runs)
    versions   = ch_versions
}
