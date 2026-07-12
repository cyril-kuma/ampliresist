//
// 06 — CSP population genetics.
//
//   VCFs ─> 17 prepare ─> 18 merge ─> 19 MAF ─┬─> 20 pi / Tajima's D / Fst / PCA
//                                             └─> 21 report
//
// A DIFFERENT question from drug resistance: csp is the RTS,S/R21 vaccine target,
// so this describes parasite population structure and diversity (haplotype
// diversity, the TH2R/TH3R network, MAF, Fst and PCA between ecological zones)
// rather than resistance alleles.
//
// It is inherently a cohort analysis — it needs every sample in one VCF — so it
// consumes the pooled cohort, like 05_drug_resistance.
//

include { CSP_PREPARE_VCFS   } from '../../modules/local/17_csp_prepare_vcfs'
include { CSP_MERGE_VCFS     } from '../../modules/local/18_csp_merge_vcfs'
include { CSP_POPULATION_MAF } from '../../modules/local/19_csp_population_maf'
include { CSP_GENEFLOW_PCA   } from '../../modules/local/20_csp_geneflow_pca'
include { CSP_REPORT         } from '../../modules/local/21_csp_report'

workflow POPULATION_GENETICS {
    take:
    vcf_gz        // bgzipped VCFs + indexes, every run
    rename_map    // path: contig rename map (plink needs chromosome-style names)
    samplesheet   // path: sample metadata .xlsx

    main:
    CSP_PREPARE_VCFS(vcf_gz, rename_map)
    CSP_MERGE_VCFS(CSP_PREPARE_VCFS.out.vcfs, samplesheet)
    CSP_POPULATION_MAF(CSP_MERGE_VCFS.out.merged, samplesheet)
    CSP_GENEFLOW_PCA(CSP_MERGE_VCFS.out.merged, CSP_POPULATION_MAF.out.maf)
    CSP_REPORT(
        CSP_MERGE_VCFS.out.merged,
        CSP_POPULATION_MAF.out.maf,
        CSP_GENEFLOW_PCA.out.geneflow,
        samplesheet
    )

    emit:
    merged   = CSP_MERGE_VCFS.out.merged
    maf      = CSP_POPULATION_MAF.out.maf
    geneflow = CSP_GENEFLOW_PCA.out.geneflow
    report   = CSP_REPORT.out.report
}
