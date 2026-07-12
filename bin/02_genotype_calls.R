#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Nextflow stages bin/ onto PATH inside the task, so resolve our own location to
# source the shared preamble (_setup.R -> _pipeline_paths.R) that sits alongside.
# ---------------------------------------------------------------------------
local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  assign("NR_BIN", dirname(normalizePath(f[1])), envir = .GlobalEnv)
})
source(file.path(NR_BIN, "_setup.R"))

########### R script: process Clair3 haploid VCFs into per-sample genotype/haplotype calls
### Current canonical run uses HAPLOID Clair3 calls:
### confident hom-alt = mutant; heterozygous-pattern calls are not counted as
### confident mutants (see the writeups; the former stage 5 feasibility analysis
### that also discussed this was retired in 2026-07).

##################### R script to read in BED files to check coverage

##### Load packages
library("tidyverse")
library("vcfR")
library("readxl")

out_variant_info_fn <- file.path(analysis_dr_2, paste0(MinION_run_name, "_variants_linked_to_metadata_intermediate_file.csv")) # for vcf_varcount_df4
out_variant_info2_fn <- file.path(analysis_dr_2, paste0(MinION_run_name, "_variants_linked_to_metadata_intermediate_file2.csv")) # for vcf_tbl_all_varinfo2
out_keysnp_allsamp_genotypes_fn <- file.path(analysis_dr_2, paste0(MinION_run_name, "_genotype_calls_samplecols.csv")) # for keysnps_allsamps_gt
out_keysnp_allsamp_genotypes_samprows_fn <- file.path(analysis_dr_2, paste0(MinION_run_name, "_genotype_calls_haplotypes_samplerows_v2.csv")) # for snp_calls_t3
out_k13_calls_info_fn <- file.path(analysis_dr_2, paste0(MinION_run_name, "_k13_variant_info.csv")) # for k13_snps_allsamps_df
out_csp_genotypes_long_fn <- file.path(analysis_dr_2, paste0(MinION_run_name, "_csp_nref_genotype_calls_long.csv")) # for vcf_tbl_all_csp
out_csp_genotypes_samps_fn <- file.path(analysis_dr_2, paste0(MinION_run_name, "_csp_genotype_calls_samplecols.csv")) # for csp_samp_calls
out_msp1_genotypes_long_fn <- file.path(analysis_dr_2, paste0(MinION_run_name, "_msp1_nref_genotype_calls_long.csv")) # for vcf_tbl_all_msp1
out_plot_msp1_genetic_distance_fn <- file.path(analysis_dr_2, paste0(MinION_run_name, "_plot_msp1_genetic_distance_NJT.png")) # for msp1 NJT


##### Import data

## list of variants to call
snp_data <- read_excel(file.path(resource_dir, "DR_variant_info_v2.xlsx"), sheet = 1)

## kelch13 info
k13_annotated <- read.csv(file.path(resource_dir, "k13_seq_annotated.csv"))
genetic_code <- read.csv(file.path(resource_dir, "genetic_code.csv"))
k13_who_catalog_path <- file.path(resource_dir, "k13_who_artemisinin_partial_resistance_markers.csv")
k13_who_catalog <- if (file.exists(k13_who_catalog_path)) {
  read.csv(k13_who_catalog_path, stringsAsFactors = FALSE)
} else {
  warning("WHO K13 marker catalogue not found: ", k13_who_catalog_path)
  data.frame(aa_change = character(), classification = character())
}

### Import the vcf files
# construct directory for importing vcf files
nflow_results_dir <- file.path(input_run_dir, "variant_calling_unzip")

if (!dir.exists(nflow_results_dir)) {
  stop("VCF directory not found: ", nflow_results_dir)
}


## Extract the file names
dataFiles <- list() # creates a list
listvcf <- list.files(nflow_results_dir, pattern = "\\.vcf$", full.names = FALSE) # all vcf files
if (length(listvcf) == 0) {
  stop("No VCF files found in: ", nflow_results_dir)
}
for (k in 1:length(listvcf)){
  dataFiles[[k]] <- read.vcfR(file.path(nflow_results_dir, listvcf[k]))
}

# Alternative method for getting all the files
#dataFiles <- lapply(Sys.glob(paste(nflow_results_dir,MinION_run_name,"_barcode*",".vcf",sep="")), read.vcfR)

## Check it out
length(dataFiles)
#dataFiles[[12]]
#dataFiles[[2]]



### Check whether any vcf files don't have any variants

n <- length(dataFiles)

for (i in 1:n) {
  
  kn <- nrow(extract.gt(dataFiles[[i]]))
  
  result <- ifelse(kn>0,
                   print(paste("PASS! Vcf file no. ", i, " contains >=1 variant", sep="")),
                   print(paste("** FAIL! Vcf file no. ", i, " contains zero variants! **", sep="")))
}



######## Create dataframe of variant counts

mylist <- list()

n <- length(dataFiles)

for (i in 1:n) {
  
  # Count number of variants in the vcf file
  kn <- nrow(extract.gt(dataFiles[[i]]))

  mylist[[i]] <- kn
  
}

vcf_varcount_df <- data.frame(do.call("rbind", mylist)) %>%
  rename("variant_count" = "do.call..rbind...mylist.")
vcf_varcount_df$vcf_fn <- listvcf
vcf_varcount_df

vcf_varcount_df$zero_variants <- ifelse(vcf_varcount_df$variant_count<1, TRUE, FALSE)



## amplicon sequence names
vcf_varcount_df2 <- within(vcf_varcount_df, {
  
  # Convert to character
  vcf_fn <- as.character(vcf_fn)
  
  # Convert to integer for the numerical value of unique ID
  amplicon_id <- as.character(as.integer(factor(vcf_fn, levels = unique(vcf_fn))))
  
  # Rename with study format
  amplicon_id <- ifelse(nchar(amplicon_id) == 1, paste0('GH22', multiplex_run, '_SEQ', '00', amplicon_id),
                        ifelse(nchar(amplicon_id) == 2, paste0('GH22', multiplex_run, '_SEQ', '0', amplicon_id),
                               ifelse(nchar(amplicon_id) == 3, paste0('GH22', multiplex_run, '_SEQ', amplicon_id), NA)
                        )
  )
})


### Barcodes
## Extract the barcodes
vcf_varcount_df2$ont_barcode <- vcf_varcount_df2$vcf_fn
vcf_varcount_df2$ont_barcode <- gsub(".*_barcode","barcode", vcf_varcount_df2$vcf_fn, perl = TRUE); head(vcf_varcount_df2)
vcf_varcount_df2$ont_barcode <- gsub("_.*","", vcf_varcount_df2$ont_barcode, perl = TRUE); head(vcf_varcount_df2)



#vcf_varcount_df2$ont_multiplex_run <- multiplex_run ### This will get imported with the metadata
#vcf_varcount_df2$ont_seq_date <- multiplex_seq_date

## Extract amplicon gene names
vcf_varcount_df2$amplicon_gene_target <- vcf_varcount_df2$vcf_fn
vcf_varcount_df2$amplicon_gene_target <- vcf_varcount_df2$amplicon_gene_target %>%
  gsub("^.*_barcode[0-9]+_", "", ., perl = TRUE) %>%
  gsub("\\\\.vcf$", "", ., perl = TRUE)
vcf_varcount_df2$amplicon_gene_target

head(vcf_varcount_df2)


## Add in the sample IDs from metadata file
colnames(metadata)
meta_ed <- metadata %>%
  select(-notes, -patient_id)

# Key on (run, barcode): barcode01 exists on every flowcell, so joining on
# ont_barcode alone would mis-assign samples once runs are pooled.
vcf_varcount_df2$ont_multiplex_group <- nr_run_from_file(vcf_varcount_df2$vcf_fn)

vcf_varcount_df3 <- merge(vcf_varcount_df2, meta_ed,
                          by=c("ont_multiplex_group","ont_barcode"),
                          all.x=T, all.y=F)

# merge() re-sorts by the join key, but `dataFiles` is in list.files() order and
# is indexed by row position below (dataFiles[-empty_vcfs]). Those two orderings
# happened to agree before; restore it explicitly rather than depend on that.
vcf_varcount_df3 <- vcf_varcount_df3[match(listvcf, vcf_varcount_df3$vcf_fn), ]
stopifnot(identical(as.character(vcf_varcount_df3$vcf_fn), as.character(listvcf)))

col_ord1 <- c("sample_id",
              "collection_site",
              "ont_seq_date",
              "ont_multiplex_group",
              "ont_barcode",
              "amplicon_id",
              "amplicon_gene_target",
              "vcf_fn",
              "variant_count",
              "zero_variants")

vcf_varcount_df3 <- vcf_varcount_df3 %>% select(all_of(col_ord1))

head(vcf_varcount_df3,7)
nrow(vcf_varcount_df3)



######## Remove vcf files from list if they have zero variants

empty_vcfs <- which(vcf_varcount_df3$zero_variants == TRUE)

dataFiles_noempty <- dataFiles[-c(empty_vcfs)];

length(dataFiles)
length(dataFiles_noempty)

# keep track of sample order
vcf_varcount_df4 <- vcf_varcount_df3 %>%
  filter(zero_variants==FALSE) %>%
  arrange(vcf_fn)

head(vcf_varcount_df4,7)
nrow(vcf_varcount_df4)



###############
list_length <- length(dataFiles_noempty)
### Testing with 1 sample
dataFiles_noempty[3]
# Extract genotypes + variant info
temp  <- vcfR2tidy(dataFiles_noempty[[list_length]], format_fields = c("GT", "DP", "AF"))
fix   <- temp$fix %>% select(-ChromKey)
gt    <- temp$gt  %>% select(-ChromKey)
# (temp$meta is VCF format metadata, not sample metadata — not needed)
# Combine into 1 table
fixgt <- left_join(fix, gt, by = "POS")



############## loooooooooooop

n2 <- length(dataFiles_noempty)

mylist2 <- list()

for (j in 1:n2) {
  
  ##### Turning vcf into formatted tibble for each gene amplicon
  # Extract genotypes + variant info
  temp <- vcfR2tidy(dataFiles_noempty[[j]], format_fields = c("GT", "DP", "AF"))
  fix <- temp$fix %>% select(-ChromKey)
  gt  <- temp$gt  %>% select(-ChromKey)
  # Combine into 1 table (temp$meta is VCF format metadata, not sample metadata — not needed here)
  fixgt <- left_join(fix, gt, by = "POS")
  # Select columns with anticipated new names
  gtcols <- c("Gene", "Pos", "Ref", "Alt", "DP", "Qual", "Filter", "Pileup_call", "Alignment_call", "GT", "AF")
  # Fiddle with formatting
  ampx_tbl_fix_sampx <- fixgt %>%
    rename(Gene=CHROM,
           Pos=POS,
           Ref=REF,
           Alt=ALT,
           Qual=QUAL,
           Filter=FILTER,
           Pileup_call=P,
           Alignment_call="F",
           GT=gt_GT,
           DP=gt_DP,
           AF=gt_AF) %>%  # rename things
    select(all_of(gtcols))  # select desired vcf components

  ampx_tbl_fix_sampx$seq_count_temp = j
  
  mylist2[[j]] <- ampx_tbl_fix_sampx
  
}

# This may have emptied some tibbles of any content, if there were no variants that make it through filtering
xtemp <- length(dataFiles_noempty)
ytemp <- length(mylist2)

ifelse(xtemp==ytemp,
       "No change in number of vcfs included in the list - all good :-D",
       "Number of vcfs has changed! Danger!")

vcf_tbl_all <- data.frame(do.call("rbind", mylist2)) %>%
  arrange(seq_count_temp) %>%
  mutate(Gene = str_replace(Gene, "ref_target_gene_cds_3D7_DR1_crt", "crt")) %>%
  mutate(Gene = str_replace(Gene, "ref_target_gene_cds_3D7_AG1_csp", "csp")) %>%
  mutate(Gene = str_replace(Gene, "ref_target_gene_cds_3D7_DR1_dhfr", "dhfr")) %>%
  mutate(Gene = str_replace(Gene, "ref_target_gene_cds_3D7_DR1_dhps", "dhps")) %>%
  mutate(Gene = str_replace(Gene, "ref_target_gene_cds_3D7_DR1_k13", "k13")) %>%
  mutate(Gene = str_replace(Gene, "ref_target_gene_cds_3D7_DR1_mdr1", "mdr1")) %>%
  mutate(Gene = str_replace(Gene, "ref_target_gene_cds_3D7_AG1_msp1", "msp1"))

# Add SNP names
vcf_tbl_all2 <- within(vcf_tbl_all, {
  snp_id = paste0(Gene,"_",Pos,"_",Ref,Alt)
})
  


###### Merge with sample metadata

### key step - need to ensure order matches ###

xtemp = n_distinct(vcf_tbl_all2$seq_count_temp)
ytemp = nrow(vcf_varcount_df4)

ifelse(xtemp==ytemp,
       "Unique seq counts matches variant count - all good :-D",
       "Number of variants has changed! Danger!")

vcf_varcount_df4$seq_count_temp <- row.names(vcf_varcount_df4)
vcf_varcount_df4$seq_count_temp <- as.numeric(vcf_varcount_df4$seq_count_temp)

vcf_tbl_all2$seq_count_temp <- as.numeric(vcf_tbl_all2$seq_count_temp)

vcf_tbl_all_varinfo <- merge(vcf_tbl_all2, vcf_varcount_df4,
                             by="seq_count_temp",
                             all.x=T, all.y=T) %>%
  arrange(seq_count_temp)

head(vcf_tbl_all_varinfo,15)

nrow(vcf_tbl_all2)
nrow(vcf_tbl_all_varinfo)

xtemp = nrow(vcf_tbl_all2)
ytemp = nrow(vcf_tbl_all_varinfo)

ifelse(xtemp==ytemp,
       "Number of variants is the same as before merging with metadata - all good :-D",
       "Number of variants has changed! Danger!")


# Check the variant gene matches the vcf call gene names - should be identical
vcf_tbl_all_varinfo$gene_check_temp <- ifelse(vcf_tbl_all_varinfo$Gene==vcf_tbl_all_varinfo$amplicon_gene_target,
       0, 1)
ifelse(sum(vcf_tbl_all_varinfo$gene_check_temp)==0,
       "Gene names from vcf calls and variant metadata all match up - all good :-D",
       "Gene names don't match! Danger!")

### Now filter SNPs by quality score, biallelic PASS status, and configured DP threshold

# CHANGED 2026-07: per-SNP depth filter is >=10x (was 50x), matching the amplicon
# coverage threshold in stage 1. Override with env var NANORAVE_MIN_DP.
min_snp_dp <- as.numeric(Sys.getenv("NANORAVE_MIN_DP", unset = "10"))
vcf_tbl_all_varinfo2 <- vcf_tbl_all_varinfo %>%
  filter(Filter=="PASS" & # quality pass
           Qual >= 1 & # Quality score >=1 (NOTE: Clair3 QUAL here maxes ~27, so this is effectively no quality filter)
           nchar(Ref) < 2 & # single SNP ref
           nchar(Alt) < 2 & # single SNP alt (ie biallelic)
           DP >= min_snp_dp) # per-SNP read depth (default 50x; was 10x)

dim(vcf_tbl_all_varinfo)
dim(vcf_tbl_all_varinfo2)


## Clean up columns
col_ord2 <- c("sample_id",
              "collection_site",
              "ont_seq_date",
              "ont_multiplex_group",
              "ont_barcode",
              "amplicon_id",
              "vcf_fn",
              "Gene",
              "Pos",
              "Ref",
              "Alt",
              "Qual",
              "Filter",
              "Pileup_call",
              "Alignment_call",
              "DP",
              "AF",
              "GT",
              "snp_id")
              
vcf_tbl_all_varinfo2 <- vcf_tbl_all_varinfo2 %>% select(all_of(col_ord2))

# ---------------------------------------------------------------------------
# COERCE THE NUMERIC FIELDS. vcfR returns the FORMAT fields as CHARACTER, and R
# then does the wrong thing quietly:
#
#   mean(c("0.65","0.66"))  -> NA        (silently!)
#   sd(c("0.65","0.66"))    -> 0.007     (works - it coerces internally)
#   "0.9091" >= 0.8         -> STRING comparison, not numeric
#
# The string comparison happens to give the right answer for AF in [0,1] purely
# because lexicographic order coincides with numeric order for same-format
# decimals - but it is wrong in principle, and it silently defeated the artefact
# rule (af_mean came out NA for every position, so nothing was ever flagged).
# ---------------------------------------------------------------------------
vcf_tbl_all_varinfo2 <- vcf_tbl_all_varinfo2 %>%
  mutate(across(any_of(c("AF", "DP", "Qual", "Pos")), ~ suppressWarnings(as.numeric(.x))))

stopifnot(is.numeric(vcf_tbl_all_varinfo2$AF), is.numeric(vcf_tbl_all_varinfo2$DP))

head(vcf_tbl_all_varinfo2,15)
dim(vcf_tbl_all_varinfo2)


### Call genotypes from the CONFIDENT HAPLOID call (CORRECTED 2026-06)
# Pf blood-stage is haploid; Clair3 was run with --haploid_precise. A real clonal
# mutation is emitted as homozygous-alt (GT "1" or "1/1") at AF ~1.0. The previous
# rule (AF > 0.51) miscalled recurrent ~0.5-AF heterozygous artifacts in low-complexity
# tracts (e.g. k13 1739 C580Y in a GT-repeat; mdr1 551 Y184F in a TA-repeat; dhps 1620
# K540N) as mutations. We therefore call a SNP "mutant" only on a confident hom-alt
# genotype, and flag 0/1 het calls separately as mixed/low-confidence (not counted as
# confident mutants). Override threshold via env NANORAVE_HOM_AF if needed.
hom_alt_gt <- c("1", "1/1", "1|1")
hom_af_min <- as.numeric(Sys.getenv("NANORAVE_HOM_AF", unset = "0.8"))

# ===========================================================================
# EXPLICIT THREE-CLASS VARIANT RULE  (replaces the implicit "hom-alt only" rule)
#
# Until now, artifact rejection happened only as a SIDE EFFECT of scoring 0/1
# calls as wild-type. That is fragile: the recurrent ~0.5-AF miscalls (k13 C580Y,
# mdr1 Y184F, dhps K540N/A581G) were still emitted by the caller and were one
# rule-change away from re-entering the results as a false "90/92 artemisinin
# resistant" finding. Make the rejection explicit, data-driven and auditable.
#
# The discriminator is the CROSS-SPECIMEN SPREAD of the allele fraction:
#   * a REAL polymorphism varies between mosquitoes - present in some, absent in
#     others, and clonal (AF ~0.9-1.0) where present;
#   * a SYSTEMATIC ERROR appears in nearly EVERY specimen at nearly the SAME
#     allele fraction (tight sd) - which no biological process produces.
#
# Anchored by a biological control: k13 C580Y is essentially absent from Africa
# (Girgis 2023: none; Nkemngo 2022: none of the WHO-validated list), yet is called
# in ~95% of specimens at AF 0.535 +/- 0.046. It CANNOT be real. That anchor
# validates the rule, which then flags mdr1 184 and dhps 540/581 by the same
# signature. The true C580Y in the KH2 control (GT=1, AF=0.985) is NOT flagged.
#
# Classes:
#   clonal   GT in {1,1/1} & AF >= hom_af_min          -> mutant (1)
#   mixed    GT 0/1 & AF in [mix_lo, mix_hi]           -> reported, not counted
#                                                          (majority-consensus, as Girgis)
#   artifact recurrent + intermediate AF + tight sd    -> NA (locus not callable)
#   ref      otherwise                                 -> wild type (0)
# ===========================================================================
# Recurrence threshold. Calibrated on the RUN23 field cohort (n=94): sweeping it
# from 0.80 down to 0.20 gives a STABLE flag set from 0.50 downward -- exactly the
# four known artefact positions (k13 1739, mdr1 551, dhps 1620, dhps 1742) -- and
# NO genuine variant is caught at any threshold, because the real ones sit at
# AF 0.84-0.95, above the 0.75 window. 0.40 sits mid-plateau with margin either
# side. The AF window and the tight sd do the discriminating; recurrence only
# separates a systematic error from a rare true variant (e.g. crt K76T, 4%).
art_min_frac <- as.numeric(Sys.getenv("NANORAVE_ARTEFACT_MIN_FRAC", unset = "0.40"))
art_af_lo    <- as.numeric(Sys.getenv("NANORAVE_ARTEFACT_AF_LO",   unset = "0.35"))
art_af_hi    <- as.numeric(Sys.getenv("NANORAVE_ARTEFACT_AF_HI",   unset = "0.75"))
art_max_sd   <- as.numeric(Sys.getenv("NANORAVE_ARTEFACT_MAX_SD",  unset = "0.10"))
art_min_n    <- as.numeric(Sys.getenv("NANORAVE_ARTEFACT_MIN_N",   unset = "5"))
mix_lo <- 0.25; mix_hi <- 0.75

# Controls are excluded from the artifact statistics: KH2's genuine C580Y (AF~1)
# would otherwise inflate the mean and sd of that position.
# Rows with a missing snp_id come from the outer join upstream and are not real
# positions; grouping them would create a phantom "NA" position and mask calls.
field_calls <- vcf_tbl_all_varinfo2 %>%
  filter(!nr_is_control(sample_id), !is.na(snp_id), !is.na(AF))
n_field_specimens <- dplyr::n_distinct(field_calls$sample_id)

artefact_stats <- field_calls %>%
  group_by(snp_id) %>%
  summarise(n_called  = dplyr::n(),
            af_mean   = mean(AF, na.rm = TRUE),
            af_sd     = stats::sd(AF, na.rm = TRUE),
            frac_het  = mean(GT %in% c("0/1", "0|1"), na.rm = TRUE),
            .groups   = "drop") %>%
  mutate(
    frac_specimens = n_called / max(n_field_specimens, 1),
    af_sd          = ifelse(is.na(af_sd), 0, af_sd),
    is_artefact    = n_field_specimens >= art_min_n &
                     frac_specimens >= art_min_frac &
                     af_mean >= art_af_lo & af_mean <= art_af_hi &
                     af_sd   <  art_max_sd
  )

artefact_ids <- artefact_stats$snp_id[artefact_stats$is_artefact]
artefact_ids <- artefact_ids[!is.na(artefact_ids)]

if (length(artefact_ids) > 0) {
  message("[stage2] SYSTEMATIC ARTEFACT positions detected (recurrent, intermediate AF, tight sd) - ",
          "these codons are marked NOT CALLABLE (NA), not wild-type: ",
          paste(artefact_ids, collapse = ", "))
} else {
  message("[stage2] no systematic-artefact positions detected (cohort n=", n_field_specimens, ").")
}

# Artefact catalogue - a first-class QC deliverable (Table T3 of the blueprint).
write.csv(artefact_stats %>% arrange(desc(frac_specimens)),
          file = file.path(analysis_dr_2, paste0(MinION_run_name, "_artefact_catalogue.csv")),
          row.names = FALSE)

vcf_tbl_all_varinfo2 <- vcf_tbl_all_varinfo2 %>%
  mutate(
    call_class = dplyr::case_when(
      !is.na(snp_id) & snp_id %in% artefact_ids                ~ "artefact",
      GT %in% hom_alt_gt & AF >= hom_af_min                    ~ "clonal",
      GT %in% c("0/1", "0|1") & AF >= mix_lo & AF <= mix_hi     ~ "mixed",
      TRUE                                                     ~ "ref"
    ),
    # Majority-consensus scoring (Girgis): only a confident clonal call is a mutant.
    # An artefact position is NOT callable -> NA (never silently scored wild-type).
    GT_majority = dplyr::case_when(
      call_class == "clonal"   ~ 1,
      call_class == "artefact" ~ NA_real_,
      TRUE                     ~ 0
    )
  )

message("[stage2] call classes: ",
        paste(names(table(vcf_tbl_all_varinfo2$call_class)),
              table(vcf_tbl_all_varinfo2$call_class), sep = "=", collapse = "  "))

## Add a column whether likely het, defined by a cutoff
het_cutoff = 0.9 # e.g. if using 90%, then non-ref AF >90%==hom-nref or <10%==hom-ref; in between == het
vcf_tbl_all_varinfo2$is_het <- ifelse( (vcf_tbl_all_varinfo2$AF < het_cutoff) &
                                         (vcf_tbl_all_varinfo2$AF > (1-het_cutoff) ),
                                       TRUE, FALSE)

summary(vcf_tbl_all_varinfo2$is_het) # Tells you how many hets

vcf_tbl_all_varinfo2 <- vcf_tbl_all_varinfo2 %>% # Rename GT_majority to genotype
  rename(genotype=GT_majority)

# Note - The AF (allelic frequency) is the percentage of the alternative allele.
# 1-AF is the reference allele. If 1-AF >0.5, the reference allele is the majority.



### Add in SNP data for the samples - note this is ignoring kelch13 which we deal with separately (below)
colnames(snp_data)

snp_data_cols <- c("snp_id",
                   "chrom",
                   "chrom_pos",
                   "pos_in_gene_coding",
                   "pos_in_codon",
                   "snp_ref",
                   "snp_alt",
                   "codon_pos_in_gene",
                   "codon_ref",
                   "codon_alt",
                   "aa_ref",
                   "aa_alt",
                   "aa_mut_name",
                   "key_snp")

snp_data_v2 <- snp_data %>%
  select(all_of(snp_data_cols)) %>%
  filter(key_snp == TRUE) %>%
  select(-key_snp)

# ---------------------------------------------------------------------------
# REFERENCE-POLARITY CORRECTION
#
# For most markers the 3D7 reference carries the WILD-TYPE allele, so
# "non-reference call" == "mutation present". That is NOT universally true.
#
# The 3D7 reference genuinely carries the RESISTANT allele at dhps codon 437
# (reference codon GGT = Gly = 437G; the marker table's own note says "present in
# 3D7 reference strain"). For that marker, a specimen matching the reference IS
# the mutant, and a non-reference call is the WILD TYPE. Scoring it like the
# others reported A437G = 0% when 437G is in fact at ~100% -- the exact opposite.
# (Girgis 2023 likewise notes the 3D7 dhps haplotype is SGKAA, i.e. 437G.)
#
# Rather than hard-code dhps 437, DERIVE the polarity: translate the reference
# codon from the marker table and compare it with the amino acid the table claims
# is the reference. Where they disagree, the reference carries the mutant and the
# genotype must be inverted. This self-corrects for any marker added later.
# ---------------------------------------------------------------------------
nr_codon_table <- c(
  TTT="F",TTC="F",TTA="L",TTG="L",CTT="L",CTC="L",CTA="L",CTG="L",
  ATT="I",ATC="I",ATA="I",ATG="M",GTT="V",GTC="V",GTA="V",GTG="V",
  TCT="S",TCC="S",TCA="S",TCG="S",CCT="P",CCC="P",CCA="P",CCG="P",
  ACT="T",ACC="T",ACA="T",ACG="T",GCT="A",GCC="A",GCA="A",GCG="A",
  TAT="Y",TAC="Y",CAT="H",CAC="H",CAA="Q",CAG="Q",AAT="N",AAC="N",
  AAA="K",AAG="K",GAT="D",GAC="D",GAA="E",GAG="E",TGT="C",TGC="C",
  TGG="W",CGT="R",CGC="R",CGA="R",CGG="R",AGT="S",AGC="S",AGA="R",
  AGG="R",GGT="G",GGC="G",GGA="G",GGG="G")

snp_data_v2 <- snp_data_v2 %>%
  mutate(
    codon_ref_aa   = unname(nr_codon_table[toupper(codon_ref)]),
    ref_is_mutant  = !is.na(codon_ref_aa) & (codon_ref_aa != aa_ref)
  )

inverted <- snp_data_v2 %>% filter(ref_is_mutant)
if (nrow(inverted) > 0) {
  message("[stage2] reference carries the MUTANT allele at ", nrow(inverted),
          " marker(s): ", paste(inverted$aa_mut_name, collapse = ", "),
          ". Raw genotypes are LEFT AS-IS here (the dhps haplotype builder below ",
          "depends on 'non-reference' semantics); the polarity is corrected when ",
          "prevalence is reported, in 03_resistance_frequencies.R.")
} else {
  warning("[stage2] no reference-carries-mutant marker detected. dhps A437G is ",
          "expected to be one (3D7 = SGKAA). Check the marker table.")
}

# IMPORTANT: `genotype` stays as "non-reference call" and is NOT inverted here.
# The dhps haplotype builder further down already encodes the correct polarity by
# hand (dhps_1310_GC == 0 -> "SGKAA" i.e. 437G; == 1 -> "SAKAA" i.e. 437A), so
# inverting the raw genotype would corrupt the haplotype -- the one output that is
# currently right. The polarity correction is applied at the REPORTING layer, in
# 03_resistance_frequencies.R, where the per-marker "A437G" prevalence is computed.
# `ref_is_mutant` is exported below so that stage 3 can do it.
write.csv(snp_data_v2 %>% select(snp_id, aa_mut_name, aa_ref, aa_alt,
                                 codon_ref, codon_ref_aa, ref_is_mutant),
          file = file.path(analysis_dr_2, paste0(MinION_run_name, "_marker_polarity.csv")),
          row.names = FALSE)

colnames(vcf_tbl_all_varinfo2)
sample_snp_calls_col <- c("sample_id",
                          "snp_id",
                          "genotype")
sample_snp_calls <- vcf_tbl_all_varinfo2 %>%
  select(all_of(sample_snp_calls_col))

### Make separate genotype calls for each sample in the batch

# Create integer for each group
sample_snp_calls <- sample_snp_calls %>%
  mutate(sample_integer = cumsum(!duplicated(sample_id)))


# loop through each sample integer to make list of genotype calls to add each one as a column to the key SNPs

maxn <- max(sample_snp_calls$sample_integer)

mylist_sample_vars <- list()

for (k in 1:maxn) {
  
  sub <- sample_snp_calls %>%
    filter(sample_integer==k) %>%
    select(c("snp_id", "sample_id", "genotype"))
  
  ncolname <- sub[1,"sample_id"]
  
  colnames(sub) <- c("snp_id","sample_id",ncolname)
  
  sub <- sub %>% select(c("snp_id", ncolname))
  
  mylist_sample_vars[[k]] <- sub
  
}

#mylist_sample_vars[[3]]
#length(mylist_sample_vars)

# get the list of SNPs to be included
snp_list <- snp_data_v2 %>%
  select("snp_id") %>%
  arrange(snp_id)
head(snp_list)


# Add the genotype calls to the snp_list for each sample, to standardise row number
  
for(m in 1:maxn) {
  mylist_sample_vars[[m]] <- left_join(snp_list, mylist_sample_vars[[m]], by="snp_id") %>%
    arrange(snp_id) %>% # make sure all in the same order so can cbind it all together
    select(-snp_id) # now drop the snp id columns to cbind all the genotypes together
}


#mylist_sample_vars[[60]]
#nrow(mylist_sample_vars[[1]]) # should be number of key_snps == 17

# Define a function to drop rows from a data frame if it has more than 17 rows
drop_rows_greater_than_17 <- function(mylist_sample_vars) {
  if (nrow(mylist_sample_vars) > 17) {
    mylist_sample_vars <- mylist_sample_vars[1:17, ]  # Keep only the first 17 rows
  }
  return(mylist_sample_vars)
}

# Apply the function to each data frame in mylist_sample_vars
mylist_sample_vars <- lapply(mylist_sample_vars, drop_rows_greater_than_17)





# cbind the genotype calls
keysnps_allsamps_gt_novars <- do.call(cbind, mylist_sample_vars)

# and cbind the snp info back in
keysnps_allsamps_gt <- cbind(snp_list, keysnps_allsamps_gt_novars)

## Replace NA with wild-type genotype call.
#
# Legitimate ONLY for the "no variant record" case: a specimen with no call at a
# codon matches the reference, i.e. wild type = 0.
#
# NOT legitimate for an ARTEFACT-masked codon. Those were deliberately set to NA
# ("not callable"), and a blanket NA -> 0 silently resurrects them as confident
# WILD-TYPE calls -- re-creating the very "missing == wild type" error this
# pipeline exists to avoid, and producing a bogus rule-out (e.g. "mdr1 Y184F =
# 0/447, 95% CI 0-0.8%") for a codon that was never measurable in the first place.
#
# So: fill NA -> 0 for the no-call case, then restore NA on the artefact codons.
keysnps_allsamps_gt <- keysnps_allsamps_gt %>% replace(is.na(.), 0)

if (length(artefact_ids) > 0) {
  artefact_rows <- keysnps_allsamps_gt$snp_id %in% artefact_ids
  if (any(artefact_rows)) {
    sample_cols <- setdiff(names(keysnps_allsamps_gt), "snp_id")
    keysnps_allsamps_gt[artefact_rows, sample_cols] <- NA
    message("[stage2] restored NA (not callable) on ", sum(artefact_rows),
            " artefact codon(s) after the wild-type fill: ",
            paste(keysnps_allsamps_gt$snp_id[artefact_rows], collapse = ", "))
  }
}


############ Check out the drug resistance (minus kelch13) genotype data!
keysnps_allsamps_gt





########### Deal with kelch13, csp and msp1
#thing <- rep(1:1000, each=3)
#write.csv(thing, file="codon_num_blank.csv", row.names=F)

#### Get the k13 variants

colnames(vcf_tbl_all_varinfo2)
col_ord3 <- c("sample_id",
              "collection_site",
              "amplicon_id",
              "snp_id",
              "Gene",
              "Pos",
              "Ref",
              "Alt",
              "Qual",
              "DP",
              "genotype")

# k13
# CORRECTED 2026-06: only confident hom-alt (genotype==1) k13 variants are eligible to
# flag artemisinin resistance. This rejects the recurrent het C580Y artifact (k13 1739
# in a GT microsatellite, AF~0.5) while retaining true hom-alt calls (e.g. KH2 control).
vcf_tbl_all_k13 <- vcf_tbl_all_varinfo2 %>%
  filter(Gene=="k13") %>%
  filter(genotype==1) %>%
  select(all_of(col_ord3)) %>%
  arrange(sample_id, snp_id)

# csp
vcf_tbl_all_csp <- vcf_tbl_all_varinfo2 %>%
  filter(Gene=="csp") %>%
  select(all_of(col_ord3)) %>%
  arrange(sample_id, snp_id)

# msp1
vcf_tbl_all_msp1 <- vcf_tbl_all_varinfo2 %>%
  filter(Gene=="msp1") %>%
  select(all_of(col_ord3)) %>%
  arrange(sample_id, snp_id)


##### Remove variants that are before or after my actual primers - unlikely to be real!!

vcf_tbl_all_k13$Pos <- as.numeric(vcf_tbl_all_k13$Pos)
vcf_tbl_all_csp$Pos <- as.numeric(vcf_tbl_all_csp$Pos)
vcf_tbl_all_msp1$Pos <- as.numeric(vcf_tbl_all_msp1$Pos)

vcf_tbl_all_k13 <- vcf_tbl_all_k13 %>%
  filter((Pos >= 1277) & (Pos <= 2145))

vcf_tbl_all_csp <- vcf_tbl_all_csp %>%
  filter((Pos >= 168) & (Pos <= 1143))  # Actual end in 3D7 is 1143

vcf_tbl_all_msp1 <- vcf_tbl_all_msp1 %>%
  filter((Pos >= 105) & (Pos <= 5099)) # Actual end in 3D7 is 5099




########### kelch13 - special analysis
### Is anything left?? Otherwise, will get errors:
ifelse(nrow(vcf_tbl_all_k13)==0,
       "No kelch13 mutations detected, no need to continue analysis",
       "At least 1 kelch13 mutation - proceed with analysis")

######### Check whether there are multiple kelch mutations affecting the same codon
# script I've written deals with each variant individually, would struggle with 2 SNPs in same codon

## Check whether there are 2 SNPs in same codon across whole dataset
k13_pos_subs <- vcf_tbl_all_k13 %>%
  select(snp_id, Pos) %>%
  rename(pos = Pos)
k13_pos_subs

# Add to the k13 seq data
k13_seq_var <- left_join(k13_annotated, k13_pos_subs, "pos") %>%
  arrange(snp_id) %>%
  drop_na %>%
  distinct(snp_id, .keep_all=TRUE)

temptest <- n_distinct(k13_seq_var$codon_num)

ifelse(temptest==nrow(k13_seq_var),
       "All of the SNPs in this dataset occur in unique codons, so it's fine",
       "Danger! Multiple SNPs within same codon")



###### If there are any kelch13 variants then proceed with this bit
if(nrow(vcf_tbl_all_k13)>0) {
  
  # Create integer for each sample
  vcf_tbl_all_k13 <- vcf_tbl_all_k13 %>%
    mutate(sample_integer = cumsum(!duplicated(sample_id)))
  
  
  k13nmax <- max(vcf_tbl_all_k13$sample_integer)
  
  mylist_k13_samps <- list()
  
  for (q in 1:k13nmax) {
    
    # subset for each sample
    vcf_tbl_k13_q <- vcf_tbl_all_k13 %>% filter(sample_integer==q)
    
    # Create integer for each SNP
    vcf_tbl_k13_q <- vcf_tbl_k13_q %>%
      mutate(snp_integer = cumsum(!duplicated(snp_id)))
    
    # start nested loop for each SNP
    k13nmaxsub <- max(vcf_tbl_k13_q$snp_integer)
    
    mylist_k13_snps <- list()
    
    for (r in 1:k13nmaxsub) {
      
      vcf_tbl_k13_q$Pos <- as.numeric(vcf_tbl_k13_q$Pos)
      vcf_tbl_k13_q$Ref <- as.character(vcf_tbl_k13_q$Ref)
      vcf_tbl_k13_q$Alt <- as.character(vcf_tbl_k13_q$Alt)
      
      ## Define position and nucleotide of each SNP in the loop
      query_pos = vcf_tbl_k13_q[r,'Pos']
      query_alt_nt = vcf_tbl_k13_q[r,'Alt']
      
      ## Extract original codon and mutant codon
      nt_ref <- k13_annotated[query_pos,3]
      nt_codon_pos <- k13_annotated[query_pos,4]
      
      ## Define the ref codon
      ref_codon <- ifelse(nt_codon_pos==1, paste(nt_ref,k13_annotated[query_pos+1,3],k13_annotated[query_pos+2,3],sep=""),
                          ifelse(nt_codon_pos==2, paste(k13_annotated[query_pos-1,3],nt_ref,k13_annotated[query_pos+1,3],sep=""),
                                 ifelse(nt_codon_pos==3, paste(k13_annotated[query_pos-2,3],k13_annotated[query_pos-1,3],nt_ref,sep=""),
                                        paste("Something is wrong")
                                 )
                          )
      )
      ## Define the alt codon
      alt_codon <- ifelse(nt_codon_pos==1, paste(query_alt_nt,k13_annotated[query_pos+1,3],k13_annotated[query_pos+2,3],sep=""),
                          ifelse(nt_codon_pos==2, paste(k13_annotated[query_pos-1,3],query_alt_nt,k13_annotated[query_pos+1,3],sep=""),
                                 ifelse(nt_codon_pos==3, paste(k13_annotated[query_pos-2,3],k13_annotated[query_pos-1,3],query_alt_nt,sep=""),
                                        paste("Something is wrong")
                                 )
                          )
      )
      
      ## Translate codons into amino acids
      ref_codon_translation_info <- genetic_code %>%
        filter(codon==ref_codon) %>%
        rename(ref_codon = codon, ref_aa_letter = aa_letter, ref_aa_code = aa_code, ref_aa_name = aa_name)
      
      alt_codon_translation_info <- genetic_code %>%
        filter(codon==alt_codon) %>%
        rename(alt_codon = codon, alt_aa_letter = aa_letter, alt_aa_code = aa_code, alt_aa_name = aa_name)
      
      ## combine codon info
      aa_change_info <- cbind(ref_codon_translation_info,alt_codon_translation_info)
      
      ## Define mutation type
      aa_change_info$mut_type <-
        ifelse(aa_change_info$ref_aa_letter==aa_change_info$alt_aa_letter, "Synonymous",
               ifelse( ((aa_change_info$ref_aa_letter!=aa_change_info$alt_aa_letter) &
                          (aa_change_info$alt_aa_letter!="Stop")), "Missense",
                       ifelse( ((aa_change_info$ref_aa_letter!=aa_change_info$alt_aa_letter) &
                                  (aa_change_info$alt_aa_letter=="Stop")), "Nonsense",
                               "Something has gone wrong")
               )
        )
      
      ## Define the codon number 
      aa_change_info$codon_count <- ifelse(query_pos %% 3 == 0, query_pos/3, round((query_pos/3)+0.5,0) )
      
      ## Populate other variant fields
      aa_change_info$gene <- "k13"
      aa_change_info$pos <- query_pos
      aa_change_info$ref_nt <- nt_ref
      aa_change_info$alt_nt <- query_alt_nt
      aa_change_info$pos_in_codon <- nt_codon_pos
      
      aa_change_info$snp_id <- paste(aa_change_info$gene, "_", aa_change_info$pos, "_", aa_change_info$ref_nt, aa_change_info$alt_nt, sep="")
      
      aa_change_info$aa_mut <- paste(aa_change_info$ref_aa_letter, aa_change_info$codon_count, aa_change_info$alt_aa_letter, sep="")
      aa_change_info$aa_mut_id <- paste(aa_change_info$gene, "_", aa_change_info$ref_aa_letter, aa_change_info$codon_count, aa_change_info$alt_aa_letter, sep="")
      
      aa_change_info$codon_count <- as.numeric(aa_change_info$codon_count)
      
      ## Is the mutation within the propeller domain/ WHO region to worry about?
      aa_change_info$WHO_ARTR_region <- ifelse(aa_change_info$codon_count >= 349, TRUE, FALSE)
      
      ## Classify against the Objective 3 WHO K13 marker catalogue.
      aa_change_info$k13_marker_classification <- k13_who_catalog$classification[
        match(aa_change_info$aa_mut, k13_who_catalog$aa_change)
      ]
      aa_change_info$k13_marker_classification <- ifelse(
        is.na(aa_change_info$k13_marker_classification),
        "not_listed",
        aa_change_info$k13_marker_classification
      )
      aa_change_info$mut_ARTR_phenotype <- aa_change_info$k13_marker_classification %in%
        c("validated", "candidate_or_associated")
      
      ## Keep track of the sample
      aa_change_info$sample_id <- vcf_tbl_k13_q[1,"sample_id"]
      
      ## Keep track of SNP quality and DP
      aa_change_info$Qual <- vcf_tbl_k13_q[r,"Qual"]
      aa_change_info$DP <- vcf_tbl_k13_q[r,"DP"]
      
      ## Organise the columns
      cols_order <- c("sample_id",
                      "snp_id",
                      "Qual",
                      "DP",
                      "gene",
                      "pos",
                      "ref_nt",
                      "alt_nt",
                      "pos_in_codon",
                      "ref_codon",
                      "alt_codon",
                      "codon_count",
                      "WHO_ARTR_region",
                      "mut_type",
                      "ref_aa_letter",
                      "alt_aa_letter",
                      "ref_aa_code",
                      "alt_aa_code",
	                      "ref_aa_name",
	                      "alt_aa_name",
	                      "aa_mut",
	                      "aa_mut_id",
	                      "k13_marker_classification",
	                      "mut_ARTR_phenotype")
      
      mutation_info <- aa_change_info %>% select(all_of(cols_order))
      
      mylist_k13_snps[[r]] <- mutation_info # Populate the list
      
    }
    
    k13_snps_df <- data.frame(do.call("rbind", mylist_k13_snps)) # Bind the list into a dataframe of SNPs for each sample
    
    mylist_k13_samps[[q]] <- k13_snps_df
    
  }
  
  k13_snps_allsamps_df <- data.frame(do.call("rbind", mylist_k13_samps)) # Bind the list into a df of samples

  k13_snps_allsamps_df

}

# A run can legitimately contain zero kelch13 variants (no ART-R sample, and no
# KH2 control on the flowcell). k13_snps_allsamps_df is only ever created inside
# the block above, but it is consumed unconditionally further down
# (k13_snps_who_artr), which crashed the stage with "object not found". Define an
# empty frame with the expected schema so the no-variant case flows through.
if (!exists("k13_snps_allsamps_df")) {
  k13_snps_allsamps_df <- data.frame(
    sample_id = character(), snp_id = character(), Qual = numeric(), DP = numeric(),
    gene = character(), pos = numeric(), ref_nt = character(), alt_nt = character(),
    pos_in_codon = numeric(), ref_codon = character(), alt_codon = character(),
    codon_count = numeric(), WHO_ARTR_region = logical(), mut_type = character(),
    ref_aa_letter = character(), alt_aa_letter = character(),
    ref_aa_code = character(), alt_aa_code = character(),
    ref_aa_name = character(), alt_aa_name = character(),
    aa_mut = character(), aa_mut_id = character(),
    k13_marker_classification = character(), mut_ARTR_phenotype = logical(),
    stringsAsFactors = FALSE
  )
  message("[stage2] No kelch13 variants in this run - continuing with an empty k13 table.")
}




########################## High level interpretations: data manipulations

####### Check out the data to be used
keysnps_allsamps_gt

####### transposed version
snp_calls <- keysnps_allsamps_gt[2:ncol(keysnps_allsamps_gt)]
var_cols <- c(keysnps_allsamps_gt$snp_id)
snp_calls_t <- snp_calls %>%
  tibble::rownames_to_column() %>%  
  pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value) 
names(snp_calls_t) <- c('sample_id', var_cols)


####### Create composite haplotype calls

## dhfr
# Options = IRNI (51I*-59R*-108N*-164I); or NRNI (51N-59R*-108N*-164I) or ICNI (51I*-59C-108N*-164I), or Other
# WT = NCSI
snp_calls_t2 <- within(snp_calls_t, {
  dhfr_haplotype <-
    ifelse( (dhfr_152_AT==1 &
               dhfr_175_TC==1 & 
               dhfr_323_GA==1 &
               dhfr_490_AT==0), "IRNI", 
            ifelse((dhfr_152_AT==0 &
                      dhfr_175_TC==1 & 
                      dhfr_323_GA==1 &
                      dhfr_490_AT==0 ), "NRNI",
                   ifelse((dhfr_152_AT==1 &
                             dhfr_175_TC==0 & 
                             dhfr_323_GA==1 &
                             dhfr_490_AT==0 ), "ICNI",
                          ifelse((dhfr_152_AT==0 &
                                    dhfr_175_TC==0 & 
                                    dhfr_323_GA==1 &
                                    dhfr_490_AT==0 ), "NCNI",
                                 ifelse((dhfr_152_AT==1 &
                                           dhfr_175_TC==1 & 
                                           dhfr_323_GA==0 &
                                           dhfr_490_AT==0 ), "IRSI",
                                        ifelse((dhfr_152_AT==0 &
                                                  dhfr_175_TC==1 & 
                                                  dhfr_323_GA==0 &
                                                  dhfr_490_AT==0 ), "NRSI",
                                               ifelse((dhfr_152_AT==1 &
                                                         dhfr_175_TC==0 & 
                                                         dhfr_323_GA==1 &
                                                         dhfr_490_AT==0 ), "ICRI",
                                                      ifelse((dhfr_152_AT==0 &
                                                                dhfr_175_TC==0 & 
                                                                dhfr_323_GA==1 &
                                                                dhfr_490_AT==0 ), "NCRI",
                                                             ifelse((dhfr_152_AT==0 &
                                                                       dhfr_175_TC==0 & 
                                                                       dhfr_323_GA==0 &
                                                                       dhfr_490_AT==0 ), "NCSI",
                                                                    ifelse((dhfr_152_AT==1 &
                                                                              dhfr_175_TC==1 & 
                                                                              dhfr_323_GA==1 &
                                                                              dhfr_490_AT==1 ), "IRNL", "Other")
                                                             )
                                                      )
                                               )
                                        )
                                 )
                          )
                   )
            )
    )
})

snp_calls_t2 %>% print(n=5, width=Inf)
snp_calls_t2 %>% count(dhfr_haplotype)
snp_calls_t2 %>% filter(dhfr_haplotype=="Other")


## dhps
# Options = AGKAA (436A*, 437G*, 540K, 581A, 613A); or AAKAA (436A*, 437A, 540K, 581A, 613A); or SGKAA (436S, 437G*, 540K, 581A, 613A); or SGEAA (436S, 437G*, 540E*, 581A, 613A)
# WT = SAKAA; 3D7 = SGKAA
#
# ---------------------------------------------------------------------------
# dhps codon 581 (dhps_1742_CG) is ARTEFACT-MASKED (NA) for the marker-prevalence
# table, because the position carries a systematic ~0.53-AF heterozygous miscall.
# But the 5-codon dhps haplotype NEEDS a value at 581, and leaving it NA would
# turn EVERY dhps haplotype into NA.
#
# Judgement call, stated explicitly: for HAPLOTYPE construction only, codon 581 is
# taken as REFERENCE (581A). Justification - Girgis 2023 reports A581G at 2.0% in
# Ghana and Amenga-Etego reports it as rare, so a genuine 581G in ~48% of specimens
# is not credible; the ALT signal is the artefact. Reference is the best estimate.
#
# The two outputs therefore say different, and both honest, things:
#   * haplotype table   -> SGKAA (581 assumed reference; flagged below)
#   * prevalence table  -> A581G = NOT CALLABLE (we cannot CONFIDENTLY RULE OUT
#                          low-frequency 581G at a corrupted position)
# ---------------------------------------------------------------------------
hap_artefact_codons <- intersect(artefact_ids, c("dhps_1742_CG", "dhps_1620_AT"))
if (length(hap_artefact_codons) > 0) {
  message("[stage2] haplotype construction: taking artefact codon(s) ",
          paste(hap_artefact_codons, collapse = ", "),
          " as REFERENCE (see comment). The marker table still reports them as NOT CALLABLE.")
}
snp_calls_t2$dhps_581_not_callable <- "dhps_1742_CG" %in% artefact_ids

# Remember which cells are artefact-NA, coalesce to reference FOR HAPLOTYPE
# CONSTRUCTION ONLY, and restore the NA afterwards (see restore block below).
# Without the restore, the coalesce would permanently overwrite the columns and the
# marker-prevalence table would report a bogus "0/413, 95% CI 0-0.9%" rule-out for
# codons that were never callable.
.hap_fill_cols <- intersect(c("dhps_1742_CG", "dhps_1620_AT"), names(snp_calls_t2))
.hap_na_mask <- lapply(snp_calls_t2[.hap_fill_cols], is.na)
snp_calls_t2 <- snp_calls_t2 %>%
  mutate(across(all_of(.hap_fill_cols), ~ dplyr::coalesce(.x, 0)))

snp_calls_t2 <- within(snp_calls_t2, {
  dhps_haplotype <-
    ifelse( (dhps_1306_TG==1 &
               dhps_1310_GC==0 & 
               dhps_1618_AG==0 &
               dhps_1742_CG==0 &
               dhps_1837_GT==0), "AGKAA", 
            ifelse( (dhps_1306_TG==1 &
                       dhps_1310_GC==1 & 
                       dhps_1618_AG==0 &
                       dhps_1742_CG==0 &
                       dhps_1837_GT==0), "AAKAA",
                    ifelse( (dhps_1306_TG==0 &
                               dhps_1310_GC==0 & 
                               dhps_1618_AG==0 &
                               dhps_1742_CG==0 &
                               dhps_1837_GT==0), "SGKAA", 
                            ifelse( (dhps_1306_TG==0 &
                                       dhps_1310_GC==1 & 
                                       dhps_1618_AG==0 &
                                       dhps_1742_CG==0 &
                                       dhps_1837_GT==0), "SAKAA",
                                    ifelse( (dhps_1306_TG==0 &
                                               dhps_1310_GC==1 & 
                                               dhps_1618_AG==0 &
                                               dhps_1742_CG==0 &
                                               dhps_1837_GT==1), "SAKAS",
                                            ifelse( (dhps_1306_TG==0 &
                                                       dhps_1310_GC==0 & 
                                                       dhps_1618_AG==0 &
                                                       dhps_1742_CG==0 &
                                                       dhps_1837_GT==1), "SGKAS",
                                                    ifelse( (dhps_1306_TG==1 &
                                                               dhps_1310_GC==0 & 
                                                               dhps_1618_AG==0 &
                                                               dhps_1742_CG==0 &
                                                               dhps_1837_GT==1), "AGKAS",
                                                            ifelse( (dhps_1306_TG==1 &
                                                                       dhps_1310_GC==1 & 
                                                                       dhps_1618_AG==0 &
                                                                       dhps_1742_CG==0 &
                                                                       dhps_1837_GT==1), "AAKAS",
                                                                    ifelse( (dhps_1306_TG==0 &
                                                                               dhps_1310_GC==0 & 
                                                                               dhps_1618_AG==0 &
                                                                               dhps_1742_CG==1 &
                                                                               dhps_1837_GT==0), "SGKGA", 
                                                                            ifelse( (dhps_1306_TG==1 &
                                                                                       dhps_1310_GC==0 & 
                                                                                       dhps_1618_AG==0 &
                                                                                       dhps_1742_CG==1 &
                                                                                       dhps_1837_GT==1), "AGKGS",   
                                                                            ifelse( (dhps_1306_TG==0 &
                                                                                       dhps_1310_GC==0 & 
                                                                                       dhps_1618_AG==0 &
                                                                                       dhps_1742_CG==1 &
                                                                                       dhps_1837_GT==1), "SGKGS",
                                                                            ifelse( (dhps_1306_TG==0 &
                                                                                       dhps_1310_GC==0 & 
                                                                                       dhps_1618_AG==1 &
                                                                                       dhps_1742_CG==0 &
                                                                                       dhps_1837_GT==0), "SGEAA", "Other")
                                                                    )
                                                            )
                                                    )
                                            )
                                    )
                            )
                    )
            ))
    ) )
})
snp_calls_t2 %>% print(n=5, width=Inf)
snp_calls_t2 %>% count(dhps_haplotype)
snp_calls_t2 %>% filter(dhps_haplotype=="Other")


## dhfr + dhps
# Options = dhfr-IRNI + dhps-AGKAA (436A*, 437G*, 540K, 581A, 613A); or AAKAA (436A*, 437A, 540K, 581A, 613A); or SGKAA (436S, 437G*, 540K, 581A, 613A)
snp_calls_t2 %>% count(dhfr_haplotype, dhps_haplotype)

snp_calls_t2 <- within(snp_calls_t2, {
  dhfr_dhps_haplotype <-
    paste0( "dhfr-", dhfr_haplotype, ", dhps-", dhps_haplotype)
})

snp_calls_t2 %>% print(n=23, width=Inf)
snp_calls_t2 %>% count(dhfr_dhps_haplotype) %>% arrange(desc(n))


########### GENETIC REPORT CARD

# Artemisinin = ART
# Chloroquine = CQ
# Pyrimethamine = PYR
# Sulfadoxine = SX
# https://www.wwarn.org/sites/default/files/drug-abbreviations.pdf

# Haplotypes are built. Restore the artefact NAs so the marker-prevalence table
# reports those codons as NOT CALLABLE rather than as a confident 0%.
for (.cc in .hap_fill_cols) {
  snp_calls_t2[[.cc]][.hap_na_mask[[.cc]]] <- NA
}
if (length(.hap_fill_cols) > 0) {
  message("[stage2] restored artefact NA on ", paste(.hap_fill_cols, collapse = ", "),
          " after haplotype construction (marker table will show NOT CALLABLE).")
}

snp_calls_t3 <- within(snp_calls_t2, {
  
  CQ <- ifelse(crt_227_AC==1,"R",
               ifelse(crt_227_AC==0,"S","Undetermined"))
  
  PYR <- ifelse(dhfr_323_GA==1,"R",
                ifelse(dhfr_323_GA==0,"S","Undetermined"))
  
  SX <- ifelse(dhps_1310_GC==0,"R",
               ifelse(dhps_1310_GC==1,"S","Undetermined"))
  
  SP.Rx <- ifelse( ((dhfr_152_AT==1) & (dhfr_175_TC==1) & (dhfr_323_GA==1)) ,"R",
                   ifelse(((dhfr_152_AT==0) | (dhfr_175_TC==0) | (dhfr_323_GA==0)) ,"S",
                          "Undetermined")
  )
  
  SP.IPTp <- ifelse( ((dhfr_152_AT==1) & (dhfr_175_TC==1) & (dhfr_323_GA==1)) &
                       ((dhps_1310_GC==0) & (dhps_1618_AG==1)) &
                       ((dhfr_490_AT==1) | (dhps_1742_CG==1) | (dhps_1837_GT==1) | (dhps_1837_GA==1)), "R",
                     "S")
})

# Sort out ART-R
# Match controls on the BASE id: pooled cohorts disambiguate repeated control
# sample_ids as e.g. KH2__LUC_DRAG1_IMPAVESRUN23_20250920, so an exact-string
# match would no longer recognise them and the C580Y-carrying control would leak
# into the list of ART-R field samples.
k13_snps_who_artr <- k13_snps_allsamps_df %>%
  filter(mut_ARTR_phenotype %in% TRUE) %>%
  filter(!nr_is_control(sample_id))

k13_resistant_samples <- unique(k13_snps_who_artr$sample_id)
snp_calls_t3$ART <- ifelse(snp_calls_t3$sample_id %in% k13_resistant_samples, "R", "S")



# Check it out
snp_calls_t3 %>% print(n=5, width=Inf)


############### create genotype format for csp

### Limit to CTD
vcf_tbl_csp_ctd <- vcf_tbl_all_csp %>%
  filter((Pos >= 816) & (Pos <= 1143))  # end of central repeat region to primer


# Create integer for each group
vcf_tbl_csp_ctd <- vcf_tbl_csp_ctd %>%
  mutate(sample_integer = cumsum(!duplicated(sample_id)))


# loop through each sample integer to make list of genotype calls to add each one as a column to the key SNPs

maxn <- max(vcf_tbl_csp_ctd$sample_integer)

mylist_csp_vars <- list()

for (t in 1:maxn) {
  
  sub <- vcf_tbl_csp_ctd %>%
    filter(sample_integer==t) %>%
    select(c("snp_id", "sample_id", "genotype"))
  
  ncolname <- sub[1,"sample_id"]
  
  colnames(sub) <- c("snp_id","sample_id",ncolname)
  
  sub <- sub %>% select(c("snp_id", ncolname))
  
  mylist_csp_vars[[t]] <- sub
  
}

# (removed a bare `mylist_csp_vars[[3]]` debug print: it hard-indexed the third
# element and crashed with "subscript out of bounds" whenever fewer than three
# samples carried csp variants -- harmless on a full flowcell, fatal on a small
# cohort.)
length(mylist_csp_vars)


## Join everything together
csp_samp_calls <- reduce(mylist_csp_vars, full_join, by = "snp_id")

## replace NA with 0
csp_samp_calls <- csp_samp_calls %>%
  replace(is.na(.), 0) %>%
  arrange(snp_id)

head(csp_samp_calls)



############### Save outputs

# Core drug resistance SNP info in various formats
write.csv(vcf_varcount_df4, file=out_variant_info_fn, row.names=FALSE)
write.csv(vcf_tbl_all_varinfo2, file=out_variant_info2_fn, row.names=FALSE)
write.csv(keysnps_allsamps_gt, file=out_keysnp_allsamp_genotypes_fn, row.names=FALSE)
write.csv(snp_calls_t3, file=out_keysnp_allsamp_genotypes_samprows_fn, row.names=FALSE)

# kelch13 variant info. Always written - an empty table with headers is a valid
# result ("no k13 variants found") and keeps the output layout consistent, whereas
# a missing file is indistinguishable from a crashed stage.
write.csv(k13_snps_allsamps_df, file=out_k13_calls_info_fn, row.names=FALSE)

# csp nref genotypes
write.csv(vcf_tbl_all_csp, file=out_csp_genotypes_long_fn, row.names=FALSE)
write.csv(csp_samp_calls, file=out_csp_genotypes_samps_fn, row.names=FALSE)
