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

############### R script for combined analysis of drug resistance genotypes for VB + DBS samples
### Clair3

##################### R script to read in BED files to check coverage

### Load packages
library("tidyverse")
library("readxl")
library("cowplot")

# 02_genotypes block: frequency/haplotype tables -> frequencies/, GRC plots -> plots/
# Resolved through the path registry (_pipeline_paths.R) rather than hardcoded
# subfolder names, so renaming a block only needs editing in one place.
plots_dir <- nr_paths$geno_plots
tables_dir <- nr_paths$geno_freq
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

### Output filenames
plot_dhfr_dhps_haplo_cnts_fn <- file.path(plots_dir, paste0("dhfr_dhps_haplotype_counts_summary_qc_",analysis_date,".png"))
plot_GRC_cnts_fn <- file.path(plots_dir, paste0("GRC_counts_summary_qc_",analysis_date,".png"))

plot_dhfr_dhps_haplo_rat_fn <- file.path(plots_dir, paste0("dhfr_dhps_haplotype_ratio_summary_qc_",analysis_date,".png"))
plot_GRC_rat_fn <- file.path(plots_dir, paste0("GRC_ratio_summary_qc_",analysis_date,".png"))

plot_combined_ratios_fn <- file.path(plots_dir, paste0("dhfr_dhps_GRC_combined_ratio_summary_qc_",analysis_date,".png"))

plot_combined_ratios_vb_fn <- file.path(plots_dir, paste0("dhfr_dhps_GRC_combined_ratio_summary_qc_vb_",analysis_date,".png"))
plot_combined_ratios_dbs_fn <- file.path(plots_dir, paste0("dhfr_dhps_GRC_combined_ratio_summary_qc_dbs_",analysis_date,".png"))
plot_combined_ratios_vb_dbs_fn <- file.path(plots_dir, paste0("dhfr_dhps_GRC_combined_ratio_summary_qc_vb_dbs_",analysis_date,".png"))

table_dr_genotypes_fn <- file.path(tables_dir, paste0("drug_resistance_variants_vb_dbs_all_",analysis_date, "_", MinION_run_name,".csv"))
table_dr_freq_fn <- file.path(tables_dir, paste0("drug_resistance_frequencies_vb_dbs_all_",analysis_date,"_",MinION_run_name,".csv"))
table_dhfr_freq_fn <- file.path(tables_dir, paste0("dhfr_haplo_frequencies_",analysis_date,"_",MinION_run_name, ".csv"))
table_dhps_freq_fn <- file.path(tables_dir, paste0("dhps_haplo_frequencies_",analysis_date,"_",MinION_run_name,".csv"))
table_dhfr_dhps_freq_fn <- file.path(tables_dir, paste0("dhfr_dhps_haplo_frequencies_",analysis_date,"_",MinION_run_name,".csv"))

### Import data
## Coverage for each sample within each MinION run
# define paths
cov_B2_raw <- read.csv(file.path(analysis_dr, paste0(MinION_run_name, "_coverage_by_run_sample", analysis_date, ".csv")))

## Genotypes for each sample within each MinION run
# define paths
geno_B2_raw <- read.csv(file.path(analysis_dr_2, paste0(MinION_run_name, "_genotype_calls_haplotypes_samplerows_v2.csv")))


### Step 1 - Mask variants below the configured locus coverage threshold
### Canonical analysis uses 50x; sensitivity configs can override this, e.g. 10x.

# Merge the coverage data into the genotype data
geno_cov_B2 <- merge(geno_B2_raw,cov_B2_raw,
                    by="sample_id",
                    all.x=T,all.y=T)
geno_cov_B2
dim(geno_cov_B2)


### Replace genotypes with NA if coverage too low for that amplicon for that run

## Run A
# crt
geno_cov_B2_v2 <- within(geno_cov_B2, {
  crt_227_AC = ifelse(crt_coverage_above_threshold=="FALSE",NA,crt_227_AC)
  
  # dhfr
  dhfr_152_AT = ifelse(dhfr_coverage_above_threshold=="FALSE",NA,dhfr_152_AT)
  dhfr_175_TC = ifelse(dhfr_coverage_above_threshold=="FALSE",NA,dhfr_175_TC)
  dhfr_323_GA = ifelse(dhfr_coverage_above_threshold=="FALSE",NA,dhfr_323_GA)
  dhfr_323_GC = ifelse(dhfr_coverage_above_threshold=="FALSE",NA,dhfr_323_GC)
  dhfr_490_AT = ifelse(dhfr_coverage_above_threshold=="FALSE",NA,dhfr_490_AT)
  dhfr_492_AG = ifelse(dhfr_coverage_above_threshold=="FALSE",NA,dhfr_492_AG)
  
  # dhps
  dhps_1306_TG = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_1306_TG)
  dhps_1307_CT = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_1307_CT)
  dhps_1310_GC = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_1310_GC)
  dhps_1618_AG = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_1618_AG)
#  dhps_1620_AG = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_1620_AG)
  dhps_1620_AT = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_1620_AT)
  dhps_1742_CG = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_1742_CG)
  dhps_1837_GA = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_1837_GA)
  dhps_1837_GT = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_1837_GT)
  
  # mdr1
  mdr1_256_AT = ifelse(mdr1_coverage_above_threshold=="FALSE",NA,mdr1_256_AT)
  mdr1_551_AT = ifelse(mdr1_coverage_above_threshold=="FALSE",NA,mdr1_551_AT)
  
  # dhfr_haplotype
  dhfr_haplotype = ifelse(dhfr_coverage_above_threshold=="FALSE",NA,dhfr_haplotype)
  
  # dhps_haplotype
  dhps_haplotype = ifelse(dhps_coverage_above_threshold=="FALSE",NA,dhps_haplotype)
  
  # dhps-dhfr haplotype
  dhfr_dhps_haplotype = ifelse((dhfr_coverage_above_threshold=="FALSE" | 
                                  dhps_coverage_above_threshold=="FALSE"),NA,dhfr_dhps_haplotype)
  
  # GRC
  CQ = ifelse(crt_coverage_above_threshold=="FALSE",NA,CQ)
  PYR = ifelse(dhfr_coverage_above_threshold=="FALSE",NA,PYR)
  SX = ifelse(dhps_coverage_above_threshold=="FALSE",NA,SX)
  SP.Rx = ifelse((dhfr_coverage_above_threshold=="FALSE" |
                    dhps_coverage_above_threshold=="FALSE"),NA,SP.Rx)
  SP.IPTp = ifelse((dhfr_coverage_above_threshold=="FALSE" |
                      dhps_coverage_above_threshold=="FALSE"),NA,SP.IPTp)
  ART = ifelse(k13_coverage_above_threshold=="FALSE",NA,ART)
})


## Specify minion_run
geno_cov_B2_v2$minion_run <- "B"


### Merge all the runs together
dim(geno_cov_B2_v2)
colnames(geno_cov_B2_v2)

#geno_cov_combined <- rbind(geno_cov_A_v2,
#                           geno_cov_B2_v2,
#                           geno_cov_C_v2,
#                           geno_cov_D_v2,
#                           geno_cov_E_v2,
#                           geno_cov_F_v2)

geno_cov_combined <- geno_cov_B2_v2

dim(geno_cov_combined)


### Check sample counts are as expected


### Fail the sample overall if has crt + either dhfr or dhps low cov
geno_cov_combined <- within(geno_cov_combined, {
  sample_qc <- ifelse((crt_coverage_above_threshold=="FALSE" &
                         (dhfr_coverage_above_threshold=="FALSE" | 
                            dhps_coverage_above_threshold=="FALSE")), "FAIL", "PASS")
})



## Shed all the coverage stats
colnames(geno_cov_combined)
cols <- c("sample_id",
          "minion_run",
          "sample_qc",
          "crt_227_AC",
          "dhfr_152_AT",
          "dhfr_175_TC",
          "dhfr_323_GA",
          "dhfr_323_GC",
          "dhfr_490_AT",
          "dhfr_492_AG",
          "dhps_1306_TG",
          "dhps_1307_CT",
          "dhps_1310_GC",
          "dhps_1618_AG",
          "dhps_1620_AT",
          "dhps_1742_CG",
          "dhps_1837_GA",
          "dhps_1837_GT",
          "mdr1_256_AT",
          "mdr1_551_AT",
          "dhfr_haplotype",
          "dhps_haplotype",
          "dhfr_dhps_haplotype",
          "CQ",
          "PYR",
          "SX",
          "SP.Rx",
          "SP.IPTp",
          "ART")
geno_cov_combined_v2 <- geno_cov_combined %>% select(all_of(cols))

head(geno_cov_combined_v2)

geno_cov_combined_v2 %>% count(sample_qc)

unique(geno_cov_combined_v2$sample_id)


### Highlight controls
# Matched on the BASE sample id: in a pooled cohort a control that was run on every
# flowcell is disambiguated as e.g. KH2__LUC_DRAG1_IMPAVESRUN23_20250920, so an
# exact-string match would silently stop recognising it and the controls would be
# counted as field samples in every frequency denominator.
control_ids <- c("Control_HB3", "Control_Dd2", "Control_KH2", "control_KH2",
                 "control_Dd2", "KH2", "NC", "PC")
geno_cov_combined_v2 <- within(geno_cov_combined_v2, {
  control <- ifelse(nr_base_sample_id(sample_id) %in% control_ids, "TRUE", "FALSE")
})
geno_cov_combined_v2 %>% count(control)


geno_qc_controls <- geno_cov_combined_v2 %>%
  filter(control=="TRUE")

geno_qc_samples_nocontrols <- geno_cov_combined_v2 %>%
  filter(control=="FALSE")

### RED-1 FIX (2026-06-26): the frequency cohort is the FULL field specimen set,
### NOT a composite sample_qc-PASS subset. Per-locus NA-masking (above) already
### removes non-callable codons, so each marker's denominator is its own callable
### count (mdr1 93 / crt 86 / dhfr 12 / dhps 52 on n_total=94). The old
### `filter(sample_qc == "PASS")` silently dropped ~8 specimens, deflating
### n_total to 86 and forcing mdr1 onto /86 instead of /93. The `sample_qc`
### column is retained (computed above) for the attrition narrative only and is
### never used as a frequency denominator. This makes the stage-3 table reconcile
### exactly with denominator_audit_marker_table.csv (the single source of truth).
geno_qc_pass_nocontrols <- geno_qc_samples_nocontrols   # full field cohort (no QC composite filter)

dim(geno_qc_pass_nocontrols)

### Remove sample duplicates
geno_qc_pass_nocontrols <- geno_qc_pass_nocontrols %>%
  arrange(minion_run) %>%
  distinct(sample_id, .keep_all=TRUE) %>%
  filter(sample_id != "SPT43954_v2") # remove the DBS sample sequenced twice

head(geno_qc_pass_nocontrols)
dim(geno_qc_pass_nocontrols)

## Add metadata for location
colnames(metadata)
meta_sub <- metadata %>%
  select(sample_id, collection_site)

geno_meta <- merge(geno_qc_pass_nocontrols, meta_sub,
                   by="sample_id",
                   all.x=T, all.y=F)
head(geno_meta)
dim(geno_meta)

geno_meta <- geno_meta %>%
  mutate(ART = str_replace(ART, ">=1 mutation detected", "R"))

geno_meta_all <- geno_meta

##### Overall drug resistance marker frequencies
geno_qc_pass_nocontrols %>% count(crt_227_AC) # Chloroquine resistance, K76T

geno_qc_pass_nocontrols %>% count(dhfr_323_GA) # Pyrimethamine resistance, S108N

geno_qc_pass_nocontrols %>% count(dhfr_490_AT) # Highly pyrimethamine resistance, I164L

geno_qc_pass_nocontrols %>% count(dhps_1310_GC) # Sulfadoxine resistance, A437G - 3D7 is resistant!

geno_qc_pass_nocontrols %>% count(dhps_1618_AG) # High-level Sulfadoxine resistance, K540E
geno_qc_pass_nocontrols %>% count(dhps_1742_CG) # High-level Sulfadoxine resistance, A581G
geno_qc_pass_nocontrols %>% count(dhps_1837_GT) # High-level Sulfadoxine resistance, A613S

geno_qc_pass_nocontrols %>% count(SP.IPTp) # Any resistance to IPTp


########## Compute DHFR + DHPS haplotype frequencies
### basic counts for DHFR haplotypes
geno_qc_pass_nocontrols %>% count(dhfr_haplotype) %>% arrange(desc(n)) %>% rename(count = n)

### basic counts for DHPS haplotypes
geno_qc_pass_nocontrols %>% count(dhps_haplotype) %>% arrange(desc(n)) %>% rename(count = n)

### basic counts for DHFR + DHPS haplotype combinations
geno_qc_pass_nocontrols %>% count(dhfr_dhps_haplotype) %>% arrange(desc(n)) %>% rename(count = n)

### basic counts for MDR1
geno_qc_pass_nocontrols %>% count(mdr1_256_AT) %>% arrange(desc(n)) %>% rename(count = n)
geno_qc_pass_nocontrols %>% count(mdr1_551_AT) %>% arrange(desc(n)) %>% rename(count = n)


##### Format drug resistance frequencies
## CRT K76T
geno_meta_all %>%
  count(crt_227_AC)
samp_count <- nrow(geno_meta_all %>% filter(crt_227_AC==1))
total_n <- nrow(geno_meta_all)
crt_76_all_nref_count <- samp_count
crt_76_all_nref_freq <- samp_count/total_n
crt_76_all_nref_pcnt <- signif(crt_76_all_nref_freq*100, 3)
crt_76_all_df <- data.frame(Nref_count = crt_76_all_nref_count,
                            Nref_freq = crt_76_all_nref_freq,
                            Nref_pcnt = crt_76_all_nref_pcnt)
crt_76_all_df$gene = "CRT"
crt_76_all_df$SNP = "crt_227_AC"
crt_76_all_df$Mutation = "K76T"
crt_76_all_df <- crt_76_all_df %>% select(gene, SNP, Mutation, Nref_count, Nref_freq, Nref_pcnt)


# DHFR_N51I
geno_meta_all %>%
  count(dhfr_152_AT)
samp_count <- nrow(geno_meta_all %>% filter(dhfr_152_AT==1))
total_n <- nrow(geno_meta_all)
dhfr_51_all_nref_count <- samp_count
dhfr_51_all_nref_freq <- samp_count/total_n
dhfr_51_all_nref_pcnt <- signif(dhfr_51_all_nref_freq*100, 3)
dhfr_51_all_df <- data.frame(Nref_count = dhfr_51_all_nref_count,
                             Nref_freq = dhfr_51_all_nref_freq,
                             Nref_pcnt = dhfr_51_all_nref_pcnt)
dhfr_51_all_df$gene = "DHFR"
dhfr_51_all_df$SNP = "dhfr_152_AT"
dhfr_51_all_df$Mutation = "N51I"
dhfr_51_all_df <- dhfr_51_all_df %>% select(gene, SNP, Mutation, Nref_count, Nref_freq, Nref_pcnt)


# DHFR_C59R
geno_meta_all %>%
  count(dhfr_175_TC)
samp_count <- nrow(geno_meta_all %>% filter(dhfr_175_TC==1))
total_n <- nrow(geno_meta_all)
dhfr_59_all_nref_count <- samp_count
dhfr_59_all_nref_freq <- samp_count/total_n
dhfr_59_all_nref_pcnt <- signif(dhfr_59_all_nref_freq*100, 3)
dhfr_59_all_df <- data.frame(Nref_count = dhfr_59_all_nref_count,
                             Nref_freq = dhfr_59_all_nref_freq,
                             Nref_pcnt = dhfr_59_all_nref_pcnt)
dhfr_59_all_df$gene = "DHFR"
dhfr_59_all_df$SNP = "dhfr_175_TC"
dhfr_59_all_df$Mutation = "C59R"
dhfr_59_all_df <- dhfr_59_all_df %>% select(gene, SNP, Mutation, Nref_count, Nref_freq, Nref_pcnt)


# DHFR_S108N
geno_meta_all %>%
  count(dhfr_323_GA)
samp_count <- nrow(geno_meta_all %>% filter(dhfr_323_GA==1))
total_n <- nrow(geno_meta_all)
dhfr_108_all_nref_count <- samp_count
dhfr_108_all_nref_freq <- samp_count/total_n
dhfr_108_all_nref_pcnt <- signif(dhfr_108_all_nref_freq*100, 3)
dhfr_108_all_df <- data.frame(Nref_count = dhfr_108_all_nref_count,
                              Nref_freq = dhfr_108_all_nref_freq,
                              Nref_pcnt = dhfr_108_all_nref_pcnt)
dhfr_108_all_df$gene = "DHFR"
dhfr_108_all_df$SNP = "dhfr_323_GA"
dhfr_108_all_df$Mutation = "S108N"
dhfr_108_all_df <- dhfr_108_all_df %>% select(gene, SNP, Mutation, Nref_count, Nref_freq, Nref_pcnt)


## DHFR_S108T - invariant, all Ref
geno_meta_all %>% count(dhfr_323_GC)

dhfr_108T_nil <- data.frame(gene = "DHFR",
                             SNP = "dhfr_323_GC",
                             Mutation = "S108T",
                             Nref_count = 0,
                             Nref_freq = 0,
                             Nref_pcnt = 0)


## DHFR_I164L
geno_meta_all %>% count(dhfr_490_AT)

samp_count <- nrow(geno_meta_all %>% filter(dhfr_490_AT==1))
total_n <- nrow(geno_meta_all)
dhfr_164_all_nref_count <- samp_count
dhfr_164_all_nref_freq <- samp_count/total_n
dhfr_164_all_nref_pcnt <- signif(dhfr_164_all_nref_freq*100, 3)
dhfr_164_nil <- data.frame(gene = "DHFR",
                           SNP = "dhfr_490_AT",
                           Mutation = "I164L",
                           Nref_count = dhfr_164_all_nref_count,
                           Nref_freq = dhfr_164_all_nref_freq,
                           Nref_pcnt = dhfr_164_all_nref_pcnt)


# DHFR_I164M - invariant, all Ref
geno_meta_all %>% count(dhfr_492_AG)

dhfr_164M_nil <- data.frame(gene = "DHFR",
                             SNP = "dhfr_492_AG",
                             Mutation = "I164M",
                             Nref_count = 0,
                             Nref_freq = 0,
                             Nref_pcnt = 0)

## DHPS_S436A
geno_meta_all %>%
  count(dhps_1306_TG)
samp_count <- nrow(geno_meta_all %>% filter(dhps_1306_TG==1))
total_n <- nrow(geno_meta_all)
dhps_436_all_nref_count <- samp_count
dhps_436_all_nref_freq <- samp_count/total_n
dhps_436_all_nref_pcnt <- signif(dhps_436_all_nref_freq*100, 3)
dhps_436_all_df <- data.frame(Nref_count = dhps_436_all_nref_count,
                              Nref_freq = dhps_436_all_nref_freq,
                              Nref_pcnt = dhps_436_all_nref_pcnt)
dhps_436_all_df$gene = "DHPS"
dhps_436_all_df$SNP = "dhps_1306_TG"
dhps_436_all_df$Mutation = "S436A"
dhps_436_all_df <- dhps_436_all_df %>% select(gene, SNP, Mutation, Nref_count, Nref_freq, Nref_pcnt)


## DHPS_S436F
geno_meta_all %>%
  count(dhps_1307_CT)
samp_count <- nrow(geno_meta_all %>% filter(dhps_1307_CT==1))
total_n <- nrow(geno_meta_all)
dhps_436F_all_nref_count <- samp_count
dhps_436F_all_nref_freq <- samp_count/total_n
dhps_436F_all_nref_pcnt <- signif(dhps_436F_all_nref_freq*100, 3)
dhps_436F_all_df <- data.frame(Nref_count = dhps_436F_all_nref_count,
                               Nref_freq = dhps_436F_all_nref_freq,
                               Nref_pcnt = dhps_436F_all_nref_pcnt)
dhps_436F_all_df$gene = "DHPS"
dhps_436F_all_df$SNP = "dhps_1307_CT"
dhps_436F_all_df$Mutation = "S436F"
dhps_436F_all_df <- dhps_436F_all_df %>% select(gene, SNP, Mutation, Nref_count, Nref_freq, Nref_pcnt)


## DHPS_A437G - note that is the 'resistant' allele in 3D7!
geno_meta_all %>%
  count(dhps_1310_GC)
samp_count <- nrow(geno_meta_all %>% filter(dhps_1310_GC==1))
total_n <- nrow(geno_meta_all)
dhps_437_all_nref_count <- samp_count
dhps_437_all_nref_freq <- samp_count/total_n
dhps_437_all_nref_pcnt <- signif(dhps_437_all_nref_freq*100, 3)
dhps_437_all_df <- data.frame(Nref_count = dhps_437_all_nref_count,
                              Nref_freq = dhps_437_all_nref_freq,
                              Nref_pcnt = dhps_437_all_nref_pcnt)
dhps_437_all_df$gene = "DHPS"
dhps_437_all_df$SNP = "dhps_1310_GC"
dhps_437_all_df$Mutation = "A437G"
dhps_437_all_df <- dhps_437_all_df %>% select(gene, SNP, Mutation, Nref_count, Nref_freq, Nref_pcnt)


# DHPS_K540E
geno_meta_all %>% count(dhps_1618_AG)

samp_count <- nrow(geno_meta_all %>% filter(dhps_1618_AG==1))
total_n <- nrow(geno_meta_all)
dhps_540E_all_nref_count <- samp_count
dhps_540E_all_nref_freq <- samp_count/total_n
dhps_540E_all_nref_pcnt <- signif(dhps_540E_all_nref_freq*100, 3)

dhps_540E_nil <- data.frame(gene = "DHPS",
                            SNP = "dhps_1618_AG",
                            Mutation = "K540E",
                            Nref_count = dhps_540E_all_nref_count,
                            Nref_freq = dhps_540E_all_nref_freq,
                            Nref_pcnt = dhps_540E_all_nref_pcnt)

# DHPS_K540N
geno_meta_all %>% count(dhps_1620_AT)

samp_count <- nrow(geno_meta_all %>% filter(dhps_1620_AT==1))
total_n <- nrow(geno_meta_all)
dhps_540N_all_nref_count <- samp_count
dhps_540N_all_nref_freq <- samp_count/total_n
dhps_540N_all_nref_pcnt <- signif(dhps_540N_all_nref_freq*100, 3)

dhps_540N_nil <- data.frame(gene = "DHPS",
                            SNP = "dhps_1620_AT",
                            Mutation = "K540N",
                            Nref_count = dhps_540N_all_nref_count,
                            Nref_freq = dhps_540N_all_nref_freq,
                            Nref_pcnt = dhps_540N_all_nref_pcnt)

# DHPS_A581G
geno_meta_all %>%
  count(dhps_1742_CG)
samp_count <- nrow(geno_meta_all %>% filter(dhps_1742_CG==1))
total_n <- nrow(geno_meta_all)
dhps_581_all_nref_count <- samp_count
dhps_581_all_nref_freq <- samp_count/total_n
dhps_581_all_nref_pcnt <- signif(dhps_581_all_nref_freq*100, 3)
dhps_581_all_df <- data.frame(Nref_count = dhps_581_all_nref_count,
                              Nref_freq = dhps_581_all_nref_freq,
                              Nref_pcnt = dhps_581_all_nref_pcnt)
dhps_581_all_df$gene = "DHPS"
dhps_581_all_df$SNP = "dhps_1742_CG"
dhps_581_all_df$Mutation = "A581G"
dhps_581_all_df <- dhps_581_all_df %>% select(gene, SNP, Mutation, Nref_count, Nref_freq, Nref_pcnt)


# DHPS_A613S
geno_meta_all %>%
  count(dhps_1837_GT)
samp_count <- nrow(geno_meta_all %>% filter(dhps_1837_GT==1))
total_n <- nrow(geno_meta_all)
dhps_613S_all_nref_count <- samp_count
dhps_613S_all_nref_freq <- samp_count/total_n
dhps_613S_all_nref_pcnt <- signif(dhps_613S_all_nref_freq*100, 3)
dhps_613S_all_df <- data.frame(Nref_count = dhps_613S_all_nref_count,
                               Nref_freq = dhps_613S_all_nref_freq,
                               Nref_pcnt = dhps_613S_all_nref_pcnt)
dhps_613S_all_df$gene = "DHPS"
dhps_613S_all_df$SNP = "dhps_1837_GT"
dhps_613S_all_df$Mutation = "A613S"
dhps_613S_all_df <- dhps_613S_all_df %>% select(gene, SNP, Mutation, Nref_count, Nref_freq, Nref_pcnt)


# DHPS_A613T - invariant, all ref
geno_meta_all %>% count(dhps_1837_GA)

dhps_613T_nil <- data.frame(gene = "DHPS",
                             SNP = "dhps_1837_GA",
                             Mutation = "A613T",
                             Nref_count = 0,
                             Nref_freq = 0,
                             Nref_pcnt = 0)

# MDR1_N86Y
geno_meta_all %>%
  count(mdr1_256_AT)
samp_count <- nrow(geno_meta_all %>% filter(mdr1_256_AT==1))
total_n <- nrow(geno_meta_all)
mdr1_86_all_nref_count <- samp_count
mdr1_86_all_nref_freq <- samp_count/total_n
mdr1_86_all_nref_pcnt <- signif(mdr1_86_all_nref_freq*100, 3)
mdr1_86_all_df <- data.frame(Nref_count = mdr1_86_all_nref_count,
                             Nref_freq = mdr1_86_all_nref_freq,
                             Nref_pcnt = mdr1_86_all_nref_pcnt)
mdr1_86_all_df$gene = "MDR1"
mdr1_86_all_df$SNP = "mdr1_256_AT"
mdr1_86_all_df$Mutation = "N86Y"
mdr1_86_all_df <- mdr1_86_all_df %>% select(gene, SNP, Mutation, Nref_count, Nref_freq, Nref_pcnt)


# MDR1_Y184F
geno_meta_all %>%
  count(mdr1_551_AT)
samp_count <- nrow(geno_meta_all %>% filter(mdr1_551_AT==1))
total_n <- nrow(geno_meta_all)
mdr1_184_all_nref_count <- samp_count
mdr1_184_all_nref_freq <- samp_count/total_n
mdr1_184_all_nref_pcnt <- signif(mdr1_184_all_nref_freq*100, 3)
mdr1_184_all_df <- data.frame(Nref_count = mdr1_184_all_nref_count,
                              Nref_freq = mdr1_184_all_nref_freq,
                              Nref_pcnt = mdr1_184_all_nref_pcnt)
mdr1_184_all_df$gene = "MDR1"
mdr1_184_all_df$SNP = "mdr1_551_AT"
mdr1_184_all_df$Mutation = "Y184F"
mdr1_184_all_df <- mdr1_184_all_df %>% select(gene, SNP, Mutation, Nref_count, Nref_freq, Nref_pcnt)


### Compile
dr_freq_table <- rbind(crt_76_all_df,
                       dhfr_51_all_df,
                       dhfr_59_all_df,
                       dhfr_108_all_df,
                       dhfr_108T_nil,
                       dhfr_164_nil,
                       dhfr_164M_nil,
                       dhps_436_all_df,
                       dhps_436F_all_df,
                       dhps_437_all_df,
                       dhps_540E_nil,
                       dhps_540N_nil,
                       dhps_581_all_df,
                       dhps_613S_all_df,
                       dhps_613T_nil,
                       mdr1_86_all_df,
                       mdr1_184_all_df)
dr_freq_table

### CORRECTED 2026-06: locus-callable denominators (not whole-cohort n).
# The blocks above divided Nref_count by nrow(geno_meta_all) (all PASS samples),
# which deflates poorly-covered loci (e.g. dhfr). Coverage-based NA is applied per
# amplicon, so the callable denominator is per gene. Recompute freq/pcnt on the
# callable n and add an exact binomial 95% CI. This matches the proposal's
# "locus-level denominators among successfully genotyped specimens".
# PER-SNP callable denominator (was: one proxy SNP per gene).
#
# The gene-level proxy is unsafe now that systematic-artefact codons are masked to
# NA. MDR1's proxy was mdr1_551_AT = Y184F, which IS the artefact: masking it would
# have driven the whole MDR1 denominator to zero and silently killed N86Y with it.
# Callability is a property of the CODON (amplicon coverage AND not-an-artefact),
# so count it per SNP.
dr_freq_table$n_callable <- vapply(
  dr_freq_table$SNP,
  function(s) if (s %in% names(geno_meta_all)) sum(!is.na(geno_meta_all[[s]])) else NA_integer_,
  integer(1))

# A codon masked as a systematic artefact is NOT CALLABLE - it is not "0% mutant".
dr_freq_table$callable <- ifelse(is.na(dr_freq_table$n_callable) | dr_freq_table$n_callable == 0,
                                 "not_callable", "callable")
if (any(dr_freq_table$callable == "not_callable")) {
  message("[stage3] markers reported as NOT CALLABLE (artefact-masked or no coverage): ",
          paste(unique(dr_freq_table$Mutation[dr_freq_table$callable == "not_callable"]),
                collapse = ", "))
}

# If NOTHING is callable at ANY locus, that is not a result - it means the join
# between coverage, genotypes and the sample sheet failed upstream. Without this
# guard the stage writes a full table of n_callable = 0 / freq = NA and the
# pipeline exits 0, which is how a silent column-ordering bug shipped a complete
# set of empty results.
if (all(is.na(dr_freq_table$n_callable)) || sum(dr_freq_table$n_callable, na.rm = TRUE) == 0) {
  stop("[stage3] every marker has n_callable = 0. Coverage/genotype/sample-sheet ",
       "joins have failed upstream - refusing to write an empty frequency table. ",
       "Check that stage 1's coverage table has real sample_ids.")
}
dr_freq_table$n_total    <- nrow(geno_meta_all)

# ---------------------------------------------------------------------------
# REFERENCE-POLARITY CORRECTION  (see 02_genotype_calls.R for the derivation)
#
# `Nref_count` counts specimens carrying a NON-REFERENCE allele. For almost every
# marker the 3D7 reference is wild-type, so non-reference == mutation present.
#
# dhps A437G is the exception: the 3D7 reference genuinely CARRIES the resistant
# allele (reference codon GGT = Gly = 437G; Girgis 2023: "3D7 = SGKAA"). There,
# `Nref_count` counts WILD-TYPE (437A) specimens, so reporting it under the label
# "A437G" states the exact opposite of the truth -- it reported A437G = 0% when
# 437G is in fact at ~100%.
#
# Carriage of the NAMED mutation = n_callable - Nref_count for such markers.
# The flag is derived (not hard-coded) in stage 2 and exported to
# <run>_marker_polarity.csv; we re-derive it here so this stage is self-contained.
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

# stage 3 does not otherwise read the marker table; load it here.
snp_data_polarity <- readxl::read_excel(file.path(resource_dir, "DR_variant_info_v2.xlsx"), sheet = 1)

polarity <- snp_data_polarity %>%
  filter(key_snp == TRUE) %>%
  transmute(SNP = snp_id,
            ref_is_mutant = !is.na(nr_codon_table[toupper(codon_ref)]) &
                            unname(nr_codon_table[toupper(codon_ref)]) != aa_ref)

dr_freq_table <- dr_freq_table %>%
  left_join(polarity, by = "SNP") %>%
  mutate(ref_is_mutant = ifelse(is.na(ref_is_mutant), FALSE, ref_is_mutant))

if (any(dr_freq_table$ref_is_mutant)) {
  flipped <- dr_freq_table$Mutation[dr_freq_table$ref_is_mutant]
  message("[stage3] reference carries the mutant allele at: ",
          paste(unique(flipped), collapse = ", "),
          " - reporting CARRIAGE of the named mutation (n_callable - Nref_count).")
  dr_freq_table$Nref_count <- ifelse(
    dr_freq_table$ref_is_mutant,
    dr_freq_table$n_callable - dr_freq_table$Nref_count,
    dr_freq_table$Nref_count)
} else {
  warning("[stage3] no reference-carries-mutant marker found; dhps A437G is ",
          "expected to be one. Check DR_variant_info_v2.xlsx.")
}

# Rename to say what it now means: carriage of the named mutation, not "non-reference".
dr_freq_table$mutant_count <- dr_freq_table$Nref_count

dr_freq_table$Nref_freq  <- ifelse(dr_freq_table$n_callable > 0,
                                    dr_freq_table$Nref_count / dr_freq_table$n_callable, NA_real_)
dr_freq_table$Nref_pcnt  <- signif(dr_freq_table$Nref_freq * 100, 3)
dr_freq_table$Nref_pcnt_95CI <- mapply(function(x, n) {
  if (is.finite(n) && n > 0) {
    bt <- binom.test(x, n)$conf.int
    paste0(signif(100 * bt[1], 3), "-", signif(100 * bt[2], 3))
  } else NA_character_
}, dr_freq_table$Nref_count, dr_freq_table$n_callable)
dr_freq_table <- dr_freq_table %>%
  select(gene, SNP, Mutation, callable, Nref_count, n_callable, n_total,
         Nref_freq, Nref_pcnt, Nref_pcnt_95CI)
dr_freq_table


############ DATA VIZ

##### Plot 1 - dhfr and dhps haplotypes
### Data manipulations

##### Summary haplotype counts
threshold_haplotype_count <- 2
max_plot_y_count = 25

### DHFR haplotype summary
dhfr_haplotype_summary1 <- geno_qc_pass_nocontrols %>%
  count(dhfr_haplotype) %>%
  arrange(desc(n)) %>%
  rename(hap_count = n)
total_n <- sum(dhfr_haplotype_summary1$hap_count)
dhfr_haplotype_summary1$hap_freq <- dhfr_haplotype_summary1$hap_count/total_n
dhfr_haplotype_summary1$hap_pcnt <- signif(dhfr_haplotype_summary1$hap_freq*100, 3)

dhfr_haplotype_summary <- dhfr_haplotype_summary1 %>% select(-hap_freq, -hap_pcnt)
dhfr_haplotype_summary$dhfr_haplotype_summary <- ifelse(dhfr_haplotype_summary$hap_count > threshold_haplotype_count,
                                                        dhfr_haplotype_summary$dhfr_haplotype,
                                                        "Other")
sub <- dhfr_haplotype_summary %>%
  filter(dhfr_haplotype_summary=="Other")
sub_n <- sum(sub$hap_count)

dhfr_haplotype_other <- data.frame(dhfr_haplotype_summary = c("Other"),
                                   hap_count = c(sub_n))

dhfr_haplotype_summary2 <- dhfr_haplotype_summary %>%
  filter(dhfr_haplotype_summary != "Other") %>%
  select(-dhfr_haplotype)

dhfr_haplotype_summary3 <- rbind(dhfr_haplotype_summary2, dhfr_haplotype_other) %>%
  rename(dhfr_haplotype = dhfr_haplotype_summary)

dhfr_haplotype_summary3$hap_freq <- dhfr_haplotype_summary3$hap_count/total_n
dhfr_haplotype_summary3$hap_pcnt <- signif(dhfr_haplotype_summary3$hap_freq*100, 3)
dhfr_haplotype_summary3

### DHPS haplotype summary
dhps_haplotype_summary1 <- geno_qc_pass_nocontrols %>%
  count(dhps_haplotype) %>%
  arrange(desc(n)) %>%
  rename(hap_count = n)
total_n <- sum(dhps_haplotype_summary1$hap_count)
dhps_haplotype_summary1$hap_freq <- dhps_haplotype_summary1$hap_count/total_n
dhps_haplotype_summary1$hap_pcnt <- signif(dhps_haplotype_summary1$hap_freq*100, 3)

dhps_haplotype_summary <- dhps_haplotype_summary1 %>% select(-hap_freq, -hap_pcnt)
dhps_haplotype_summary$dhps_haplotype_summary <- ifelse(dhps_haplotype_summary$hap_count > threshold_haplotype_count,
                                                        dhps_haplotype_summary$dhps_haplotype,
                                                        "Other")
sub <- dhps_haplotype_summary %>%
  filter(dhps_haplotype_summary=="Other")
sub_n <- sum(sub$hap_count)

dhps_haplotype_other <- data.frame(dhps_haplotype_summary = c("Other"),
                                   hap_count = c(sub_n))

dhps_haplotype_summary2 <- dhps_haplotype_summary %>%
  filter(dhps_haplotype_summary != "Other") %>%
  select(-dhps_haplotype)

dhps_haplotype_summary3 <- rbind(dhps_haplotype_summary2, dhps_haplotype_other) %>%
  rename(dhps_haplotype = dhps_haplotype_summary)

dhps_haplotype_summary3$hap_freq <- dhps_haplotype_summary3$hap_count/total_n
dhps_haplotype_summary3$hap_pcnt <- signif(dhps_haplotype_summary3$hap_freq*100, 3)
dhps_haplotype_summary3


### DHFR + DHPS combined haplotype summary
dhfr_dhps_hap_freq <- geno_qc_pass_nocontrols %>%
  count(dhfr_dhps_haplotype) %>%
  arrange(desc(n)) %>%
  rename(hap_count = n) 
total_n <- sum(dhfr_dhps_hap_freq$hap_count)
dhfr_dhps_hap_freq$hap_freq <- dhfr_dhps_hap_freq$hap_count/total_n
dhfr_dhps_hap_freq$hap_pcnt <- signif(dhfr_dhps_hap_freq$hap_freq*100, 3)
dhfr_dhps_hap_freq

dhfr_dhps_hap_freq_v2 <- dhfr_dhps_hap_freq
dhfr_dhps_hap_freq_v2$dhfr_dhps_haplotype <- ifelse(dhfr_dhps_hap_freq_v2$hap_count > threshold_haplotype_count,
                                                    dhfr_dhps_hap_freq_v2$dhfr_dhps_haplotype,
                                                    "Other")

dhfr_dhps_hap_freq_v2_sub <- dhfr_dhps_hap_freq_v2 %>%
  filter(dhfr_dhps_haplotype=="Other")
dhfr_dhps_hap_count_other <- sum(dhfr_dhps_hap_freq_v2_sub$hap_count)

dhfr_dhps_hap_count_other_tab <- data.frame(dhfr_dhps_haplotype = c("Other"),
                                            hap_count = c(dhfr_dhps_hap_count_other))

dhfr_dhps_hap_freq_v3 <- dhfr_dhps_hap_freq_v2 %>%
  filter(dhfr_dhps_haplotype != "Other") %>%
  select(-hap_freq, -hap_pcnt)

dhfr_dhps_hap_freq_short <- rbind(dhfr_dhps_hap_freq_v3, dhfr_dhps_hap_count_other_tab)

dhfr_dhps_hap_freq_short$hap_freq <- dhfr_dhps_hap_freq_short$hap_count/total_n
dhfr_dhps_hap_freq_short$hap_pcnt <- signif(dhfr_dhps_hap_freq_short$hap_freq*100, 3)
dhfr_dhps_hap_freq_short

dhfr_haplotype_summary3
dhps_haplotype_summary3
dhfr_dhps_hap_freq_short



##### Prepare for plotting
colnames(geno_qc_pass_nocontrols)

geno_qc_pass_nocontrols$crt_227_AC <- as.character(geno_qc_pass_nocontrols$crt_227_AC)
geno_qc_pass_nocontrols$mdr1_256_AT <- as.character(geno_qc_pass_nocontrols$mdr1_256_AT)
geno_qc_pass_nocontrols$mdr1_551_AT <- as.character(geno_qc_pass_nocontrols$mdr1_551_AT)

cols_geno_plot <- c("sample_id",
                    "dhfr_haplotype",
                    "dhps_haplotype",
                    "dhfr_dhps_haplotype")

cols_GRC_plot <- c("sample_id",
                   "CQ",
                   "SX",
                   "PYR",
                   "SP.Rx",
                   "SP.IPTp",
                   "ART")

data_geno_plot <- geno_qc_pass_nocontrols %>%
  select(all_of(cols_geno_plot)) %>%
  rename(DHFR = dhfr_haplotype,
         DHPS = dhps_haplotype,
         DHFR.DHPS = dhfr_dhps_haplotype) %>%
  mutate(DHFR.DHPS = str_replace(DHFR.DHPS, "dhfr-IRNI, dhps-AGKAA", "IRNI.AGKAA")) %>%
  mutate(DHFR.DHPS = str_replace(DHFR.DHPS, "dhfr-IRNI, dhps-SGKAA", "IRNI.SGKAA")) %>%
  mutate(DHFR.DHPS = str_replace(DHFR.DHPS, "dhfr-IRNI, dhps-AGKAS", "IRNI.AGKAS"))
data_geno_plot$DHFR.DHPS <- ifelse( (data_geno_plot$DHFR.DHPS == "IRNI.AGKAA" |
                                       data_geno_plot$DHFR.DHPS == "IRNI.SGKAA" |
                                       data_geno_plot$DHFR.DHPS == "IRNI.AGKAS" ),
                                    data_geno_plot$DHFR.DHPS, "Other")

data_geno_plot$DHFR <- ifelse( (data_geno_plot$DHFR == "IRNI" |
                                  data_geno_plot$DHFR == "NRNI" |
                                  data_geno_plot$DHFR == "NCSI" ),
                               data_geno_plot$DHFR, "Other")

data_geno_plot$DHPS <- ifelse( (data_geno_plot$DHPS == "AGKAA" |
                                  data_geno_plot$DHPS == "SGKAA" |
                                  data_geno_plot$DHPS == "AGKAS" |
                                  data_geno_plot$DHPS == "AAKAA"),
                               data_geno_plot$DHPS, "Other")

data_GRC_plot <- geno_qc_pass_nocontrols %>%
  select(all_of(cols_GRC_plot))


##### Genotypes
data_geno_plot2 <- data_geno_plot %>%
  gather(Genotype, Status, -c(sample_id)) %>% # go from wide to long format
  group_by(Genotype) %>%
  count(Status) %>%
  rename(Count = n)
data_geno_plot2 <- drop_na(data_geno_plot2)
head(data_geno_plot2,10)

### Control plot variables
unique(data_geno_plot2$Status)

legend_labels <- c("IRNI", # dhfr triple mutant
                   "NRNI", # dhfr double mutant
                   "NCSI",
                   "AGKAA", # dhps double mutant
                   "SGKAA", # dhps single mutant
                   "AGKAS", # dhps triple mutant
                   "AAKAA",
                   "IRNI.AGKAA", # dhfr double + dhps double
                   "IRNI.SGKAA", # dhfr triple + dhps single
                   "IRNI.AGKAS",
                   "Other")

break_ord <- c("IRNI",
               "NRNI",
               "NCSI",
               "AGKAA",
               "SGKAA",
               "AGKAS",
               "AAKAA",
               "IRNI.AGKAA",
               "IRNI.SGKAA",
               "IRNI.AGKAS",
               "Other")

order <- c("DHFR",
           "DHPS",
           "DHFR.DHPS")

# Colour-blind friendly palette (up to 8):
n_distinct(data_geno_plot2$Status) # Need 11 colours

#https://thenode.biologists.com/data-visualization-with-flying-colors/research/
#https://zenodo.org/record/3381072
#col_geno <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#d0d0d0")
#col_geno <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#0072B2", "#CC79A7", "#d5d5d5")
#col_geno <- c("#E69F00", "#AA3377", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#0072B2", "#CC79A7", '#EE6677', "#000000", "#d5d5d5")

col_geno <- c("#56B4E9", "#AA3377", "#E69F00", "#117733", "#F0E442", "#D55E00", "#000000", "#0072B2", '#EE6677', "#44AA99", "#d5d5d5")


### Generate the plots

unique(data_geno_plot2$Status)
n_distinct(data_geno_plot2$Status)

## plot - counts
data_geno_plot2$Status = factor(data_geno_plot2$Status, levels=break_ord)
p1_1 <- data_geno_plot2 %>%
  ggplot(aes(x=Genotype, y=Count, fill=Status)) +
  scale_x_discrete(limits = order) +
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("") +
  ylab("Count") + 
  ylim(0,max_plot_y_count) +
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=30, hjust=1),
        axis.title=element_text(size=11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        panel.background = element_rect(fill = 'white', colour = NA),
        plot.background = element_rect(fill="white", colour = NA)) +
  scale_fill_manual(values = col_geno, name = "Genotype", labels = legend_labels, breaks = break_ord)
p1_1

## plot - ratios
p1_2 <- data_geno_plot2 %>%
  ggplot(aes(x=Genotype, y=Count, fill=Status)) +
  scale_x_discrete(limits = order) +
  geom_bar(position="fill", stat="identity") +
  theme_minimal() +
  xlab("") +
  ylab("Proportion") + 
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=30, hjust=1),
        axis.title=element_text(size=11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        panel.background = element_rect(fill = 'white', colour = NA),
        plot.background = element_rect(fill="white", colour = NA)) +
  scale_fill_manual(values = col_geno, name = "Genotype", labels = legend_labels, breaks = break_ord)
p1_2

## plot - ratios - V2 for combined plot
p1_2_v2 <- data_geno_plot2 %>%
  ggplot(aes(x=Genotype, y=Count, fill=Status)) +
  scale_x_discrete(limits = order) +
  geom_bar(position="fill", stat="identity") +
  theme_minimal() +
  xlab("") +
  ylab("Proportion") + 
  theme(axis.text=element_text(size=8),
        axis.text.x = element_text(angle=30, hjust=1),
        axis.title=element_text(size=9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        panel.background = element_rect(fill = 'white', colour = NA),
        plot.background = element_rect(fill="white", colour = NA)) +
  scale_fill_manual(values = col_geno, name = "Genotype", labels = legend_labels, breaks = break_ord)
p1_2_v2



##### Drug resistance predictions

data_GRC_plot <- data_GRC_plot %>%
  mutate(ART = str_replace(ART, ">=1 mutation detected", "R"))

data_GRC_plot2 <- data_GRC_plot %>%
  gather(Drug, Status, -c(sample_id)) %>% # go from wide to long format
  group_by(Drug) %>%
  count(Status) %>%
  rename(Count = n)

data_GRC_plot2 <- drop_na(data_GRC_plot2)

head(data_GRC_plot2,10)

### Control plot variables
legend_labels_dr <- c("Susceptible",
                      "Resistant")

break_ord_dr <- c("S",
                  "R")

order_dr <- c("CQ",
              "SX",
              "PYR",
              "SP.Rx",
              "SP.IPTp",
              "ART")

cols_dr <- c("#0072B2", # Blue; Sensitive
             "#F0E442") # Yellow; Resistant
# If need an 'Unknown' then use Green, #009E73

## plot - counts
p2_1 <- data_GRC_plot2 %>%
  ggplot(aes(x=Drug, y=Count, fill=Status)) +
  scale_x_discrete(limits = order_dr) +
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("") +
  ylab("Count") + 
  ylim(0,max_plot_y_count) +
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=30, hjust=1),
        axis.title=element_text(size=11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        panel.background = element_rect(fill = 'white', colour = NA),
        plot.background = element_rect(fill="white", colour = NA)) +
  scale_fill_manual(values = cols_dr, name = "Status", labels = legend_labels_dr, breaks = break_ord_dr)
p2_1

## plot - ratios
p2_2 <- data_GRC_plot2 %>%
  ggplot(aes(x=Drug, y=Count, fill=Status)) +
  scale_x_discrete(limits = order_dr) +
  geom_bar(stat="identity", position="fill") +
  theme_minimal() +
  xlab("") +
  ylab("Proportion") + 
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(angle=30, hjust=1),
        axis.title=element_text(size=11),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        panel.background = element_rect(fill = 'white', colour = NA),
        plot.background = element_rect(fill="white", colour = NA)) +
  scale_fill_manual(values = cols_dr, name = "Status", labels = legend_labels_dr, breaks = break_ord_dr)
p2_2

## plot - ratios; ALT for combined plot, no y-axis
p2_2_v2 <- data_GRC_plot2 %>%
  ggplot(aes(x=Drug, y=Count, fill=Status)) +
  scale_x_discrete(limits = order_dr) +
  geom_bar(stat="identity", position="fill") +
  theme_minimal() +
  xlab("") +
  ylab("") + 
  theme(axis.text.x=element_text(size=8, angle=30, hjust=1),
        axis.text.y=element_text(size=8),
        axis.title=element_text(size=9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        panel.background = element_rect(fill = 'white', colour = NA),
        plot.background = element_rect(fill="white", colour = NA)) +
  scale_fill_manual(values = cols_dr, name = "Status", labels = legend_labels_dr, breaks = break_ord_dr)
p2_2_v2


###### Combined plot
plot_combine <- plot_grid(p1_2_v2,
                          p2_2_v2,
                          ncol = 2,
                          rel_widths = c(1, 1.3),
                          labels = c('A', 'B'), 
                          label_size = 10)
plot_combine


######### Save outputs
ggsave(plot=p1_1, file=plot_dhfr_dhps_haplo_cnts_fn, dpi=200, height=4, width=7, units="in")
ggsave(plot=p2_1, file=plot_GRC_cnts_fn, dpi=200, height=4, width=7, units="in")

ggsave(plot=p1_2, file=plot_dhfr_dhps_haplo_rat_fn, dpi=200, height=4, width=7, units="in")
ggsave(plot=p2_2, file=plot_GRC_rat_fn, dpi=200, height=4, width=7, units="in")

ggsave(plot=plot_combine, file=plot_combined_ratios_fn, dpi=200, height=4, width=8, units="in")

write.csv(geno_meta_all, table_dr_genotypes_fn, row.names = F)
write.csv(dr_freq_table, table_dr_freq_fn, row.names = F)
write.csv(dhfr_haplotype_summary1, table_dhfr_freq_fn, row.names = F)
write.csv(dhps_haplotype_summary1, table_dhps_freq_fn, row.names = F)
write.csv(dhfr_dhps_hap_freq, table_dhfr_dhps_freq_fn, row.names = F)
