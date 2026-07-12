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

##### Load packages
library("tidyverse")
library("readxl")

### Output file names
overall_cov_fn <- file.path(
  analysis_dr,
  paste0(MinION_run_name, "_coverage_amplicon_run_summary_", analysis_date, ".csv")
)
cov_by_sample_fn <- file.path(
  analysis_dr,
  paste0(MinION_run_name, "_coverage_by_run_sample", analysis_date, ".csv")
)


### Import data
## Identify file paths
# find BED files
bedfile_dir <- file.path(input_run_dir, "genome_coverage")

# amplicon position data
amp_pos_raw_fn <- file.path(resource_dir, "amplicon_gene_positions.csv")

if (!dir.exists(bedfile_dir)) {
  stop("Coverage bedGraph directory not found: ", bedfile_dir)
}

if (!file.exists(amp_pos_raw_fn)) {
  stop("Amplicon position file not found: ", amp_pos_raw_fn)
}

# Nanopore sequence run metadata
## Read in the data
amp_pos_raw <- read.csv(amp_pos_raw_fn) 

### Set coverage threshold
# CORRECTED 2026-06: minimum median amplicon coverage required to call markers.
# CHANGED 2026-07: 10x (was 50x). The 50x bar came from Girgis et al. (2023) and
# was used for the feasibility phase, which has now been answered; at 50x, dhfr
# was callable in only 27.7% of samples. Override with env var NANORAVE_MIN_COV.
min_cov_threshold <- as.numeric(Sys.getenv("NANORAVE_MIN_COV", unset = "10"))

### Manipulations
## Prepare amplicon position data
amp_pos_crt <- amp_pos_raw %>%
  filter(gene=="crt")

amp_pos_dhfr <- amp_pos_raw %>%
  filter(gene=="dhfr")

amp_pos_dhps <- amp_pos_raw %>%
  filter(gene=="dhps")

amp_pos_mdr1 <- amp_pos_raw %>%
  filter(gene=="mdr1")

amp_pos_kelch13 <- amp_pos_raw %>%
  filter(gene=="kelch13")

amp_pos_csp <- amp_pos_raw %>%
  filter(gene=="csp")

#amp_pos_msp1 <- amp_pos_raw %>%
#  filter(gene=="msp1")

## ONT metadata
#samp_names <- meta %>% select(ont_barcode, sample_id, patient_id, ont_multiplex_group, ont_seq_date)

samp_names <- meta %>%
  select(ont_barcode, sample_id, patient_id, ont_multiplex_group, ont_seq_date) 


########## CRT
# read in bed file
filenames <- list.files(bedfile_dir, pattern="_crt.bedGraph", full.names=TRUE)
if (length(filenames) == 0) {
  stop("No CRT bedGraph files found in: ", bedfile_dir)
}
ldf <- lapply(filenames, read.table)

#barcode = str_sub(filenames, 142, 150)
barcode <- str_extract(filenames, "barcode\\d+")
run_of_file <- nr_run_from_file(filenames)   # barcode01 exists on every flowcell

max_count <- length(ldf)


### Filter variants to only include those within the amplicon target region
list_crt <- list()
for (i in 1:max_count) {
  
  ldf[[i]] <- ldf[[i]] %>%
    rename(gene_id = V1,
           start = V2,
           end = V3,
           depth = V4)
  
  ldf[[i]]$start <- as.numeric(ldf[[i]]$start)
  ldf[[i]]$end <- as.numeric(ldf[[i]]$end)
  ldf[[i]]$depth <- as.numeric(ldf[[i]]$depth)
  
  ldf[[i]] <- ldf[[i]] %>%
    filter(start > amp_pos_crt$amplicon_start) %>%
    filter(end < amp_pos_crt$amplicon_end)
  
  list_crt[[i]] <- ldf[[i]]
  
}

### Now extract coverage statistics for each amplicon
list_crt_cov <- list()

for (j in 1:max_count) {
  
  summ <- summary(list_crt[[j]]$depth)
  summdf <- data.frame(summ=matrix(summ, ncol=6))
  colnames(summdf) <- names(summ)
  
  ## Reproduce barcode ID
  id <- barcode[j]
  
  summdf$ont_barcode <- id
  summdf$ont_multiplex_group <- run_of_file[j]
  summdf$gene_target <- "crt"
  
  summdf2 <- merge(summdf, samp_names, by=c("ont_multiplex_group","ont_barcode"),
                   all.x=T,all.y=F)
  
  # Rename the coverage statistics BY NAME, not by position. merge() places the
  # `by` columns first, so a positional `colnames(summdf2) <- colns` silently
  # mislabels every column the moment the join key changes -- which is exactly
  # what happened when the key became (ont_multiplex_group, ont_barcode):
  # gene_target was relabelled "sample_id", so nothing ever matched a sample and
  # every marker came out with n_callable = 0.
  summdf2 <- summdf2 %>%
    rename(coverage_min    = `Min.`,
           coverage_Q1     = `1st Qu.`,
           coverage_median = `Median`,
           coverage_mean   = `Mean`,
           coverage_Q3     = `3rd Qu.`,
           coverage_max    = `Max.`)
  
  col_ord <- c("sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date", "ont_barcode",
               "gene_target",
               "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean", "coverage_Q3", "coverage_max")
  
  summdf2 <- summdf2 %>%
    select(all_of(col_ord))
  
  list_crt_cov[[j]] <- summdf2
  
}

cov_df_crt <- data.frame(do.call("rbind", list_crt_cov))

# summary stats
summary_cov_crt <- summary(cov_df_crt$coverage_median)

median_median_cov_crt <- median(cov_df_crt$coverage_median)
median_Q1_cov_crt <- round(as.numeric(summary(cov_df_crt$coverage_median)[2]),0)
median_Q3_cov_crt <- round(as.numeric(summary(cov_df_crt$coverage_median)[5]),0)
median_IQR_cov_crt <- median_Q3_cov_crt - median_Q1_cov_crt

lower_cov_5pcnt_crt <- round(as.numeric(quantile(cov_df_crt$coverage_median, probs = 0.05, na.rm = TRUE)), 0)
#lower_cov_10pcnt_mdr1 <- round(as.numeric(quantile(cov_df_mdr1$coverage_median, probs=.1)),0)
lower_cov_10pcnt_crt <- round(as.numeric(quantile(cov_df_crt$coverage_median, probs = 0.10, na.rm = TRUE)), 0)



########## dhfr
# read in bed file
filenames <- list.files(bedfile_dir, pattern="_dhfr.bedGraph", full.names=TRUE)
if (length(filenames) == 0) {
  stop("No DHFR bedGraph files found in: ", bedfile_dir)
}
ldf <- lapply(filenames, read.table)

#barcode = str_sub(filenames, 142, 150)
barcode <- str_extract(filenames, "barcode\\d+")
run_of_file <- nr_run_from_file(filenames)   # barcode01 exists on every flowcell
max_count <- length(ldf)

### Filter variants to only include those within the amplicon target region
list_dhfr <- list()
for (i in 1:max_count) {
  
  ldf[[i]] <- ldf[[i]] %>%
    rename(gene_id = V1,
           start = V2,
           end = V3,
           depth = V4)
  
  ldf[[i]]$start <- as.numeric(ldf[[i]]$start)
  ldf[[i]]$end <- as.numeric(ldf[[i]]$end)
  ldf[[i]]$depth <- as.numeric(ldf[[i]]$depth)
  
  ldf[[i]] <- ldf[[i]] %>%
    filter(start > amp_pos_dhfr$amplicon_start) %>%
    filter(end < amp_pos_dhfr$amplicon_end)
  
  list_dhfr[[i]] <- ldf[[i]]
  
}

### Now extract coverage statistics for each amplicon
list_dhfr_cov <- list()

for (j in 1:max_count) {
  
  summ <- summary(list_dhfr[[j]]$depth)
  summdf <- data.frame(summ=matrix(summ, ncol=6))
  colnames(summdf) <- names(summ)
  

  ## Reproduce barcode ID
  id <- barcode[j]
  
  summdf$ont_barcode <- id
  summdf$ont_multiplex_group <- run_of_file[j]
  summdf$gene_target <- "dhfr"
  
  summdf2 <- merge(summdf, samp_names, by=c("ont_multiplex_group","ont_barcode"),
                   all.x=T,all.y=F)
  
  # Rename the coverage statistics BY NAME, not by position. merge() places the
  # `by` columns first, so a positional `colnames(summdf2) <- colns` silently
  # mislabels every column the moment the join key changes -- which is exactly
  # what happened when the key became (ont_multiplex_group, ont_barcode):
  # gene_target was relabelled "sample_id", so nothing ever matched a sample and
  # every marker came out with n_callable = 0.
  summdf2 <- summdf2 %>%
    rename(coverage_min    = `Min.`,
           coverage_Q1     = `1st Qu.`,
           coverage_median = `Median`,
           coverage_mean   = `Mean`,
           coverage_Q3     = `3rd Qu.`,
           coverage_max    = `Max.`)
  
  col_ord <- c("sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date", "ont_barcode",
               "gene_target",
               "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean", "coverage_Q3", "coverage_max")
  
  summdf2 <- summdf2 %>%
    select(all_of(col_ord))
  
  list_dhfr_cov[[j]] <- summdf2
  
}

cov_df_dhfr <- data.frame(do.call("rbind", list_dhfr_cov))

# summary stats
summary_cov_dhfr <- summary(cov_df_dhfr$coverage_median)

median_median_cov_dhfr <- median(cov_df_dhfr$coverage_median)
median_Q1_cov_dhfr <- round(as.numeric(summary(cov_df_dhfr$coverage_median)[2]),0)
median_Q3_cov_dhfr <- round(as.numeric(summary(cov_df_dhfr$coverage_median)[5]),0)
median_IQR_cov_dhfr <- median_Q3_cov_dhfr - median_Q1_cov_dhfr
lower_cov_5pcnt_dhfr <- round(as.numeric(quantile(cov_df_dhfr$coverage_median, probs=.05 , na.rm = TRUE)),0)
lower_cov_10pcnt_dhfr <- round(as.numeric(quantile(cov_df_dhfr$coverage_median, probs=.1, na.rm = TRUE)),0)



########## dhps
# read in bed file
filenames <- list.files(bedfile_dir, pattern="_dhps.bedGraph", full.names=TRUE)
if (length(filenames) == 0) {
  stop("No DHPS bedGraph files found in: ", bedfile_dir)
}
ldf <- lapply(filenames, read.table)

#barcode = str_sub(filenames, 142, 150)
barcode <- str_extract(filenames, "barcode\\d+")
run_of_file <- nr_run_from_file(filenames)   # barcode01 exists on every flowcell
max_count <- length(ldf)

### Filter variants to only include those within the amplicon target region
list_dhps <- list()
for (i in 1:max_count) {
  
  ldf[[i]] <- ldf[[i]] %>%
    rename(gene_id = V1,
           start = V2,
           end = V3,
           depth = V4)
  
  ldf[[i]]$start <- as.numeric(ldf[[i]]$start)
  ldf[[i]]$end <- as.numeric(ldf[[i]]$end)
  ldf[[i]]$depth <- as.numeric(ldf[[i]]$depth)
  
  ldf[[i]] <- ldf[[i]] %>%
    filter(start > amp_pos_dhps$amplicon_start) %>%
    filter(end < amp_pos_dhps$amplicon_end)
  
  list_dhps[[i]] <- ldf[[i]]
  
}

### Now extract coverage statistics for each amplicon
list_dhps_cov <- list()

for (j in 1:max_count) {
  
  summ <- summary(list_dhps[[j]]$depth)
  summdf <- data.frame(summ=matrix(summ, ncol=6))
  colnames(summdf) <- names(summ)
  
  ## Reproduce barcode ID
  id <- barcode[j]
  
  summdf$ont_barcode <- id
  summdf$ont_multiplex_group <- run_of_file[j]
  summdf$gene_target <- "dhps"
  
  summdf2 <- merge(summdf, samp_names, by=c("ont_multiplex_group","ont_barcode"),
                   all.x=T,all.y=F)
  
  # Rename the coverage statistics BY NAME, not by position. merge() places the
  # `by` columns first, so a positional `colnames(summdf2) <- colns` silently
  # mislabels every column the moment the join key changes -- which is exactly
  # what happened when the key became (ont_multiplex_group, ont_barcode):
  # gene_target was relabelled "sample_id", so nothing ever matched a sample and
  # every marker came out with n_callable = 0.
  summdf2 <- summdf2 %>%
    rename(coverage_min    = `Min.`,
           coverage_Q1     = `1st Qu.`,
           coverage_median = `Median`,
           coverage_mean   = `Mean`,
           coverage_Q3     = `3rd Qu.`,
           coverage_max    = `Max.`)
  
  col_ord <- c("sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date", "ont_barcode",
               "gene_target",
               "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean", "coverage_Q3", "coverage_max")
  
  summdf2 <- summdf2 %>%
    select(all_of(col_ord))
  
  list_dhps_cov[[j]] <- summdf2
  
}

cov_df_dhps <- data.frame(do.call("rbind", list_dhps_cov))

# summary stats
summary_cov_dhps <- summary(cov_df_dhps$coverage_median)

median_median_cov_dhps <- median(cov_df_dhps$coverage_median)
median_Q1_cov_dhps <- round(as.numeric(summary(cov_df_dhps$coverage_median)[2]),0)
median_Q3_cov_dhps <- round(as.numeric(summary(cov_df_dhps$coverage_median)[5]),0)
median_IQR_cov_dhps <- median_Q3_cov_dhps - median_Q1_cov_dhps
lower_cov_5pcnt_dhps <- round(as.numeric(quantile(cov_df_dhps$coverage_median, probs=.05, na.rm = TRUE)),0)
lower_cov_10pcnt_dhps <- round(as.numeric(quantile(cov_df_dhps$coverage_median, probs=.1, na.rm = TRUE)),0)



########## mdr1
# read in bed file
filenames <- list.files(bedfile_dir, pattern="_mdr1.bedGraph", full.names=TRUE)
if (length(filenames) == 0) {
  stop("No MDR1 bedGraph files found in: ", bedfile_dir)
}
ldf <- lapply(filenames, read.table)

#barcode = str_sub(filenames, 142, 150)
barcode <- str_extract(filenames, "barcode\\d+")
run_of_file <- nr_run_from_file(filenames)   # barcode01 exists on every flowcell
max_count <- length(ldf)

### Filter variants to only include those within the amplicon target region
list_mdr1 <- list()
for (i in 1:max_count) {
  
  ldf[[i]] <- ldf[[i]] %>%
    rename(gene_id = V1,
           start = V2,
           end = V3,
           depth = V4)
  
  ldf[[i]]$start <- as.numeric(ldf[[i]]$start)
  ldf[[i]]$end <- as.numeric(ldf[[i]]$end)
  ldf[[i]]$depth <- as.numeric(ldf[[i]]$depth)
  
  ldf[[i]] <- ldf[[i]] %>%
    filter(start > amp_pos_mdr1$amplicon_start) %>%
    filter(end < amp_pos_mdr1$amplicon_end)
  
  list_mdr1[[i]] <- ldf[[i]]
  
}

### Now extract coverage statistics for each amplicon
list_mdr1_cov <- list()

for (j in 1:max_count) {
  
  summ <- summary(list_mdr1[[j]]$depth)
  summdf <- data.frame(summ=matrix(summ, ncol=6))
  colnames(summdf) <- names(summ)
  
  ## Reproduce barcode ID
  id <- barcode[j]
  
  summdf$ont_barcode <- id
  summdf$ont_multiplex_group <- run_of_file[j]
  summdf$gene_target <- "mdr1"
  
  summdf2 <- merge(summdf, samp_names, by=c("ont_multiplex_group","ont_barcode"),
                   all.x=T,all.y=F)
  
  # Rename the coverage statistics BY NAME, not by position. merge() places the
  # `by` columns first, so a positional `colnames(summdf2) <- colns` silently
  # mislabels every column the moment the join key changes -- which is exactly
  # what happened when the key became (ont_multiplex_group, ont_barcode):
  # gene_target was relabelled "sample_id", so nothing ever matched a sample and
  # every marker came out with n_callable = 0.
  summdf2 <- summdf2 %>%
    rename(coverage_min    = `Min.`,
           coverage_Q1     = `1st Qu.`,
           coverage_median = `Median`,
           coverage_mean   = `Mean`,
           coverage_Q3     = `3rd Qu.`,
           coverage_max    = `Max.`)
  
  col_ord <- c("sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date", "ont_barcode",
               "gene_target",
               "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean", "coverage_Q3", "coverage_max")
  
  summdf2 <- summdf2 %>%
    select(all_of(col_ord))
  
  list_mdr1_cov[[j]] <- summdf2
  
}

cov_df_mdr1 <- data.frame(do.call("rbind", list_mdr1_cov))

# summary stats
summary_cov_mdr1 <- summary(cov_df_mdr1$coverage_median)

median_median_cov_mdr1 <- median(cov_df_mdr1$coverage_median)
median_Q1_cov_mdr1 <- round(as.numeric(summary(cov_df_mdr1$coverage_median)[2]),0)
median_Q3_cov_mdr1 <- round(as.numeric(summary(cov_df_mdr1$coverage_median)[5]),0)
median_IQR_cov_mdr1 <- median_Q3_cov_mdr1 - median_Q1_cov_mdr1
#lower_cov_5pcnt_mdr1 <- round(as.numeric(quantile(cov_df_mdr1$coverage_median, probs=.05)),0)
#lower_cov_5pcnt_mdr1 <- round(as.numeric(quantile(cov_df_mdr1$coverage_median, probs = 0.05, na.rm = TRUE)), 0)
lower_cov_5pcnt_mdr1 <- round(as.numeric(quantile(cov_df_mdr1$coverage_median, probs = 0.05, na.rm = TRUE)), 0)

#lower_cov_10pcnt_mdr1 <- round(as.numeric(quantile(cov_df_mdr1$coverage_median, probs=.1)),0)
lower_cov_10pcnt_mdr1 <- round(as.numeric(quantile(cov_df_mdr1$coverage_median, probs = 0.10, na.rm = TRUE)), 0)



########## kelch13
# read in bed file
filenames <- list.files(bedfile_dir, pattern="_k13.bedGraph", full.names=TRUE)
if (length(filenames) == 0) {
  stop("No K13 bedGraph files found in: ", bedfile_dir)
}
ldf <- lapply(filenames, read.table)

#barcode = str_sub(filenames, 142, 150)
barcode <- str_extract(filenames, "barcode\\d+")
run_of_file <- nr_run_from_file(filenames)   # barcode01 exists on every flowcell
max_count <- length(ldf)

### Filter variants to only include those within the amplicon target region
list_k13 <- list()
for (i in 1:max_count) {
  
  ldf[[i]] <- ldf[[i]] %>%
    rename(gene_id = V1,
           start = V2,
           end = V3,
           depth = V4)
  
  ldf[[i]]$start <- as.numeric(ldf[[i]]$start)
  ldf[[i]]$end <- as.numeric(ldf[[i]]$end)
  ldf[[i]]$depth <- as.numeric(ldf[[i]]$depth)
  
  ldf[[i]] <- ldf[[i]] %>%
    filter(start > amp_pos_kelch13$amplicon_start) %>%
    filter(end < amp_pos_kelch13$amplicon_end)
  
  list_k13[[i]] <- ldf[[i]]
  
}

### Now extract coverage statistics for each amplicon
list_k13_cov <- list()

for (j in 1:max_count) {
  
  summ <- summary(list_k13[[j]]$depth)
  summdf <- data.frame(summ=matrix(summ, ncol=6))
  colnames(summdf) <- names(summ)
  
  ## Reproduce barcode ID
  id <- barcode[j]
  
  summdf$ont_barcode <- id
  summdf$ont_multiplex_group <- run_of_file[j]
  summdf$gene_target <- "kelch13"
  
  summdf2 <- merge(summdf, samp_names, by=c("ont_multiplex_group","ont_barcode"),
                   all.x=T,all.y=F)
  
  # Rename the coverage statistics BY NAME, not by position. merge() places the
  # `by` columns first, so a positional `colnames(summdf2) <- colns` silently
  # mislabels every column the moment the join key changes -- which is exactly
  # what happened when the key became (ont_multiplex_group, ont_barcode):
  # gene_target was relabelled "sample_id", so nothing ever matched a sample and
  # every marker came out with n_callable = 0.
  summdf2 <- summdf2 %>%
    rename(coverage_min    = `Min.`,
           coverage_Q1     = `1st Qu.`,
           coverage_median = `Median`,
           coverage_mean   = `Mean`,
           coverage_Q3     = `3rd Qu.`,
           coverage_max    = `Max.`)
  
  col_ord <- c("sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date", "ont_barcode",
               "gene_target",
               "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean", "coverage_Q3", "coverage_max")
  
  summdf2 <- summdf2 %>%
    select(all_of(col_ord))
  
  list_k13_cov[[j]] <- summdf2
  
}

cov_df_k13 <- data.frame(do.call("rbind", list_k13_cov))

# summary stats
summary_cov_k13 <- summary(cov_df_k13$coverage_median)

median_median_cov_k13 <- median(cov_df_k13$coverage_median)
median_Q1_cov_k13 <- round(as.numeric(summary(cov_df_k13$coverage_median)[2]),0)
median_Q3_cov_k13 <- round(as.numeric(summary(cov_df_k13$coverage_median)[5]),0)
median_IQR_cov_k13 <- median_Q3_cov_k13 - median_Q1_cov_k13
lower_cov_5pcnt_k13 <- round(as.numeric(quantile(cov_df_k13$coverage_median, probs=.05, na.rm = TRUE)),0)
lower_cov_10pcnt_k13 <- round(as.numeric(quantile(cov_df_k13$coverage_median, probs=.1, na.rm = TRUE)),0)



########## csp
# read in bed file
filenames <- list.files(bedfile_dir, pattern="_csp.bedGraph", full.names=TRUE)
if (length(filenames) == 0) {
  stop("No CSP bedGraph files found in: ", bedfile_dir)
}
ldf <- lapply(filenames, read.table)

#barcode = str_sub(filenames, 142, 150)
barcode <- str_extract(filenames, "barcode\\d+")
run_of_file <- nr_run_from_file(filenames)   # barcode01 exists on every flowcell
max_count <- length(ldf)

### Filter variants to only include those within the amplicon target region
list_csp <- list()
for (i in 1:max_count) {
  
  ldf[[i]] <- ldf[[i]] %>%
    rename(gene_id = V1,
           start = V2,
           end = V3,
           depth = V4)
  
  ldf[[i]]$start <- as.numeric(ldf[[i]]$start)
  ldf[[i]]$end <- as.numeric(ldf[[i]]$end)
  ldf[[i]]$depth <- as.numeric(ldf[[i]]$depth)
  
  ldf[[i]] <- ldf[[i]] %>%
    filter(start > amp_pos_csp$amplicon_start) %>%
    filter(end < amp_pos_csp$amplicon_end)
  
  list_csp[[i]] <- ldf[[i]]
  
}

### Now extract coverage statistics for each amplicon
list_csp_cov <- list()

for (j in 1:max_count) {
  
  summ <- summary(list_csp[[j]]$depth)
  summdf <- data.frame(summ=matrix(summ, ncol=6))
  colnames(summdf) <- names(summ)
  
  ## Reproduce barcode ID
  id <- barcode[j]
  
  summdf$ont_barcode <- id
  summdf$ont_multiplex_group <- run_of_file[j]
  summdf$gene_target <- "csp"
  
  summdf2 <- merge(summdf, samp_names, by=c("ont_multiplex_group","ont_barcode"),
                   all.x=T,all.y=F)
  
  # Rename the coverage statistics BY NAME, not by position. merge() places the
  # `by` columns first, so a positional `colnames(summdf2) <- colns` silently
  # mislabels every column the moment the join key changes -- which is exactly
  # what happened when the key became (ont_multiplex_group, ont_barcode):
  # gene_target was relabelled "sample_id", so nothing ever matched a sample and
  # every marker came out with n_callable = 0.
  summdf2 <- summdf2 %>%
    rename(coverage_min    = `Min.`,
           coverage_Q1     = `1st Qu.`,
           coverage_median = `Median`,
           coverage_mean   = `Mean`,
           coverage_Q3     = `3rd Qu.`,
           coverage_max    = `Max.`)
  
  col_ord <- c("sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date", "ont_barcode",
               "gene_target",
               "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean", "coverage_Q3", "coverage_max")
  
  summdf2 <- summdf2 %>%
    select(all_of(col_ord))
  
  list_csp_cov[[j]] <- summdf2
  
}

cov_df_csp <- data.frame(do.call("rbind", list_csp_cov))

# summary stats
summary_cov_csp <- summary(cov_df_csp$coverage_median)

median_median_cov_csp <- median(cov_df_csp$coverage_median)
median_Q1_cov_csp <- round(as.numeric(summary(cov_df_csp$coverage_median)[2]),0)
median_Q3_cov_csp <- round(as.numeric(summary(cov_df_csp$coverage_median)[5]),0)
median_IQR_cov_csp <- median_Q3_cov_csp - median_Q1_cov_csp
lower_cov_5pcnt_csp <- round(as.numeric(quantile(cov_df_csp$coverage_median, probs=.05, na.rm = TRUE)),0)
lower_cov_10pcnt_csp <- round(as.numeric(quantile(cov_df_csp$coverage_median, probs=.1, na.rm = TRUE)),0)



########## msp1
# read in bed file
#filenames <- list.files(bedfile_dir, pattern="_msp1.bedGraph", full.names=TRUE)
#ldf <- lapply(filenames, read.table)
#
#max_count <- length(ldf)
#
### Filter variants to only include those within the amplicon target region
#list_msp1 <- list()
#for (i in 1:max_count) {
#  
#  ldf[[i]] <- ldf[[i]] %>%
#    rename(gene_id = V1,
#           start = V2,
#           end = V3,
#           depth = V4)
#  
#  ldf[[i]]$start <- as.numeric(ldf[[i]]$start)
#  ldf[[i]]$end <- as.numeric(ldf[[i]]$end)
#  ldf[[i]]$depth <- as.numeric(ldf[[i]]$depth)
#  
#  ldf[[i]] <- ldf[[i]] %>%
#    filter(start > amp_pos_msp1$amplicon_start) %>%
#    filter(end < amp_pos_msp1$amplicon_end)
#  
#  list_msp1[[i]] <- ldf[[i]]
#  
#}
#
### Now extract coverage statistics for each amplicon
#list_msp1_cov <- list()
#
#for (j in 1:max_count) {
#  
#  summ <- summary(list_msp1[[j]]$depth)
#  summdf <- data.frame(summ=matrix(summ, ncol=6))
#  colnames(summdf) <- names(summ)
#  
#  ## Reproduce barcode ID
#  id <- ifelse(nchar(j)==1, paste0('barcode0', as.character(j)),
#               ifelse(nchar(j)==2, paste0('barcode', as.character(j)), NA)
#  )
#  
#  summdf$ont_barcode <- id
  summdf$ont_multiplex_group <- run_of_file[j]
#  summdf$gene_target <- "msp1"
#  
#  summdf2 <- merge(summdf, samp_names, by=c("ont_multiplex_group","ont_barcode"),
#                   all.x=T,all.y=F)
#  
#  colns <- c("ont_barcode", "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean",
#             "coverage_Q3", "coverage_max", "gene_target", "sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date")
#  
#  colnames(summdf2) <- colns
#  
#  col_ord <- c("sample_id", "patient_id", "ont_multiplex_group", "ont_seq_date", "ont_barcode",
#               "gene_target",
#               "coverage_min", "coverage_Q1", "coverage_median", "coverage_mean", "coverage_Q3", "coverage_max")
#  
#  summdf2 <- summdf2 %>%
#    select(all_of(col_ord))
#  
#  list_msp1_cov[[j]] <- summdf2
#  
#}
#
#cov_df_msp1 <- data.frame(do.call("rbind", list_msp1_cov))
#
## summary stats
#summary_cov_msp1 <- summary(cov_df_msp1$coverage_median)
#
#median_median_cov_msp1 <- median(cov_df_msp1$coverage_median)
#median_Q1_cov_msp1 <- round(as.numeric(summary(cov_df_msp1$coverage_median)[2]),0)
#median_Q3_cov_msp1 <- round(as.numeric(summary(cov_df_msp1$coverage_median)[5]),0)
#median_IQR_cov_msp1 <- median_Q3_cov_msp1 - median_Q1_cov_msp1
#lower_cov_5pcnt_msp1 <- round(as.numeric(quantile(cov_df_msp1$coverage_median, probs=.05)),0)
#lower_cov_10pcnt_msp1 <- round(as.numeric(quantile(cov_df_msp1$coverage_median, probs=.1)),0)
#


### Summary statistics for the whole run
gene_target <- c("crt", "dhfr", "dhps", "mdr1", "kelch13", "csp") #msp1

coverage_median <- c(median_median_cov_crt, median_median_cov_dhfr, median_median_cov_dhps, median_median_cov_mdr1, median_median_cov_k13, median_median_cov_csp) #median_median_cov_msp1
coverage_Q1 <- c(median_Q1_cov_crt, median_Q1_cov_dhfr, median_Q1_cov_dhps, median_Q1_cov_mdr1, median_Q1_cov_k13, median_Q1_cov_csp) #median_Q1_cov_msp1
coverage_Q3 <- c(median_Q3_cov_crt, median_Q3_cov_dhfr, median_Q3_cov_dhps, median_Q3_cov_mdr1, median_Q3_cov_k13, median_Q3_cov_csp) #median_Q3_cov_msp1
coverage_IQR <- c(median_IQR_cov_crt, median_IQR_cov_dhfr, median_IQR_cov_dhps, median_IQR_cov_mdr1, median_IQR_cov_k13, median_IQR_cov_csp) #median_IQR_cov_msp1
coverage_lower_5pcnt <- c(lower_cov_5pcnt_crt, lower_cov_5pcnt_dhfr, lower_cov_5pcnt_dhps, lower_cov_5pcnt_mdr1, lower_cov_5pcnt_k13, lower_cov_5pcnt_csp) #lower_cov_5pcnt_msp1
coverage_lower_10pcnt <- c(lower_cov_10pcnt_crt, lower_cov_10pcnt_dhfr, lower_cov_10pcnt_dhps, lower_cov_10pcnt_mdr1, lower_cov_10pcnt_k13, lower_cov_10pcnt_csp) #lower_cov_10pcnt_msp1

median_cov_stats <- data.frame(gene_target, coverage_median, coverage_Q1, coverage_Q3, coverage_IQR, coverage_lower_5pcnt, coverage_lower_10pcnt)

median_cov_stats



#### Identify low coverage samples - RELATIVE lower 5% or 10%
cov_df_crt$cov_below_5pcnt <-
  ifelse(cov_df_crt$coverage_median < lower_cov_5pcnt_crt,
         TRUE, FALSE)
cov_df_crt$cov_below_10pcnt <-
  ifelse(cov_df_crt$coverage_median < lower_cov_10pcnt_crt,
         TRUE, FALSE)

cov_df_dhfr$cov_below_5pcnt <-
  ifelse(cov_df_dhfr$coverage_median < lower_cov_5pcnt_dhfr,
         TRUE, FALSE)
cov_df_dhfr$cov_below_10pcnt <-
  ifelse(cov_df_dhfr$coverage_median < lower_cov_10pcnt_dhfr,
         TRUE, FALSE)

cov_df_dhps$cov_below_5pcnt <-
  ifelse(cov_df_dhps$coverage_median < lower_cov_5pcnt_dhps,
         TRUE, FALSE)
cov_df_dhps$cov_below_10pcnt <-
  ifelse(cov_df_dhps$coverage_median < lower_cov_10pcnt_dhps,
         TRUE, FALSE)

cov_df_mdr1$cov_below_5pcnt <-
  ifelse(cov_df_mdr1$coverage_median < lower_cov_5pcnt_mdr1,
         TRUE, FALSE)
cov_df_mdr1$cov_below_10pcnt <-
  ifelse(cov_df_mdr1$coverage_median < lower_cov_10pcnt_mdr1,
         TRUE, FALSE)

cov_df_k13$cov_below_5pcnt <-
  ifelse(cov_df_k13$coverage_median < lower_cov_5pcnt_k13,
         TRUE, FALSE)
cov_df_k13$cov_below_10pcnt <-
  ifelse(cov_df_k13$coverage_median < lower_cov_10pcnt_k13,
         TRUE, FALSE)

cov_df_csp$cov_below_5pcnt <-
  ifelse(cov_df_csp$coverage_median < lower_cov_5pcnt_csp,
         TRUE, FALSE)
cov_df_csp$cov_below_10pcnt <-
  ifelse(cov_df_csp$coverage_median < lower_cov_10pcnt_csp,
         TRUE, FALSE)

#cov_df_msp1$cov_below_5pcnt <-
#  ifelse(cov_df_msp1$coverage_median < lower_cov_5pcnt_msp1,
#         TRUE, FALSE)
#cov_df_msp1$cov_below_10pcnt <-
#  ifelse(cov_df_msp1$coverage_median < lower_cov_10pcnt_msp1,
#         TRUE, FALSE)



#### Identify low coverage samples - ABSOLUTE coverage cutoff
#### (min_cov_threshold is set above; canonical run uses >=50x, sensitivity configs may differ.)
	cov_df_crt$coverage_above_threshold <-
	  ifelse(cov_df_crt$coverage_median >= min_cov_threshold,
	         "TRUE", "FALSE")
	
	cov_df_dhfr$coverage_above_threshold <-
	  ifelse(cov_df_dhfr$coverage_median >= min_cov_threshold,
	         "TRUE", "FALSE")
	
	cov_df_dhps$coverage_above_threshold <-
	  ifelse(cov_df_dhps$coverage_median >= min_cov_threshold,
	         "TRUE", "FALSE")
	
	cov_df_mdr1$coverage_above_threshold <-
	  ifelse(cov_df_mdr1$coverage_median >= min_cov_threshold,
	         "TRUE", "FALSE")
	
	cov_df_k13$coverage_above_threshold <-
	  ifelse(cov_df_k13$coverage_median >= min_cov_threshold,
	         "TRUE", "FALSE")
	
	cov_df_csp$coverage_above_threshold <-
	  ifelse(cov_df_csp$coverage_median >= min_cov_threshold,
	         "TRUE", "FALSE")

#cov_df_msp1$coverage_above_threshold <-
#  ifelse(cov_df_msp1$coverage_median > min_cov_threshold,
#         "TRUE", "FALSE")



#### Merge all the data
cov_df_sub_crt <- cov_df_crt %>%
  select(sample_id, coverage_median, cov_below_5pcnt, cov_below_10pcnt, coverage_above_threshold) %>%
  rename(crt_coverage_median = coverage_median,
         crt_coverage_below_5pcnt = cov_below_5pcnt,
         crt_coverage_below_10pcnt = cov_below_10pcnt,
         crt_coverage_above_threshold = coverage_above_threshold)

cov_df_sub_dhfr <- cov_df_dhfr %>%
  select(sample_id, coverage_median, cov_below_5pcnt, cov_below_10pcnt, coverage_above_threshold) %>%
  rename(dhfr_coverage_median = coverage_median,
         dhfr_coverage_below_5pcnt = cov_below_5pcnt,
         dhfr_coverage_below_10pcnt = cov_below_10pcnt,
         dhfr_coverage_above_threshold = coverage_above_threshold)

cov_df_sub_dhps <- cov_df_dhps %>%
  select(sample_id, coverage_median, cov_below_5pcnt, cov_below_10pcnt, coverage_above_threshold) %>%
  rename(dhps_coverage_median = coverage_median,
         dhps_coverage_below_5pcnt = cov_below_5pcnt,
         dhps_coverage_below_10pcnt = cov_below_10pcnt,
         dhps_coverage_above_threshold = coverage_above_threshold)

cov_df_sub_mdr1 <- cov_df_mdr1 %>%
  select(sample_id, coverage_median, cov_below_5pcnt, cov_below_10pcnt, coverage_above_threshold) %>%
  rename(mdr1_coverage_median = coverage_median,
         mdr1_coverage_below_5pcnt = cov_below_5pcnt,
         mdr1_coverage_below_10pcnt = cov_below_10pcnt,
         mdr1_coverage_above_threshold = coverage_above_threshold)

cov_df_sub_k13 <- cov_df_k13 %>%
  select(sample_id, coverage_median, cov_below_5pcnt, cov_below_10pcnt, coverage_above_threshold) %>%
  rename(k13_coverage_median = coverage_median,
         k13_coverage_below_5pcnt = cov_below_5pcnt,
         k13_coverage_below_10pcnt = cov_below_10pcnt,
         k13_coverage_above_threshold = coverage_above_threshold)

cov_df_sub_csp <- cov_df_csp %>%
  select(sample_id, coverage_median, cov_below_5pcnt, cov_below_10pcnt, coverage_above_threshold) %>%
  rename(csp_coverage_median = coverage_median,
         csp_coverage_below_5pcnt = cov_below_5pcnt,
         csp_coverage_below_10pcnt = cov_below_10pcnt,
         csp_coverage_above_threshold = coverage_above_threshold)

#cov_df_sub_msp1 <- cov_df_msp1 %>%
#  select(sample_id, coverage_median, cov_below_5pcnt, cov_below_10pcnt, coverage_above_threshold) %>%
#  rename(msp1_coverage_median = coverage_median,
#         msp1_coverage_below_5pcnt = cov_below_5pcnt,
#         msp1_coverage_below_10pcnt = cov_below_10pcnt,
#         msp1_coverage_above_threshold = coverage_above_threshold)

# merge (inner join to keep only samples present in both CRT and DHFR coverage)
cov_merged_df1 <- merge(cov_df_sub_crt, cov_df_sub_dhfr, by = "sample_id")



cov_merged_df2 <- merge(cov_merged_df1,cov_df_sub_dhps,
                        by="sample_id",
                        all.x=T, all.y=T)

cov_merged_df3 <- merge(cov_merged_df2,cov_df_sub_mdr1,
                        by="sample_id",
                        all.x=T, all.y=T)

cov_merged_df4 <- merge(cov_merged_df3,cov_df_sub_k13,
                        by="sample_id",
                        all.x=T, all.y=T)

cov_merged_df_final <- merge(cov_merged_df4,cov_df_sub_csp,
                             by="sample_id",
                             all.x=T, all.y=T)

#cov_merged_df_final <- merge(cov_merged_df5,cov_df_sub_msp1,
#                        by="sample_id",
#                        all.x=T, all.y=T)

head(cov_merged_df_final)
median_cov_stats


##### Sanity check before writing
#
# This stage joins bedGraph coverage to the sample sheet. If that join silently
# fails, every marker downstream comes out with n_callable = 0 and the pipeline
# still exits 0 -- which is exactly what happened when a positional
# `colnames() <- ...` mislabelled the columns after the join key changed. Fail
# loudly instead of shipping a plausible-looking but empty result.
if (!"sample_id" %in% names(cov_merged_df_final)) {
  stop("[stage1] coverage table has no sample_id column - the join to the sample sheet broke.")
}

matched <- sum(cov_merged_df_final$sample_id %in% meta$sample_id, na.rm = TRUE)
if (matched == 0) {
  stop("[stage1] not one coverage row matched a sample_id in the sample sheet. ",
       "The bedGraph -> sample-sheet join is broken (join key, or column ordering). ",
       "Saw sample_id values like: ",
       paste(utils::head(unique(cov_merged_df_final$sample_id), 5), collapse = ", "))
}

# The final table merges six per-gene tables by sample_id, so it must have at most
# one row per sample. If sample_id is not unique (controls are re-sequenced on
# every flowcell), each merge multiplies the duplicates and the table explodes --
# 5 flowcells x 6 merges produced 15,625 rows for a single control.
if (nrow(cov_merged_df_final) > nrow(meta)) {
  offenders <- names(which(table(cov_merged_df_final$sample_id) > 1))
  stop(sprintf(paste0("[stage1] coverage table has %d rows but the cohort has only %d samples - ",
                      "the per-gene merges have fanned out on a non-unique sample_id. Offenders: %s"),
               nrow(cov_merged_df_final), nrow(meta),
               paste(utils::head(offenders, 5), collapse = ", ")))
}

message(sprintf("[stage1] coverage rows: %d (cohort samples: %d); matched to a sample_id: %d (%.0f%%)",
                nrow(cov_merged_df_final), nrow(meta), matched,
                100 * matched / nrow(cov_merged_df_final)))

##### Save outputs
write.csv(median_cov_stats, file=overall_cov_fn, row.names=FALSE)
write.csv(cov_merged_df_final, file=cov_by_sample_fn, row.names=FALSE)
