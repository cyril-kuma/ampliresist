#!/usr/bin/env Rscript
# =============================================================================
# 07_integrate_metadata.R
#
# Builds the analysis-ready specimen table by joining Objective 3 (this pipeline)
# to Objective 2 (qPCR: parasite density as Cq) and Objective 1 (bloodmeal host
# metabarcoding). All three key on specimen_id.
#
# It also runs the CALLABILITY-BIAS TEST, which is the reason this join is not
# optional.
#
# THE PROBLEM IT EXISTS TO ANSWER
# -------------------------------
# dhfr is callable in only ~28% of specimens at 50x. Those are not a random 28%:
# a specimen is callable because it had enough parasite DNA, so the callable set
# is enriched for high parasitaemia BY CONSTRUCTION. If resistance genotype
# correlates with parasite density -- and there is every reason to think it might
# -- then the dhfr frequencies are computed on a biased subsample and are wrong
# in an unknown direction.
#
# Objective 2 gives us the instrument to test this: Cq is an inverse proxy for
# parasite density (low Cq = high density). So:
#
#   Compare Cq between dhfr-CALLABLE and dhfr-UNCALLABLE specimens.
#
#   If they differ  -> callability is density-dependent, the dhfr frequencies
#                      carry a selection bias, and that must be stated as a
#                      limitation (and ideally corrected/weighted).
#   If they do not  -> callability is missing-at-random with respect to density,
#                      and the frequencies stand.
#
# This is a real test with a real chance of failing. It is not a formality.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
samplesheet_xlsx <- args[1]   # Obj3 assets/samplesheet.xlsx
qpcr_csv         <- args[2]   # Obj2 results/processed/mosquito_level.csv
host_tsv         <- args[3]   # Obj1 05_endpoint_files/host_call_table.tsv
coverage_csv     <- args[4]   # Obj3 stage-1 coverage table (per-specimen callability)
coi_tsv          <- args[5]   # Obj3 06_complexity_of_infection.R output
out_tsv          <- args[6]

meta <- readxl::read_excel(samplesheet_xlsx)

# --- The join key, and why it is not specimen_id ----------------------------
# specimen_id is NOT unique across the study: six ids are reused by two
# genuinely different mosquitoes, collected at different sites in different
# years (e.g. MOS1071 is both FT027 -- Bekwai, 2024, Pf-negative -- and FT275 --
# Prang, 2025, Pf-positive). Joining on it silently fans out the table AND
# attaches the wrong mosquito's Cq to a genotype.
#
# The two objectives also SWAP their column names, which is its own trap:
#
#            bare code (FT275)      concatenation (FT275MOS1071)
#   Obj2:    sample_code            sample_id
#   Obj3:    sample_id              sample_code
#
# The concatenation embeds both ids, so it is unique and it disambiguates the
# collisions. That is the key we use. Obj2's processed table only carries the
# two parts, so rebuild it here.
qpcr <- readr::read_csv(qpcr_csv, show_col_types = FALSE) %>%
  mutate(join_key = paste0(sample_code, specimen_id)) %>%
  select(join_key, Pf_call, Pf_cq, qpcr_collection_date = collection_date)

if (any(duplicated(qpcr$join_key)))
  stop("[meta] qPCR join key is not unique - a join on it would fan out")

# --- Objective 1: bloodmeal host --------------------------------------------
# host_call_table is MARKER-level (several rows per specimen, one per marker and
# candidate host). Collapse to the dominant host: the assignment carrying the
# largest read fraction. Keep the fraction so a near-tie stays visible rather
# than being silently resolved.
host <- readr::read_tsv(host_tsv, show_col_types = FALSE) %>%
  filter(is.na(control_status) | control_status == "sample") %>%
  group_by(sample_id) %>%
  arrange(desc(host_fraction)) %>%
  summarise(host_primary       = dplyr::first(host_assignment),
            host_fraction      = dplyr::first(host_fraction),
            host_n_assignments = dplyr::n_distinct(host_assignment),
            .groups = "drop") %>%
  mutate(host_is_mixed = host_n_assignments > 1)

n_meta <- nrow(meta)
master <- meta %>%
  left_join(qpcr, by = c("sample_code" = "join_key")) %>%   # Obj3 sample_code IS the concatenation
  left_join(host, by = "sample_id")                          # Obj3 sample_id IS the bare code

# A left join must never ADD rows. If it does, a key was not unique and some
# specimen has just been silently duplicated -- which is how a cohort of 457
# becomes 463 and every downstream denominator is wrong.
if (nrow(master) != n_meta)
  stop("[meta] join fanned out: ", n_meta, " specimens in, ", nrow(master),
       " out. A join key is not unique.")

# --- Objective 3: callability + complexity ----------------------------------
cov <- readr::read_csv(coverage_csv, show_col_types = FALSE)
cov_cols <- grep("_coverage_above_threshold$", names(cov), value = TRUE)
if (length(cov_cols) > 0) {
  master <- master %>%
    left_join(cov %>% select(sample_id, all_of(cov_cols)), by = "sample_id")
}

if (file.exists(coi_tsv)) {
  coi <- readr::read_tsv(coi_tsv, show_col_types = FALSE) %>%
    select(ont_multiplex_group = run_name, ont_barcode = barcode,
           n_het_real, infection_class)
  master <- master %>% left_join(coi, by = c("ont_multiplex_group", "ont_barcode"))
}

readr::write_tsv(master, out_tsv)
message("[meta] wrote ", out_tsv, "  (", nrow(master), " specimens)")
message("[meta] with Cq: ",   sum(!is.na(master$Pf_cq)),
        " | with host: ",     sum(!is.na(master$host_primary)))

# =============================================================================
# CALLABILITY-BIAS TEST
# =============================================================================
cat("\n=====================================================================\n")
cat(" Callability-bias test: is a marker callable because the specimen had\n")
cat(" more parasite DNA?  (Cq is an INVERSE proxy: low Cq = high density)\n")
cat("=====================================================================\n\n")

for (cc in cov_cols) {
  gene <- sub("_coverage_above_threshold$", "", cc)
  d <- master %>% filter(!is.na(Pf_cq), !is.na(.data[[cc]]))
  if (nrow(d) < 10) next

  callable <- d$Pf_cq[d[[cc]] %in% c(TRUE, "TRUE", 1, "yes")]
  uncall   <- d$Pf_cq[!d[[cc]] %in% c(TRUE, "TRUE", 1, "yes")]
  if (length(callable) < 3 || length(uncall) < 3) {
    cat(sprintf("  %-8s  n_callable=%-4d n_uncallable=%-4d  (too few to test)\n",
                gene, length(callable), length(uncall)))
    next
  }

  w <- suppressWarnings(stats::wilcox.test(callable, uncall))
  flag <- if (w$p.value < 0.05) "  <-- BIASED" else ""
  cat(sprintf("  %-8s  callable n=%-4d Cq=%5.1f  |  uncallable n=%-4d Cq=%5.1f  |  dCq=%+5.1f  p=%.2e%s\n",
              gene, length(callable), median(callable, na.rm = TRUE),
              length(uncall),  median(uncall,  na.rm = TRUE),
              median(callable, na.rm = TRUE) - median(uncall, na.rm = TRUE),
              w$p.value, flag))
}

cat("\n  A significantly LOWER Cq in the callable group means callability tracks\n")
cat("  parasite density: that marker's frequencies are computed on a\n")
cat("  high-parasitaemia subsample and the selection bias must be reported.\n")
