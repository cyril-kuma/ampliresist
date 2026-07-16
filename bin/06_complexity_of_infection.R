#!/usr/bin/env Rscript
# =============================================================================
# 06_complexity_of_infection.R
#
# Flags polygenomic (multi-clone) infections, so that multilocus drug-resistance
# genotypes can be restricted to mono-infections.
#
# WHY NOT msp1
# ------------
# The canonical marker for complexity of infection is msp1 block 2 allelic-family
# typing. It is not available here: msp1 is in the reference manifest and is
# aligned against, but it carries ZERO reads in every specimen (mean max depth 0
# across the cohort; 0/40 bedGraphs with any coverage). The amplicon was never
# generated -- the primer pair is either absent from the multiplex or failed
# outright. No amount of reanalysis recovers it; it needs a wet-lab fix.
#
# WHAT WE DO INSTEAD
# ------------------
# P. falciparum is haploid in the host. A heterozygous call therefore does not
# mean a diploid genotype -- it means MORE THAN ONE CLONE is present. So
# within-sample heterozygosity across the surviving amplicons is a direct
# (if coarser) read on complexity of infection:
#
#     >= 1 heterozygous call at a non-artefact site  ->  polygenomic
#     0 heterozygous calls                           ->  consistent with mono-clonal
#
# "Consistent with" is doing real work in that second line: absence of a het call
# is weak evidence, since a second clone that is identical across our 6 amplicons
# is invisible. This flag is therefore CONSERVATIVE for polygenomic and NOT
# conclusive for monogenomic. It is a filter, not a COI estimate, and it must not
# be reported as one.
#
# EXCLUDING THE ARTEFACTS IS NOT OPTIONAL
# ---------------------------------------
# The four systematic-artefact positions (see docs/CALLABILITY.md) are called
# heterozygous in ~100% of specimens. Counting them would classify essentially
# the entire cohort as polygenomic -- the number would be an artefact of the
# assay, not a property of the parasites.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
per_call_files <- strsplit(args[1], ",")[[1]]   # one or more per_call.tsv
artefact_csv   <- args[2]                        # <cohort>_artefact_catalogue.csv
min_cov        <- as.numeric(args[3])            # depth threshold, e.g. 50
out_tsv        <- args[4]

# Read all per-call files with one stable schema. AF is deliberately character:
# biallelic records contain one value (e.g. 0.48), while multiallelic records
# contain one value per ALT allele (e.g. 0.46,0.27). Letting readr infer each
# file independently makes the former double and the latter character, which
# cannot be combined with bind_rows().
calls <- bind_rows(lapply(per_call_files, function(f)
  readr::read_tsv(
    f,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  ))) %>%
  mutate(
    pos = as.numeric(pos),
    DP  = as.numeric(DP)
  )

if (nrow(calls) == 0) stop("[coi] no calls read from per_call.tsv")

# snp_id in the artefact catalogue is <gene>_<pos>_<REF><ALT>; rebuild it here so
# the two tables key on the same thing.
calls <- calls %>%
  mutate(snp_id = paste0(gene, "_", pos, "_", ref, alt))

# Return whether GT contains at least two distinct called alleles. This covers
# 0/1 and 0|1 as before, plus multiallelic genotypes such as 1/2.
is_heterozygous <- function(gt) {
  if (is.na(gt) || gt == "") return(FALSE)
  alleles <- strsplit(gt, "[/|]")[[1]]
  alleles <- alleles[alleles != "" & alleles != "."]
  length(alleles) >= 2 && length(unique(alleles)) >= 2
}

# FORMAT/AF contains one frequency per ALT allele; the reference frequency is
# the remainder. For a heterozygous genotype, select the frequencies of its
# called alleles and return the smaller one. Thus 0/1 uses min(1-AF1, AF1),
# while 1/2 uses min(AF1, AF2) instead of discarding the multiallelic call.
minor_allele_fraction <- function(gt, af) {
  if (!is_heterozygous(gt) || is.na(af) || af == "") return(NA_real_)

  alt_af <- suppressWarnings(as.numeric(strsplit(af, ",", fixed = TRUE)[[1]]))
  if (length(alt_af) == 0 || anyNA(alt_af)) return(NA_real_)

  allele_af <- c(max(0, 1 - sum(alt_af)), alt_af)
  allele_ids <- suppressWarnings(as.integer(strsplit(gt, "[/|]")[[1]]))
  if (anyNA(allele_ids) || any(allele_ids + 1 > length(allele_af))) return(NA_real_)

  min(allele_af[allele_ids + 1])
}

calls <- calls %>%
  mutate(
    is_het = vapply(GT, is_heterozygous, logical(1)),
    minor_af = mapply(minor_allele_fraction, GT, AF)
  )

artefacts <- readr::read_csv(artefact_csv, show_col_types = FALSE) %>%
  filter(is_artefact) %>%
  pull(snp_id)

message("[coi] excluding ", length(artefacts), " artefact positions: ",
        paste(artefacts, collapse = ", "))

coi <- calls %>%
  filter(DP >= min_cov) %>%
  group_by(run_name, barcode) %>%
  summarise(
    n_sites_called   = dplyr::n(),
    n_het_all        = sum(is_het),
    n_het_real       = sum(is_het & !snp_id %in% artefacts),
    genes_with_het   = paste(sort(unique(gene[is_het &
                                              !snp_id %in% artefacts])), collapse = ";"),
    # At a het site the minor allele fraction indexes how balanced the clones are.
    # min(AF, 1-AF) so it is symmetric regardless of which allele is 'alt'.
    mean_minor_af    = ifelse(n_het_real > 0,
                              mean(minor_af[is_het & !snp_id %in% artefacts],
                                   na.rm = TRUE),
                              NA_real_),
    .groups = "drop") %>%
  mutate(
    infection_class = ifelse(n_het_real > 0, "polygenomic", "mono_compatible"),
    # Inflation caused by the artefacts, kept visible so the effect is auditable
    # rather than silently corrected away.
    n_het_artefactual = n_het_all - n_het_real
  ) %>%
  arrange(desc(n_het_real))

readr::write_tsv(coi, out_tsv)

n_poly <- sum(coi$infection_class == "polygenomic")
message("[coi] specimens: ", nrow(coi))
message("[coi] polygenomic (>=1 real het call): ", n_poly,
        sprintf(" (%.1f%%)", 100 * n_poly / nrow(coi)))
message("[coi] mono-compatible: ", nrow(coi) - n_poly)
message("[coi] NOTE: het calls attributable to artefacts and excluded: ",
        sum(coi$n_het_artefactual))
message("[coi] wrote ", out_tsv)
