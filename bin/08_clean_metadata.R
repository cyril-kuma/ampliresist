#!/usr/bin/env Rscript
# ============================================================================
# 08_clean_metadata.R — generate the cleaned ecological metadata consumed by
# stage 18 (publication figures) directly from the authoritative workbook.
#
# The metadata WORKBOOK (samplesheet.xlsx, sheet 'samplesheet') is the single
# source of truth. samplesheet_clean_metadata.csv is a DERIVED product, produced
# here inside the workflow and propagated to stage 18 by channel — it is never a
# separately maintained input.
#
# Sample identifiers are taken from the pipeline's own genotype<->metadata map
# (produced by stage 15), so they match the analysed cohort exactly — including
# controls disambiguated as KH2__<run> — which guarantees the one-to-one join
# the figure script requires. Nothing is silently dropped: analysed samples with
# no ecological row in the workbook are retained with NA ecology and reported.
#
# Usage:
#   08_clean_metadata.R <summary_block_dir> <samplesheet.xlsx> <sheet> <out.csv>
# ============================================================================
suppressPackageStartupMessages({library(readxl); library(dplyr); library(readr)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: 08_clean_metadata.R <summary_block_dir> <samplesheet.xlsx> <sheet> <out.csv>")
}
summary_block <- args[[1]]; xlsx <- args[[2]]; sheet <- args[[3]]; out <- args[[4]]

map_path <- file.path(summary_block, "processed_genotypes",
                      "summary_variants_linked_to_metadata_intermediate_file.csv")
if (!file.exists(map_path)) stop("Cannot find genotype-metadata map: ", map_path)
if (!file.exists(xlsx))     stop("Cannot find metadata workbook: ", xlsx)

map <- readr::read_csv(map_path, show_col_types = FALSE)
key <- c("ont_multiplex_group", "ont_barcode")
for (col in c("sample_id", key)) {
  if (!col %in% names(map)) stop("Genotype-metadata map is missing column: ", col)
}
map <- map |> select(all_of(c("sample_id", key))) |> distinct()
if (any(duplicated(map$sample_id))) {
  stop("Analysed sample_id is not unique in the map file; cannot build a 1:1 clean metadata.")
}

req_eco <- c("bioclimatic_zone", "sibling_species", "collection_site",
             "latitude", "longitude", "host_feeding_type", "pf_ct")
ss <- readxl::read_xlsx(xlsx, sheet = sheet)
missing_cols <- setdiff(c(key, req_eco, "sample_status"), names(ss))
if (length(missing_cols) > 0) {
  stop("Workbook '", basename(xlsx), "' (sheet '", sheet,
       "') is missing required columns: ", paste(missing_cols, collapse = ", "))
}

eco <- ss |>
  select(all_of(c(key, req_eco, "sample_status"))) |>
  distinct(across(all_of(key)), .keep_all = TRUE)

clean <- map |> left_join(eco, by = key)

unmatched <- clean |> filter(is.na(bioclimatic_zone))
if (nrow(unmatched) > 0) {
  message(sprintf(
    "[clean_metadata] %d analysed sample(s) have no ecological row in the workbook (kept with NA ecology): %s",
    nrow(unmatched), paste(head(unmatched$sample_id, 10), collapse = ", ")))
}

clean <- clean |> select(sample_id, all_of(req_eco), sample_status)
readr::write_csv(clean, out)
message(sprintf(
  "[clean_metadata] wrote %s: %d row(s), %d with ecology, from workbook %s (sheet '%s')",
  out, nrow(clean), sum(!is.na(clean$bioclimatic_zone)), basename(xlsx), sheet))
