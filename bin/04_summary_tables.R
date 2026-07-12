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

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(purrr)
  library(tibble)
  library(viridis)
  library(reshape2)
})

# The run output directory is supplied by _setup.R (run_paths / run_output_dir),
# which Nextflow populates from NANORAVE_RUN_OUTPUT_DIR. The old filesystem-walk
# discovery (nr_resolve_run_dir + home_dir) is gone.
run_name_current <- MinION_run_name

run_paths <- run_paths[file.exists(run_paths)]
if (length(run_paths) == 0) {
  stop("Run output directory does not exist for: ", run_name_current)
}

# Output directories (04_summary block)
summary_root <- nr_block_paths(run_paths[[1]])$geno_summarised
summary_cov_dir <- file.path(summary_root, "coverage")
summary_cov_pivot_dir <- file.path(summary_cov_dir, "pivoted")
summary_cov_plot_dir <- file.path(summary_cov_dir, "plots")
summary_tab_dir <- file.path(summary_root, "tables")
summary_tab_plot_dir <- file.path(summary_tab_dir, "plots")
summary_geno_dir <- file.path(summary_root, "processed_genotypes")

for (d in c(summary_cov_dir, summary_cov_pivot_dir, summary_cov_plot_dir,
            summary_tab_dir, summary_tab_plot_dir, summary_geno_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

merged_data_amplicon <- data.frame()
merged_data_sample <- data.frame()
merged_dhfr_dhps_haplo_frequencies <- data.frame()
merged_dhfr_haplo_frequencies_ <- data.frame()
merged_dhps_haplo_frequencies_ <- data.frame()
merged_drug_resistance_frequencies_ <- data.frame()
merged_drug_resistance_variants_ <- data.frame()
merged_csp_nref_genotype_calls_long <- data.frame()
merged_genotype_calls_haplotypes_samplerows_v2 <- data.frame()
merged_k13_variant_info <- data.frame()
merged_variants_linked_to_metadata_intermediate_file <- data.frame()
merged_variants_linked_to_metadata_intermediate_file2 <- data.frame()

for (run_path in run_paths) {
  # Label the rows with the cohort name, not the output directory's basename.
  # Nextflow stages the run output as a directory literally called "run", so
  # basename() produced the string "run" in every run_name column.
  run_name_only <- MinION_run_name

  in_paths <- nr_block_paths(run_path)
  coverage_summary <- list.files(in_paths$quality, full.names = TRUE)
  processed_genotypes_summary <- list.files(in_paths$geno_calls, full.names = TRUE)
  tables_summary <- list.files(in_paths$geno_freq, full.names = TRUE)

  # -------------------------
  # Coverage summaries
  # -------------------------
  for (file in coverage_summary) {
    file_sec <- basename(file)

    if (grepl("_coverage_amplicon_run_summary_", file_sec)) {
      file_data <- read.csv(file, stringsAsFactors = FALSE)
      file_data$run_name <- run_name_only
      merged_data_amplicon <- bind_rows(merged_data_amplicon, file_data)

      write.csv(
        merged_data_amplicon,
        file.path(summary_cov_dir, "summary_coverage_amplicon_run_summary_.csv"),
        row.names = FALSE
      )

      pivot_data <- merged_data_amplicon %>%
        pivot_wider(names_from = run_name,
                    values_from = coverage_median:coverage_lower_10pcnt)
      write.csv(
        pivot_data,
        file.path(summary_cov_pivot_dir, "summary_coverage_amplicon_run_summary_.csv"),
        row.names = FALSE
      )

      plot_df <- merged_data_amplicon %>%
        mutate(
          coverage_median = suppressWarnings(as.numeric(coverage_median)),
          gene_target = as.character(gene_target)
        ) %>%
        filter(!is.na(coverage_median))

      if (nrow(plot_df) > 0) {
        p <- ggplot(plot_df, aes(x = run_name, y = coverage_median, fill = gene_target)) +
          geom_col(position = "dodge") +
          labs(x = "Run", y = "Median coverage", title = "Coverage amplicon run summary") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))

        ggsave(
          file.path(summary_cov_plot_dir, "coverage_amplicon_run_summary.jpeg"),
          plot = p,
          width = 10,
          height = 6
        )
      }
    } else if (grepl("_coverage_by_run_sample", file_sec)) {
      file_data <- read.csv(file, stringsAsFactors = FALSE) %>%
        mutate(
          sample_id = as.character(sample_id),
          run_name = run_name_only
        )

      merged_data_sample <- bind_rows(merged_data_sample, file_data)
      write.csv(
        merged_data_sample,
        file.path(summary_cov_dir, "summary_coverage_by_run_sample.csv"),
        row.names = FALSE
      )
    }
  }

  # -------------------------
  # Tables summaries
  # -------------------------
  for (file in tables_summary) {
    file_sec <- basename(file)

    if (grepl("dhfr_dhps_haplo_frequencies_", file_sec)) {
      file_data <- read.csv(file, stringsAsFactors = FALSE)
      file_data$run_name <- run_name_only
      merged_dhfr_dhps_haplo_frequencies <- bind_rows(merged_dhfr_dhps_haplo_frequencies, file_data)
      write.csv(
        merged_dhfr_dhps_haplo_frequencies,
        file.path(summary_tab_dir, "summary_dhfr_dhps_haplo_frequencies.csv"),
        row.names = FALSE
      )
    } else if (grepl("dhfr_haplo_frequencies_", file_sec)) {
      file_data <- read.csv(file, stringsAsFactors = FALSE)
      file_data$run_name <- run_name_only
      merged_dhfr_haplo_frequencies_ <- bind_rows(merged_dhfr_haplo_frequencies_, file_data)
      write.csv(
        merged_dhfr_haplo_frequencies_,
        file.path(summary_tab_dir, "summary_dhfr_haplo_frequencies.csv"),
        row.names = FALSE
      )
    } else if (grepl("dhps_haplo_frequencies_", file_sec)) {
      file_data <- read.csv(file, stringsAsFactors = FALSE)
      file_data$run_name <- run_name_only
      merged_dhps_haplo_frequencies_ <- bind_rows(merged_dhps_haplo_frequencies_, file_data)
      write.csv(
        merged_dhps_haplo_frequencies_,
        file.path(summary_tab_dir, "summary_dhps_haplo_frequencies.csv"),
        row.names = FALSE
      )
    } else if (grepl("drug_resistance_frequencies_", file_sec)) {
      file_data <- read.csv(file, stringsAsFactors = FALSE)
      file_data$run_name <- run_name_only
      merged_drug_resistance_frequencies_ <- bind_rows(merged_drug_resistance_frequencies_, file_data)
      write.csv(
        merged_drug_resistance_frequencies_,
        file.path(summary_tab_dir, "summary_drug_resistance_frequencies_.csv"),
        row.names = FALSE
      )
    } else if (grepl("drug_resistance_variants_", file_sec)) {
      file_data <- read.csv(file, stringsAsFactors = FALSE)
      file_data$run_name <- run_name_only
      merged_drug_resistance_variants_ <- bind_rows(merged_drug_resistance_variants_, file_data)
      write.csv(
        merged_drug_resistance_variants_,
        file.path(summary_tab_dir, "summary_drug_resistance_variants_.csv"),
        row.names = FALSE
      )

      if (all(c("CQ", "PYR", "SX", "SP.Rx", "SP.IPTp", "ART") %in% colnames(merged_drug_resistance_variants_))) {
        variant_counts <- merged_drug_resistance_variants_ %>%
          select(CQ, PYR, SX, SP.Rx, SP.IPTp, ART) %>%
          na.omit()

        if (nrow(variant_counts) > 0) {
          count_df <- data.frame(
            Feature = colnames(variant_counts),
            R = colSums(variant_counts == "R", na.rm = TRUE),
            S = colSums(variant_counts == "S", na.rm = TRUE)
          )

          data_melted <- reshape2::melt(count_df, id.vars = "Feature")
          colnames(data_melted)[2] <- "status"

          p <- ggplot(data_melted, aes(x = Feature, y = value, fill = status)) +
            geom_bar(stat = "identity") +
            labs(x = "", y = "Count", title = "All sample drug-response call summary") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))

          ggsave(
            file.path(summary_tab_plot_dir, "count_Summary_1_plot.jpeg"),
            plot = p,
            width = 8,
            height = 7
          )
        }
      }
    }
  }

  # -------------------------
  # Processed genotype summaries
  # -------------------------
  for (file in processed_genotypes_summary) {
    file_sec <- basename(file)

    if (grepl("_csp_nref_genotype_calls_long", file_sec)) {
      file_data <- read.csv(file, stringsAsFactors = FALSE) %>%
        mutate(
          sample_id = as.character(sample_id),
          run_name = run_name_only
        )
      merged_csp_nref_genotype_calls_long <- bind_rows(merged_csp_nref_genotype_calls_long, file_data)
      write.csv(
        merged_csp_nref_genotype_calls_long,
        file.path(summary_geno_dir, "summary_csp_nref_genotype_calls_long.csv"),
        row.names = FALSE
      )
    } else if (grepl("genotype_calls_haplotypes_samplerows_v2", file_sec)) {
      file_data <- read.csv(file, stringsAsFactors = FALSE)
      file_data$run_name <- run_name_only
      merged_genotype_calls_haplotypes_samplerows_v2 <- bind_rows(merged_genotype_calls_haplotypes_samplerows_v2, file_data)
      write.csv(
        merged_genotype_calls_haplotypes_samplerows_v2,
        file.path(summary_geno_dir, "summary_genotype_calls_haplotypes_samplerows_v2.csv"),
        row.names = FALSE
      )
    } else if (grepl("k13_variant_info", file_sec)) {
      file_data <- read.csv(file, stringsAsFactors = FALSE) %>%
        mutate(
          sample_id = as.character(sample_id),
          alt_nt = as.character(alt_nt),
          run_name = run_name_only
        )
      merged_k13_variant_info <- bind_rows(merged_k13_variant_info, file_data)
      write.csv(
        merged_k13_variant_info,
        file.path(summary_geno_dir, "summary_k13_variant_info.csv"),
        row.names = FALSE
      )
    } else if (grepl("variants_linked_to_metadata_intermediate_file\\.csv", file_sec)) {
      file_data <- read.csv(file, stringsAsFactors = FALSE)
      file_data$run_name <- run_name_only
      merged_variants_linked_to_metadata_intermediate_file <- bind_rows(merged_variants_linked_to_metadata_intermediate_file, file_data)
      write.csv(
        merged_variants_linked_to_metadata_intermediate_file,
        file.path(summary_geno_dir, "summary_variants_linked_to_metadata_intermediate_file.csv"),
        row.names = FALSE
      )
    } else if (grepl("variants_linked_to_metadata_intermediate_file2\\.csv", file_sec)) {
      file_data <- read.csv(file, stringsAsFactors = FALSE)
      file_data$run_name <- run_name_only
      merged_variants_linked_to_metadata_intermediate_file2 <- bind_rows(merged_variants_linked_to_metadata_intermediate_file2, file_data)
      write.csv(
        merged_variants_linked_to_metadata_intermediate_file2,
        file.path(summary_geno_dir, "summary_variants_linked_to_metadata_intermediate_file2.csv"),
        row.names = FALSE
      )
    }
  }
}
