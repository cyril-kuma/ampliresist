#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(tidyr)
  library(vcfR)
  library(ggplot2)
  library(pegas)
})

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    if (i == length(args)) stop("Missing value for argument: ", key)
    out[[sub("^--", "", key)]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

required_arg <- function(x, name) {
  if (is.null(x) || !nzchar(x)) stop("Missing required argument --", name)
  x
}

read_metadata <- function(path, sheet) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    as.data.frame(readxl::read_xlsx(path, sheet = sheet))
  }
}

norm_group <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x == "Middlebelt" ~ "Middle Belt",
    x == "Middle Belt" ~ "Middle Belt",
    x == "Coastal" ~ "Coastal",
    x == "Savannah" ~ "Savannah",
    TRUE ~ x
  )
}

safe_gt_to_allele <- function(gt, ref, alt) {
  if (is.na(gt) || gt == "" || gt == "." || gt == "./." || gt == ".|.") return(NA_character_)
  gt_first <- strsplit(gt, ":", fixed = TRUE)[[1]][1]
  if (grepl("\\.", gt_first)) return(NA_character_)
  tokens <- strsplit(gt_first, "[/|]")[[1]]
  idx <- suppressWarnings(as.integer(tokens))
  if (any(is.na(idx))) return(NA_character_)
  if (all(idx == 0)) return(ref)
  nz <- idx[idx > 0]
  if (length(nz) == 0) return(ref)
  alt_values <- strsplit(alt, ",", fixed = TRUE)[[1]]
  alt_idx <- nz[[1]]
  if (alt_idx <= length(alt_values)) return(alt_values[[alt_idx]])
  "ALT"
}

hap_diversity <- function(h) {
  h <- h[!is.na(h) & h != ""]
  n <- length(h)
  if (n <= 1) return(NA_real_)
  p <- table(h) / n
  (n / (n - 1)) * (1 - sum(p^2))
}

make_haplotype <- function(row_vec) {
  if (length(row_vec) == 0) return("No_polymorphism")
  if (all(is.na(row_vec))) return("Missing")
  row_vec[is.na(row_vec)] <- "N"
  paste(row_vec, collapse = "")
}

make_diagnostic_plot <- function(path, title, detail) {
  png(path, width = 1200, height = 1200)
  plot.new()
  title(main = title, col.main = "black")
  text(0.5, 0.6, detail, cex = 1.2)
  text(0.5, 0.45, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), cex = 0.9)
  dev.off()
}

make_diagnostic_pdf <- function(path, title, detail) {
  grDevices::cairo_pdf(path, width = 9, height = 7)
  plot.new()
  title(main = title, col.main = "black")
  text(0.5, 0.6, detail, cex = 1.2)
  text(0.5, 0.45, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), cex = 0.9)
  dev.off()
}

plot_pca <- function(df, xcol, ycol, out_png, out_pdf, title) {
  if (nrow(df) == 0) {
    make_diagnostic_plot(out_png, title, "No PCA rows available after group join")
    make_diagnostic_pdf(out_pdf, title, "No PCA rows available after group join")
    return(invisible(NULL))
  }

  gp <- ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]], color = group)) +
    geom_point(alpha = 0.85, size = 2) +
    scale_color_manual(values = c("Savannah" = "gray20", "Coastal" = "#1f77b4", "Middle Belt" = "#d62728", "Unknown" = "#7f7f7f")) +
    theme_classic() +
    labs(title = title, x = xcol, y = ycol, color = "Group")

  ggsave(out_png, gp, width = 9, height = 7, dpi = 300)
  ggsave(out_pdf, gp, width = 9, height = 7, device = cairo_pdf)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

merged_vcf <- required_arg(args[["merged-vcf"]], "merged-vcf")
metadata_file <- required_arg(args[["metadata-file"]], "metadata-file")
metadata_sheet <- ifelse(is.null(args[["metadata-sheet"]]), "Sheet1", args[["metadata-sheet"]])
maf_dir <- required_arg(args[["maf-dir"]], "maf-dir")
geneflow_dir <- required_arg(args[["geneflow-dir"]], "geneflow-dir")
report_dir <- required_arg(args[["output-report-dir"]], "output-report-dir")
run_name <- ifelse(is.null(args[["run-name"]]), "unknown_run", args[["run-name"]])

th2r_start <- as.integer(ifelse(is.null(args[["th2r-start"]]), 940L, args[["th2r-start"]]))
th2r_end <- as.integer(ifelse(is.null(args[["th2r-end"]]), 988L, args[["th2r-end"]]))
th3r_start <- as.integer(ifelse(is.null(args[["th3r-start"]]), 1030L, args[["th3r-start"]]))
th3r_end <- as.integer(ifelse(is.null(args[["th3r-end"]]), 1099L, args[["th3r-end"]]))

if (!file.exists(merged_vcf)) stop("Merged VCF not found: ", merged_vcf)
if (!dir.exists(maf_dir)) stop("MAF directory not found: ", maf_dir)
if (!dir.exists(geneflow_dir)) stop("Gene-flow directory not found: ", geneflow_dir)
if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)

data_dir <- file.path(report_dir, "data")
figures_dir <- file.path(report_dir, "figures")
tables_dir <- file.path(report_dir, "tables")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

sample_group_map_path <- file.path(data_dir, "csp_sample_group_map.tsv")
variant_matrix_wide_path <- file.path(data_dir, "csp_variant_matrix_wide.tsv")
variant_matrix_meta_path <- file.path(data_dir, "csp_variant_matrix_with_metadata.tsv")
haplotype_profile_path <- file.path(data_dir, "csp_haplotype_profiles.tsv")
site_catalog_path <- file.path(data_dir, "csp_window_site_catalog.tsv")
summary_metrics_path <- file.path(tables_dir, "csp_haplotype_diversity_metrics.tsv")
hap_freq_group_path <- file.path(tables_dir, "csp_combined_haplotype_frequency_by_group.tsv")
network_png_path <- file.path(figures_dir, "csp_th2r_th3r_haplotype_network.png")
maf_hist_png <- file.path(figures_dir, "csp_maf_distribution_by_group.png")
pca21_png <- file.path(figures_dir, "pca_pc2_vs_pc1.png")
pca21_pdf <- file.path(figures_dir, "pca_pc2_vs_pc1.pdf")
pca41_png <- file.path(figures_dir, "pca_pc4_vs_pc1.png")
pca41_pdf <- file.path(figures_dir, "pca_pc4_vs_pc1.pdf")
pi_png <- file.path(figures_dir, "nucleotide_diversity_pi_by_group.png")
tajima_png <- file.path(figures_dir, "tajimas_d_by_group.png")
fst_png <- file.path(figures_dir, "pairwise_fst_distribution.png")

meta <- read_metadata(metadata_file, metadata_sheet)
if (!("sample_id" %in% names(meta))) stop("metadata must contain sample_id")
if (!("collection_region" %in% names(meta))) stop("metadata must contain collection_region")
meta$sample_id <- as.character(meta$sample_id)

filtered_samples_path <- file.path(maf_dir, "filtered_samples.tsv")
if (file.exists(filtered_samples_path)) {
  fs <- readr::read_tsv(filtered_samples_path, show_col_types = FALSE)
  if (!("sample_id" %in% names(fs))) stop("filtered_samples.tsv missing sample_id")
  group_col <- if ("group" %in% names(fs)) "group" else if ("collection_region" %in% names(fs)) "collection_region" else NA_character_
  if (is.na(group_col)) stop("filtered_samples.tsv missing group and collection_region")
  sample_group <- fs %>%
    mutate(sample_id = as.character(sample_id), group = norm_group(.data[[group_col]])) %>%
    select(sample_id, group) %>%
    distinct(sample_id, .keep_all = TRUE)
} else {
  sample_group <- meta %>%
    transmute(sample_id = as.character(sample_id), group = norm_group(collection_region)) %>%
    distinct(sample_id, .keep_all = TRUE)
}

readr::write_tsv(sample_group, sample_group_map_path)

message("Reading merged VCF: ", merged_vcf)
vcf <- vcfR::read.vcfR(merged_vcf, verbose = FALSE)
fix <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
gt <- vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE)

if (is.null(dim(gt)) || nrow(gt) == 0 || ncol(gt) == 0) {
  stop("No genotype matrix in merged VCF.")
}

samples <- colnames(gt)
contigs <- unique(fix$CHROM)
csp_contigs <- contigs[grepl("csp", contigs, ignore.case = TRUE)]
if (length(csp_contigs) == 0) stop("No CSP-like contig found in merged VCF CHROM values")
csp_contig <- csp_contigs[[1]]

idx_csp <- which(fix$CHROM == csp_contig)
if (length(idx_csp) == 0) stop("No records found for contig: ", csp_contig)

csp_fix <- fix[idx_csp, , drop = FALSE]
csp_gt <- gt[idx_csp, , drop = FALSE]
csp_fix$POS <- as.integer(csp_fix$POS)

site_window <- dplyr::case_when(
  csp_fix$POS >= th2r_start & csp_fix$POS <= th2r_end ~ "TH2R",
  csp_fix$POS >= th3r_start & csp_fix$POS <= th3r_end ~ "TH3R",
  TRUE ~ "Outside"
)

idx_window <- which(site_window %in% c("TH2R", "TH3R"))
if (length(idx_window) == 0) {
  make_diagnostic_plot(
    network_png_path,
    "CSP TH2R/TH3R Haplotype Network",
    "No SNPs found in TH2R/TH3R coordinate windows"
  )
  empty_summary <- tibble(
    run_name = run_name,
    csp_contig = csp_contig,
    n_samples = length(samples),
    n_window_sites = 0L,
    n_polymorphic_window_sites = 0L,
    n_unique_haplotypes = 0L,
    overall_haplotype_diversity = NA_real_
  )
  readr::write_tsv(empty_summary, summary_metrics_path)
  quit(save = "no", status = 0)
}

win_fix <- csp_fix[idx_window, , drop = FALSE]
win_gt <- csp_gt[idx_window, , drop = FALSE]
win_label <- site_window[idx_window]

is_poly <- apply(win_gt, 1, function(g) any(grepl("[1-9]", g)))
poly_fix <- win_fix[is_poly, , drop = FALSE]
poly_gt <- win_gt[is_poly, , drop = FALSE]
poly_label <- win_label[is_poly]

# The merged cohort VCF can carry more than one record for the same site
# (identical POS/REF/ALT). Sites become matrix column names below
# (csp_<POS>_<REF><ALT>), so duplicates produce duplicate columns and the
# left_join() on the variant matrix fails with "Input columns must be unique".
# They are the same variant, so collapse them.
dup_site <- duplicated(paste(poly_fix$POS, poly_fix$REF, poly_fix$ALT))
if (any(dup_site)) {
  message(sprintf("[csp] collapsing %d duplicate site record(s) in the merged VCF",
                  sum(dup_site)))
  poly_fix <- poly_fix[!dup_site, , drop = FALSE]
  poly_gt <- poly_gt[!dup_site, , drop = FALSE]
  poly_label <- poly_label[!dup_site]
}

if (nrow(poly_fix) == 0) {
  make_diagnostic_plot(
    network_png_path,
    "CSP TH2R/TH3R Haplotype Network",
    "No polymorphic SNPs in TH2R/TH3R windows"
  )

  empty_summary <- tibble(
    run_name = run_name,
    csp_contig = csp_contig,
    n_samples = length(samples),
    n_window_sites = nrow(win_fix),
    n_polymorphic_window_sites = 0L,
    n_unique_haplotypes = 0L,
    overall_haplotype_diversity = NA_real_
  )
  readr::write_tsv(empty_summary, summary_metrics_path)
  quit(save = "no", status = 0)
}

allele_matrix <- matrix(NA_character_, nrow = length(samples), ncol = nrow(poly_fix),
                        dimnames = list(samples, paste0("csp_", poly_fix$POS, "_", poly_fix$REF, poly_fix$ALT)))

for (i in seq_len(nrow(poly_fix))) {
  ref <- poly_fix$REF[[i]]
  alt <- poly_fix$ALT[[i]]
  for (j in seq_along(samples)) {
    allele_matrix[j, i] <- safe_gt_to_allele(poly_gt[i, j], ref, alt)
  }
}

site_catalog <- tibble(
  site_id = colnames(allele_matrix),
  contig = csp_contig,
  position = poly_fix$POS,
  ref = poly_fix$REF,
  alt = poly_fix$ALT,
  window = poly_label
) %>% arrange(position)

readr::write_tsv(site_catalog, site_catalog_path)

mat_df <- as.data.frame(allele_matrix, stringsAsFactors = FALSE)
mat_df$sample_id <- rownames(allele_matrix)
mat_df <- mat_df %>%
  left_join(sample_group, by = "sample_id") %>%
  mutate(group = ifelse(is.na(group) | group == "", "Unknown", group)) %>%
  relocate(sample_id, group)

readr::write_tsv(mat_df, variant_matrix_wide_path)

meta_keep_cols <- intersect(
  c("sample_id", "collection_region", "collection_site", "bioclimatic_zone", "ont_barcode"),
  names(meta)
)

mat_meta <- mat_df %>%
  left_join(meta %>% select(all_of(meta_keep_cols)), by = "sample_id")

relocate_cols <- intersect(
  c("sample_id", "ont_barcode", "group", "collection_region", "collection_site", "bioclimatic_zone"),
  names(mat_meta)
)
mat_meta <- mat_meta %>% relocate(all_of(relocate_cols))

readr::write_tsv(mat_meta, variant_matrix_meta_path)

th2r_cols <- site_catalog %>% filter(window == "TH2R") %>% pull(site_id)
th3r_cols <- site_catalog %>% filter(window == "TH3R") %>% pull(site_id)

hap_profile <- mat_df %>%
  rowwise() %>%
  mutate(
    th2r_haplotype = make_haplotype(c_across(all_of(th2r_cols))),
    th3r_haplotype = make_haplotype(c_across(all_of(th3r_cols))),
    combined_haplotype = paste0(th2r_haplotype, "|", th3r_haplotype)
  ) %>%
  ungroup()

readr::write_tsv(hap_profile, haplotype_profile_path)

freq_group <- hap_profile %>%
  count(group, combined_haplotype, name = "n") %>%
  group_by(group) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup() %>%
  arrange(group, desc(n))

readr::write_tsv(freq_group, hap_freq_group_path)

metrics <- hap_profile %>%
  group_by(group) %>%
  summarise(
    n_samples = n(),
    n_unique_haplotypes = n_distinct(combined_haplotype),
    haplotype_diversity = hap_diversity(combined_haplotype),
    .groups = "drop"
  )

overall <- tibble(
  group = "Overall",
  n_samples = nrow(hap_profile),
  n_unique_haplotypes = n_distinct(hap_profile$combined_haplotype),
  haplotype_diversity = hap_diversity(hap_profile$combined_haplotype)
)
metrics <- bind_rows(metrics, overall)
readr::write_tsv(metrics, summary_metrics_path)

hap_matrix <- hap_profile %>%
  select(all_of(site_catalog$site_id)) %>%
  mutate(across(everything(), ~ ifelse(is.na(.x), "N", .x)))

unique_hap <- n_distinct(hap_profile$combined_haplotype)

if (nrow(hap_matrix) < 2 || ncol(hap_matrix) < 1 || unique_hap < 2) {
  make_diagnostic_plot(
    network_png_path,
    "CSP TH2R/TH3R Haplotype Network",
    paste0("Insufficient haplotype diversity for network (unique haplotypes = ", unique_hap, ")")
  )
} else {
  tryCatch({
    h <- pegas::haplotype(as.matrix(hap_matrix))
    net <- pegas::haploNet(h)

    idx <- attr(h, "index")
    hap_name <- rep(names(idx), lengths(idx))
    sample_idx <- unlist(idx)
    grp <- hap_profile$group[sample_idx]
    group_levels <- sort(unique(hap_profile$group))

    group_palette <- c(
      "Savannah" = "gray20",
      "Coastal" = "#1f77b4",
      "Middle Belt" = "#d62728",
      "Middlebelt" = "#d62728",
      "Unknown" = "#7f7f7f"
    )
    if (!all(group_levels %in% names(group_palette))) {
      missing_groups <- setdiff(group_levels, names(group_palette))
      extra_cols <- grDevices::rainbow(length(missing_groups))
      names(extra_cols) <- missing_groups
      group_palette <- c(group_palette, extra_cols)
    }

    pie <- table(
      factor(hap_name, levels = names(idx)),
      factor(grp, levels = group_levels)
    )

    png(network_png_path, width = 1200, height = 1200)
    par(mar = c(5.1, 4.1, 5.1, 2.1))
    plot(
      net,
      size = attr(net, "freq"),
      pie = pie,
      bg = unname(group_palette[group_levels]),
      scale.ratio = 0.5,
      cex = 0.9
    )
    title(main = "CSP TH2R/TH3R SNP Haplotype Network")
    mtext(
      side = 1,
      line = 3,
      cex = 0.9,
      text = paste0(
        "Samples=", nrow(hap_profile),
        " | Polymorphic sites=", nrow(poly_fix),
        " | Unique haplotypes=", unique_hap
      )
    )
    legend(
      "topleft",
      legend = group_levels,
      fill = unname(group_palette[group_levels]),
      border = unname(group_palette[group_levels]),
      bty = "n",
      cex = 0.9
    )
    dev.off()
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    make_diagnostic_plot(network_png_path, "CSP TH2R/TH3R Haplotype Network", paste("Network plotting failed:", conditionMessage(e)))
  })
}

maf_file <- file.path(maf_dir, "combined_maf.tsv")
if (file.exists(maf_file)) {
  maf_df <- readr::read_tsv(maf_file, show_col_types = FALSE) %>%
    mutate(Group = norm_group(Group))

  p_maf <- ggplot(maf_df, aes(x = MAF, fill = Group)) +
    geom_histogram(binwidth = 0.05, position = "dodge", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = c("Savannah" = "gray20", "Coastal" = "#1f77b4", "Middle Belt" = "#d62728")) +
    theme_classic() +
    labs(title = "MAF Distribution by Population Group", x = "Minor Allele Frequency", y = "Variant Count")

  ggsave(maf_hist_png, p_maf, width = 9, height = 7, dpi = 300)
}

pca_file <- file.path(geneflow_dir, "combined_pca.eigenvec")
if (file.exists(pca_file)) {
  pca <- read.table(pca_file, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(pca) >= 7) {
    colnames(pca)[1:7] <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5")
    pca <- pca %>%
      mutate(IID = as.character(IID)) %>%
      left_join(sample_group %>% mutate(group = norm_group(group)), by = c("IID" = "sample_id")) %>%
      mutate(group = ifelse(is.na(group) | group == "", "Unknown", group))

    plot_pca(pca, "PC2", "PC1", pca21_png, pca21_pdf, "PCA: PC2 vs PC1")
    plot_pca(pca, "PC4", "PC1", pca41_png, pca41_pdf, "PCA: PC4 vs PC1")
  }
}

pi_files <- list.files(geneflow_dir, pattern = "_pi\\.windowed\\.pi$", full.names = TRUE)
if (length(pi_files) > 0) {
  pi_df <- bind_rows(lapply(pi_files, function(f) {
    x <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
    x$group <- norm_group(sub("_pi\\.windowed\\.pi$", "", basename(f)))
    x
  }))

  if (all(c("BIN_START", "PI", "group") %in% names(pi_df))) {
    p_pi <- ggplot(pi_df, aes(x = BIN_START, y = PI, color = group)) +
      geom_line(linewidth = 0.8) +
      theme_classic() +
      scale_color_manual(values = c("savannah" = "gray20", "coastal" = "#1f77b4", "middlebelt" = "#d62728", "Savannah" = "gray20", "Coastal" = "#1f77b4", "Middle Belt" = "#d62728")) +
      labs(title = "Nucleotide Diversity (Pi)", x = "Position", y = "Pi", color = "Group")
    ggsave(pi_png, p_pi, width = 10, height = 6, dpi = 300)
  }
}

tajima_files <- list.files(geneflow_dir, pattern = "_tajimaD\\.Tajima\\.D$", full.names = TRUE)
if (length(tajima_files) > 0) {
  tj_df <- bind_rows(lapply(tajima_files, function(f) {
    x <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
    x$group <- norm_group(sub("_tajimaD\\.Tajima\\.D$", "", basename(f)))
    x
  }))

  if (all(c("BIN_START", "TajimaD", "group") %in% names(tj_df))) {
    p_tj <- ggplot(tj_df, aes(x = BIN_START, y = TajimaD, color = group)) +
      geom_line(linewidth = 0.8) +
      theme_classic() +
      scale_color_manual(values = c("savannah" = "gray20", "coastal" = "#1f77b4", "middlebelt" = "#d62728", "Savannah" = "gray20", "Coastal" = "#1f77b4", "Middle Belt" = "#d62728")) +
      labs(title = "Tajima's D", x = "Position", y = "Tajima's D", color = "Group")
    ggsave(tajima_png, p_tj, width = 10, height = 6, dpi = 300)
  }
}

fst_files <- c(
  file.path(geneflow_dir, "savannah_vs_coastal_fst.weir.fst"),
  file.path(geneflow_dir, "savannah_vs_middlebelt_fst.weir.fst"),
  file.path(geneflow_dir, "coastal_vs_middlebelt_fst.weir.fst")
)

if (all(file.exists(fst_files))) {
  fst_df <- bind_rows(
    read.table(fst_files[[1]], header = TRUE, stringsAsFactors = FALSE) %>% mutate(comparison = "Savannah vs Coastal"),
    read.table(fst_files[[2]], header = TRUE, stringsAsFactors = FALSE) %>% mutate(comparison = "Savannah vs Middle Belt"),
    read.table(fst_files[[3]], header = TRUE, stringsAsFactors = FALSE) %>% mutate(comparison = "Coastal vs Middle Belt")
  ) %>%
    filter(!is.na(WEIR_AND_COCKERHAM_FST))

  p_fst <- ggplot(fst_df, aes(x = comparison, y = WEIR_AND_COCKERHAM_FST, fill = comparison)) +
    geom_violin(trim = FALSE, alpha = 0.8) +
    geom_boxplot(width = 0.1, outlier.size = 0.4, color = "black") +
    theme_classic() +
    scale_fill_manual(values = c("Savannah vs Coastal" = "#1f77b4", "Savannah vs Middle Belt" = "gray40", "Coastal vs Middle Belt" = "#d62728")) +
    labs(title = "Pairwise Fst Distributions", x = "Comparison", y = "Weir-Cockerham Fst", fill = "Comparison")

  ggsave(fst_png, p_fst, width = 10, height = 6, dpi = 300)
}

run_summary <- tibble(
  run_name = run_name,
  csp_contig = csp_contig,
  n_samples = nrow(hap_profile),
  n_window_sites = nrow(win_fix),
  n_polymorphic_window_sites = nrow(poly_fix),
  n_unique_haplotypes = n_distinct(hap_profile$combined_haplotype),
  overall_haplotype_diversity = hap_diversity(hap_profile$combined_haplotype)
)

readr::write_tsv(run_summary, file.path(tables_dir, "csp_report_run_summary.tsv"))

message("Prepared population genetics CSP report assets in: ", report_dir)
