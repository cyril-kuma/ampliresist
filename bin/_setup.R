# ============================================================================
# _setup.R  —  shared preamble for the drug-resistance analysis stages.
#
# Not an entrypoint. Each stage script (01_..04_) sources this first; it builds
# the shared context (paths, sample metadata, thresholds) from environment
# variables that Nextflow sets. There is no config file and no bash orchestration
# — Nextflow is the orchestrator.
#
# Contract (set by the AMPLICON_COVERAGE / GENOTYPE_CALLS / RESISTANCE_FREQUENCIES
# / SUMMARY_TABLES processes):
#
#   NANORAVE_RUN_NAME         run id; must equal ont_multiplex_group in the sheet
#   NANORAVE_ANALYSIS_DATE    stamp used in output filenames (YYYYMMDD)
#   NANORAVE_MIN_COV          coverage / per-SNP depth threshold (default 10x)
#   NANORAVE_INPUT_DIR        holds genome_coverage/ and variant_calling_unzip/
#                             (stages 1-2 only)
#   NANORAVE_RESOURCE_DIR     flat dir of analysis resources (assets/resources)
#   NANORAVE_METADATA_FILE    sample sheet workbook (.xlsx)
#   NANORAVE_MULTIPLEX_SHEET  sheet name inside that workbook
#   NANORAVE_RUN_OUTPUT_DIR   run output dir holding the numbered blocks
#
# NR_BIN is set by the calling stage script (it resolves its own location, since
# Nextflow puts bin/ on PATH rather than in the task working directory).
# ============================================================================
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

nr_env <- function(var, default = NULL) {
  val <- Sys.getenv(var, unset = "")
  if (!nzchar(val)) {
    if (is.null(default)) stop("Required environment variable not set: ", var)
    return(default)
  }
  val
}

if (!exists("NR_BIN")) {
  stop("NR_BIN is not set. _setup.R must be sourced from one of the numbered stage scripts.")
}

# Label used in output filenames. For a pooled cohort this is the cohort name
# (e.g. DRAG1_cohort), not a flowcell id.
MinION_run_name <- nr_env("NANORAVE_RUN_NAME")

# The sequencing runs whose samples belong to this analysis. Pipe-separated;
# defaults to the single run named above. Pooling every flowcell into one cohort
# is the normal mode — resistance-allele frequencies are a property of the sample
# set, not of a flowcell.
nanorave_runs <- strsplit(nr_env("NANORAVE_RUNS", MinION_run_name), "\\|")[[1]]
nanorave_runs <- trimws(nanorave_runs[nzchar(nanorave_runs)])

analysis_date   <- nr_env("NANORAVE_ANALYSIS_DATE", format(Sys.Date(), "%Y%m%d"))
multiplex_run   <- nr_env("NANORAVE_MULTIPLEX_SHEET", "samplesheet")

# ---------------------------------------------------------------------------
# Recover the sequencing run from a pipeline output filename.
#
# Files are named <run>_<barcodeNN>_<gene>.{bedGraph,vcf} and run names contain
# underscores (LUC_DRAG1_IMPAVESRUN23_20250920), so strip from "_barcode" right.
#
# This matters: barcode01 exists on EVERY flowcell. The analysis used to join
# samples on ont_barcode alone, which is only safe when the metadata is filtered
# to a single run. Pooling runs without also keying on the run would silently
# multiply rows and attach the wrong sample to a genotype.
# ---------------------------------------------------------------------------
nr_run_from_file <- function(paths) {
  sub("_barcode\\d+_.*$", "", basename(as.character(paths)))
}

# Coverage threshold. Lowered from 50x to 10x in 2026-07: the 50x bar came from
# the retired feasibility phase (Girgis et al.), where dhfr was callable in only
# 27.7% of samples.
min_cov_threshold <- as.numeric(nr_env("NANORAVE_MIN_COV", "10"))

# Stages 1-2 read the pipeline's alignment output; stages 3-4 only read the
# earlier blocks, so the input dir is optional here and each stage validates what
# it actually needs.
input_run_dir <- nr_env("NANORAVE_INPUT_DIR", "")
if (nzchar(input_run_dir)) input_run_dir <- normalizePath(input_run_dir, mustWork = FALSE)

# Flat resources directory (assets/resources): amplicon_gene_positions.csv,
# DR_variant_info_v2.xlsx, genetic_code.csv, k13_seq_annotated.csv and
# k13_who_artemisinin_partial_resistance_markers.csv.
resource_dir  <- normalizePath(nr_env("NANORAVE_RESOURCE_DIR"), mustWork = TRUE)
metadata_file <- normalizePath(nr_env("NANORAVE_METADATA_FILE"), mustWork = TRUE)

run_output_dir <- nr_env("NANORAVE_RUN_OUTPUT_DIR")
dir.create(run_output_dir, recursive = TRUE, showWarnings = FALSE)
run_output_dir <- normalizePath(run_output_dir, mustWork = TRUE)

# --- output block layout -----------------------------------------------------
source(file.path(NR_BIN, "_pipeline_paths.R"))
nr_paths <- nr_block_paths(run_output_dir)
nr_create_blocks(nr_paths)

# Stage-facing aliases, so the analysis code needs no path edits.
analysis_dr   <- nr_paths$quality      # stage 1 - amplicon coverage
analysis_dr_2 <- nr_paths$geno_calls   # stage 2 - genotype calls
analysis_dr_3 <- nr_paths$geno_freq    # stage 3 - frequencies (plots/ beneath)

# Stage 4 resolves the run dir itself; hand it the one we were given.
run_paths <- run_output_dir

# --- metadata ----------------------------------------------------------------
meta <- read_xlsx(metadata_file, sheet = multiplex_run) %>%
  filter(ont_multiplex_group %in% nanorave_runs)

if (nrow(meta) == 0) {
  stop("No sample-sheet rows matched ont_multiplex_group in {",
       paste(nanorave_runs, collapse = ", "), "} in ", metadata_file,
       " (sheet '", multiplex_run, "'). Each --input run_name must appear in the workbook.")
}

missing_runs <- setdiff(nanorave_runs, unique(meta$ont_multiplex_group))
if (length(missing_runs) > 0) {
  stop("These runs have no rows in the sample sheet: ",
       paste(missing_runs, collapse = ", "),
       ". Fix the workbook or --input; a run with no metadata would be silently dropped.")
}

# ---------------------------------------------------------------------------
# Controls (KH2, NC, PC) are re-sequenced on EVERY flowcell, so sample_id is not
# unique once the runs are pooled into one cohort.
#
# Every downstream stage groups and merges by sample_id, so a repeated id is
# corrosive:
#   * stage 1 merges six per-gene coverage tables by sample_id, so one control
#     present on 5 runs becomes 5^6 = 15,625 rows (a cartesian blow-up);
#   * stage 2's per-sample loops collapse the 5 instances into one and mangle the
#     variant list -- this is how kelch13 C580Y vanished from the KH2 control.
#
# Give every sequencing instance a unique id. The original is kept in
# sample_base, and nr_is_control() still recognises controls by their base id.
# ---------------------------------------------------------------------------
NR_CONTROL_IDS <- c("Control_KH2", "control_KH2", "KH2", "PC", "NC")

nr_base_sample_id <- function(x) sub("__.*$", "", as.character(x))
nr_is_control     <- function(x) nr_base_sample_id(x) %in% NR_CONTROL_IDS

meta$sample_base <- meta$sample_id
dup_ids <- names(which(table(meta$sample_id) > 1))
if (length(dup_ids) > 0) {
  message(sprintf("[setup] sample_id repeats across runs (%s) - disambiguating each sequencing instance with its run id",
                  paste(dup_ids, collapse = ", ")))
  idx <- meta$sample_id %in% dup_ids
  meta$sample_id[idx] <- paste0(meta$sample_id[idx], "__", meta$ont_multiplex_group[idx])
}

if (any(duplicated(meta$sample_id))) {
  stop("[setup] sample_id is still not unique after disambiguation: ",
       paste(unique(meta$sample_id[duplicated(meta$sample_id)]), collapse = ", "))
}

metadata <- meta

message(sprintf("[setup] cohort=%s  runs=%d (%s)  samples=%d  min_cov=%sx  out=%s",
                MinION_run_name, length(nanorave_runs),
                paste(nanorave_runs, collapse = ","),
                nrow(meta), min_cov_threshold, run_output_dir))
