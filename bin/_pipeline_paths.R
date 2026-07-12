# ============================================================================
# _pipeline_paths.R  —  SINGLE SOURCE OF TRUTH for the downstream output layout.
#
# Every stage resolves its output/input subdirectories through nr_block_paths().
# To rename or renumber a block, edit it HERE ONCE; the whole pipeline follows.
#
# Layout inside a run output directory — one numbered block per stage, so the
# directory order mirrors the order the stages run in:
#
#   01_coverage_qc/      stage 1 - per-amplicon coverage summaries
#   02_genotype_calls/   stage 2 - DR genotype + k13 variant calls
#   03_frequencies/      stage 3 - DR + haplotype frequency tables
#     plots/                       and their plots
#   04_summary/          stage 4 - merged summary tables
#     coverage/ processed_genotypes/ tables/
#
# REMOVED 2026-07: the feasibility block (formerly 03_feasibility). It existed to
# answer "can resistance alleles be genotyped from Pf-positive mosquito blood
# meals at all?" — now settled, so the stage and its outputs are gone. The
# publication-figure and 10x-sensitivity blocks went with it; their scripts were
# never present in this repo and always logged [SKIP].
# ============================================================================

nr_block_paths <- function(run_dir) {
  list(
    run             = run_dir,
    manifest        = file.path(run_dir, "00_run_manifest.txt"),

    # 01 - coverage QC (stage 1). No plots/ block: stage 1 writes only tables, and
    # creating the directory left an empty plots/ in every run's output.
    quality         = file.path(run_dir, "01_coverage_qc"),

    # 02 - genotype calls (stage 2)
    geno_calls      = file.path(run_dir, "02_genotype_calls"),

    # 03 - frequencies + plots (stage 3)
    geno_freq       = file.path(run_dir, "03_frequencies"),
    geno_plots      = file.path(run_dir, "03_frequencies", "plots"),

    # 04 - merged summaries (stage 4)
    geno_summarised = file.path(run_dir, "04_summary")
  )
}

# Create all leaf output directories for a run (idempotent).
nr_create_blocks <- function(p) {
  dirs <- c(p$quality,
            p$geno_calls,
            p$geno_freq, p$geno_plots,
            p$geno_summarised)
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(p)
}
