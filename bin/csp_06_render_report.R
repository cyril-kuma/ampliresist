#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rmarkdown)
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

args <- parse_args(commandArgs(trailingOnly = TRUE))

report_dir <- required_arg(args[["output-report-dir"]], "output-report-dir")
run_name <- ifelse(is.null(args[["run-name"]]), "unknown_run", args[["run-name"]])
script_dir <- ifelse(is.null(args[["script-dir"]]), getwd(), args[["script-dir"]])

template_rmd <- file.path(script_dir, "csp_07_report_template.Rmd")
if (!file.exists(template_rmd)) stop("Template not found: ", template_rmd)

data_dir <- file.path(report_dir, "data")
figures_dir <- file.path(report_dir, "figures")
tables_dir <- file.path(report_dir, "tables")

if (!dir.exists(data_dir)) stop("Expected data directory not found: ", data_dir)
if (!dir.exists(figures_dir)) stop("Expected figures directory not found: ", figures_dir)
if (!dir.exists(tables_dir)) stop("Expected tables directory not found: ", tables_dir)

out_file <- "population_genetics_csp_report.html"

rmarkdown::render(
  input = template_rmd,
  output_file = out_file,
  output_dir = report_dir,
  # rmarkdown knits in the TEMPLATE's directory by default. The template lives in
  # bin/ (which Nextflow puts on PATH, not in the task working directory), so the
  # report's relative paths - report/tables/... - would resolve against bin/ and
  # nothing would be found. Knit where the data actually is.
  knit_root_dir = normalizePath(getwd()),
  params = list(
    run_name = run_name,
    report_dir = report_dir,
    data_dir = data_dir,
    figures_dir = figures_dir,
    tables_dir = tables_dir
  ),
  envir = new.env(parent = globalenv())
)

message("Rendered report: ", file.path(report_dir, out_file))
