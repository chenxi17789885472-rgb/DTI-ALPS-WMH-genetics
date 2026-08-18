#!/usr/bin/env Rscript

execute <- "--execute" %in% commandArgs(trailingOnly = TRUE)
bundle_root <- normalizePath(getwd(), mustWork = TRUE)
project_root <- file.path(bundle_root, "ALPS_GENETIC_ARCHITECTURE")
if (!dir.exists(project_root)) stop("Run this script from the reproducibility-package root.")
Sys.setenv(ALPS_PROJECT_ROOT = project_root)

analysis_scripts <- file.path(project_root, "06_scripts", c(
  "run_gate5b_coloc.R",
  "run_gate5c_independent_wmh_coloc.R",
  "run_gate5c_susie_coloc.R",
  "extract_gate5d_credible_variants.R",
  "run_malps_to_nonukb_wmh_mr.R"
))
figure_scripts <- file.path(project_root, "09_manuscript", "scripts", c(
  "build_figure1_sci_minimal.R",
  "build_main_figures_tables.R"
))
scripts <- c(analysis_scripts, figure_scripts)

missing_scripts <- scripts[!file.exists(scripts)]
if (length(missing_scripts)) stop("Missing scripts: ", paste(missing_scripts, collapse = ", "))
for (script in scripts) parse(script)

packages <- c(
  "coloc", "data.table", "ggplot2", "ggrepel", "jsonlite", "mr.raps",
  "MRPRESSO", "openxlsx", "patchwork", "png", "ragg", "scales",
  "susieR", "svglite", "TwoSampleMR"
)
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing locked packages; run renv::restore(): ", paste(missing_packages, collapse = ", "))
}

required_inputs <- c(
  file.path(project_root, "Gate5B", "02_extracted", "mALPS_Gate5B_regions.tsv.gz"),
  file.path(project_root, "Gate5B", "02_extracted", "WMH_chr16.tsv.gz"),
  file.path(project_root, "Gate5C", "01_WMH_independent", "Traylor2016_WMH_chr16.tsv.gz"),
  file.path(project_root, "Gate5C", "05_finemapping", "susie_harmonized_common.csv.gz"),
  file.path(project_root, "Gate5C", "05_finemapping", "1000G_EUR_chr16_susie.unphased.vcor1"),
  file.path(project_root, "Gate5C", "05_finemapping", "1000G_EUR_chr16_susie.unphased.vcor1.vars"),
  file.path(project_root, "Gate5C", "10_mALPS_to_WMH_MR", "02_harmonization", "mALPS_Traylor2016_WMH_harmonized.tsv")
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) stop("Missing frozen inputs: ", paste(missing_inputs, collapse = ", "))

message("Syntax, package, and frozen-input validation: PASS")
if (!execute) {
  message("Validation-only mode completed. Add --execute to rerun analyses and figures.")
  quit(status = 0)
}

rscript <- file.path(R.home("bin"), "Rscript")
log_dir <- file.path(bundle_root, "reproduction_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
for (script in scripts) {
  message("Running: ", basename(script))
  log_file <- file.path(log_dir, paste0(tools::file_path_sans_ext(basename(script)), ".log"))
  status <- system2(rscript, script, stdout = log_file, stderr = log_file)
  if (!identical(status, 0L)) stop("Script failed with status ", status, ": ", script)
}
message("Core R reproduction sequence: PASS")
