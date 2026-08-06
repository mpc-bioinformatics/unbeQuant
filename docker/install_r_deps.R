#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Install the R dependencies required by ProtStatsWF::workflow_QC, then install
# ProtStatsWF itself from the vendored submodule at third_party/ProtStatsWF.
#
# Run during the Docker build. Exits non-zero if anything failed to install, so
# a broken image fails at build time rather than halfway through a pipeline run.
#
# Usage: Rscript install_r_deps.R [/path/to/ProtStatsWF]
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
protstats_path <- if (length(args) >= 1) args[1] else "/opt/ProtStatsWF"

# ProtStatsWF Imports available from CRAN (served as jammy binaries by P3M).
cran_pkgs <- c(
  "amap", "car", "circlize", "dendextend", "dplyr", "ggplot2", "ggrepel",
  "gtools", "magrittr", "matrixStats", "multcomp", "nlme", "openxlsx",
  "pbapply", "reshape2", "scales", "stringr", "tidyr",
  "BiocManager"
)

# ProtStatsWF Imports that live on Bioconductor. P3M does not mirror
# Bioconductor, so these are compiled from source (the slow part of the build).
bioc_pkgs <- c("affy", "ComplexHeatmap", "limma", "vsn")

message("== Installing CRAN dependencies ==")
install.packages(cran_pkgs)

message("== Installing Bioconductor dependencies ==")
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

message("== Installing ProtStatsWF from ", protstats_path, " ==")
if (!dir.exists(protstats_path)) {
  stop("Vendored ProtStatsWF not found at: ", protstats_path,
       "\nDid the submodule get checked out? Try: git submodule update --init --recursive")
}
install.packages(protstats_path, repos = NULL, type = "source")

# ── Verify everything landed, and fail the build if not ────────────────────
required <- c(cran_pkgs, bioc_pkgs, "ProtStatsWF")
installed <- rownames(installed.packages())
missing <- setdiff(required, installed)

if (length(missing) > 0) {
  message("\nERROR: the following packages failed to install:")
  for (p in missing) message("  - ", p)
  quit(status = 1)
}

# ComplexHeatmap >= 2.24.1 is an explicit ProtStatsWF requirement.
ch_version <- packageVersion("ComplexHeatmap")
if (ch_version < "2.24.1") {
  message("\nERROR: ComplexHeatmap ", ch_version, " is older than the required 2.24.1")
  quit(status = 1)
}

# Loading the package exercises the whole dependency graph.
suppressPackageStartupMessages(library(ProtStatsWF))
if (!exists("workflow_QC")) {
  message("\nERROR: ProtStatsWF loaded but workflow_QC is not available")
  quit(status = 1)
}

message("\n== All R dependencies installed ==")
message("  R              : ", as.character(getRversion()))
message("  ProtStatsWF    : ", as.character(packageVersion("ProtStatsWF")))
message("  ComplexHeatmap : ", as.character(ch_version))
message("  workflow_QC    : available")
