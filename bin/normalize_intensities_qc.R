#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Normalize an unbeQuant quantification table with ProtStatsWF::workflow_QC.
#
# Input is the table produced by bin/prepare_qc_input.py: an ID column
# ('openms_ceid') followed by the intensity columns, renamed to
# <CONDITION>_<REPLICATE> so that workflow_QC can recover the sample grouping
# (it derives the group as the text before the first underscore).
#
# The workflow_QC call below follows the package README verbatim, except for
# the data-form arguments that must describe our input/output:
#   filetype           "xlsx" -> "tsv"   (our table is tab-separated text)
#   sep                ","    -> "\t"
#   output_type        "xlsx" -> "tsv"   (xlsx is write-only in this repo and
#                                         caps at 1,048,576 rows)
#   data_path / output_path / intensity_columns  -> supplied by the caller
# Every other argument is left at its package default, which already matches
# the README (do_log_transformation = TRUE, log_base = 2, zero_to_NA = TRUE,
# use_groups = TRUE, plot = TRUE).
#
# Outputs written to --output_path by workflow_QC:
#   D_norm_ID<suffix>.tsv    ID column + normalized, log2-transformed intensities
#   D_norm_wide<suffix>.csv, D_norm_long<suffix>.csv
#   boxplot / PCA / valid-value plots (plot_device = "png")
#
# Usage:
#   Rscript normalize_intensities_qc.R \
#       --data_path=qc_input.tsv \
#       --output_path=qc_out \
#       --intensity_col_start=2 --intensity_col_end=26
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(ProtStatsWF))

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}

parse_logical <- function(name, default) {
  raw <- parse_arg(name, NULL)
  if (is.null(raw)) return(default)
  val <- as.logical(raw)
  if (is.na(val)) {
    stop("--", name, " must be TRUE or FALSE, got: ", raw)
  }
  val
}

parse_int <- function(name) {
  raw <- parse_arg(name, NULL)
  if (is.null(raw)) return(NA_integer_)
  val <- suppressWarnings(as.integer(raw))
  if (is.na(val)) stop("--", name, " must be an integer, got: ", raw)
  val
}

data_path   <- parse_arg("data_path")
output_path <- parse_arg("output_path")
int_start   <- parse_int("intensity_col_start")
int_end     <- parse_int("intensity_col_end")
method      <- parse_arg("normalization_method", "loess")
suffix      <- parse_arg("suffix", "")
gen_ma      <- parse_logical("generate_MAplots", TRUE)

# ── Validate ───────────────────────────────────────────────────────────────
if (is.null(data_path) || is.null(output_path) ||
    is.na(int_start) || is.na(int_end)) {
  stop("Required arguments: --data_path, --output_path, ",
       "--intensity_col_start, --intensity_col_end")
}
if (!file.exists(data_path)) {
  stop("Input table not found: ", data_path)
}
if (int_start < 2) {
  stop("--intensity_col_start must be >= 2 (column 1 is the ID column), got: ",
       int_start)
}
if (int_end < int_start) {
  stop("--intensity_col_end (", int_end, ") must be >= --intensity_col_start (",
       int_start, ")")
}

# workflow_QC appends `suffix` to every output file name and requires it to
# begin with an underscore when non-empty.
if (nzchar(suffix) && !startsWith(suffix, "_")) {
  suffix <- paste0("_", suffix)
}

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

message("======================================================================")
message("ProtStatsWF::workflow_QC normalization")
message("======================================================================")
message("  data_path            : ", data_path)
message("  output_path          : ", output_path)
message("  intensity_columns    : ", int_start, ":", int_end,
        "  (", int_end - int_start + 1, " samples)")
message("  normalization_method : ", method)
message("  generate_MAplots     : ", gen_ma)
message("  suffix               : ", if (nzchar(suffix)) suffix else "(none)")

# ── Run ────────────────────────────────────────────────────────────────────
res <- workflow_QC(
  data_path            = data_path,
  filetype             = "tsv",
  sep                  = "\t",
  dec                  = ".",
  header               = TRUE,
  sheet                = 1,
  output_path          = output_path,
  output_type          = "tsv",
  intensity_columns    = int_start:int_end,
  normalization_method = method,
  plot_device          = "png",
  PCA_label            = TRUE,
  suffix               = suffix,
  generate_MAplots     = gen_ma
)

if (!is.null(res) && !is.null(res$message)) {
  message(paste(res$message, collapse = "\n"))
}

# ── Confirm the normalized table exists, so callers fail loudly ────────────
expected <- file.path(output_path, paste0("D_norm_ID", suffix, ".tsv"))
if (!file.exists(expected)) {
  message("ERROR: expected normalized table not found: ", expected)
  message("Files written to ", output_path, ":")
  message(paste(" ", list.files(output_path), collapse = "\n"))
  quit(status = 1)
}

message("✓ Normalized table: ", expected)
