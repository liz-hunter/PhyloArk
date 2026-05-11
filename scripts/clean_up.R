#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "  Rscript scripts/clean_up.R <input_list.tsv> [output_cleaned_list.tsv] [--remove-organellar] ...",
  "",
  "Options:",
  "  --remove-organellar   Remove rows flagged as organellar (instead of only warning/listing).",
  "  --quiet-organellar    Do not print the organellar warning/list (still detects/removes if enabled).",
  "  --no-report           Do not print cleanup summary / removed-row reports.",
  "",
  "Example:",
  "  Rscript scripts/clean_up.R ncbi_assemblies.tsv assemblies_cleaned.tsv --remove-organellar",
  sep = "\n"
)

if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  cat(usage, "\n")
  quit(status = 0)
}

# Split positional vs flags
positional <- args[!str_starts(args, "--")]
flags <- args[str_starts(args, "--")]

if (!(length(positional) %in% c(1, 2))) {
  cat("ERROR: Expected 1 or 2 positional arguments.\n\n", usage, "\n", sep = "")
  quit(status = 2)
}

infile <- positional[1]

# Default output: <input_basename>_cleaned.tsv (in current working directory)
if (length(positional) == 2) {
  outfile <- positional[2]
} else {
  base <- tools::file_path_sans_ext(basename(infile))
  outfile <- file.path(dirname(infile), paste0(base, "_cleaned.tsv"))
}

remove_organellar <- "--remove-organellar" %in% flags
quiet_organellar  <- "--quiet-organellar"  %in% flags
report            <- !("--no-report" %in% flags)

# Validate flags
known_flags <- c("--remove-organellar", "--quiet-organellar", "--no-report")
unknown <- setdiff(flags, known_flags)
if (length(unknown) > 0) {
  cat("ERROR: Unknown option(s): ", paste(unknown, collapse = ", "), "\n\n", usage, "\n", sep = "")
  quit(status = 2)
}

if (!file.exists(infile)) {
  cat("ERROR: Input file does not exist: ", infile, "\n\n", usage, "\n", sep = "")
  quit(status = 2)
}

# Source lib relative to this script's location
get_script_dir <- function() {
  # Works when invoked via Rscript; falls back to working directory
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  match <- grep(file_arg, cmd_args, value = TRUE)
  if (length(match) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", match[1]), winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir()
stage_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
lib_path <- file.path(stage_dir, "lib", "clean_up_lib.R")

if (!file.exists(lib_path)) {
  cat("ERROR: Could not find library file: ", lib_path, "\n", sep = "")
  quit(status = 2)
}

source(lib_path)

# Run
df_out <- clean_assembly_list(
  infile = infile,
  warn_organellar = TRUE,
  remove_organellar = remove_organellar,
  report = report,
  out = stdout()
)

# If quiet requested, suppress only the organellar listing (summary still prints unless --no-report)
if (quiet_organellar) {
}

readr::write_tsv(df_out, outfile)
cat("Wrote ", nrow(df_out), " rows to ", outfile, "\n", sep = "")
