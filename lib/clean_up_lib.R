suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(janitor)
})

required_columns <- c(
  "assembly_name",
  "assembly_accession",
  "assembly_paired_assembly_accession",
  "organism_name",
  "organism_taxonomic_id"
)

# Reads in genome metadata table from NCBI genome browser & cleans
read_and_validate <- function(path, needed = required_columns) {
  if (!file.exists(path)) {
    stop("Input file does not exist: ", path, call. = FALSE)
  }

  df <- read_tsv(path, show_col_types = FALSE) %>%
    janitor::clean_names()

  missing <- setdiff(needed, names(df))
  if (length(missing) > 0) {
    stop("Missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  df %>% select(all_of(needed))
}

# Initial processing to find downstream issues
find_issues <- function(df) {
  out <- df %>%
    mutate(
      assembly_paired_assembly_accession = na_if(str_trim(assembly_paired_assembly_accession), ""),
      key_source = coalesce(assembly_paired_assembly_accession, assembly_accession),
      pair_key = key_source %>%
        str_remove("^GC[AF]_") %>%
        str_remove("\\..*$"),
      is_refseq = str_starts(assembly_accession, "GCF_"),
      is_organellar = str_detect(
        str_to_lower(assembly_name),
        "mitochondria|plastid|chloroplast|apicoplast"
      )
    )

  if (any(is.na(out$pair_key)) || any(out$pair_key == "")) {
    bad <- out %>% filter(is.na(pair_key) | pair_key == "")
    stop(
      "pair_key contains NA/empty values (n=",
      nrow(bad),
      "). Check assembly_accession / assembly_paired_assembly_accession.\n",
      "Example rows:\n",
      paste(utils::capture.output(print(utils::head(bad, 10))), collapse = "\n"),
      call. = FALSE
    )
  }

  out
}

# Flags organellar genome entries for potential removal
find_organellar <- function(df_annotated, out = stdout(), verbose = TRUE) {
  org_hits <- df_annotated %>% filter(is_organellar)

  if (verbose && nrow(org_hits) > 0) {
    cat(
      "\nWARNING: Potential organellar assemblies detected based on assembly_name.\n",
      "         (matching: mitochondria|plastid|chloroplast|apicoplast)\n",
      "Count: ", nrow(org_hits), "\n\n",
      sep = ""
    )

    org_hits %>%
      transmute(
        assembly_accession,
        organism_name,
        organism_taxonomic_id,
        assembly_name
      ) %>%
      write_tsv(out)

    cat("\n")
  }

  invisible(org_hits)
}

# Remove organellar assemblies from the table (optional)
remove_organellar_genomes <- function(df_annotated) {
  removed <- df_annotated %>% filter(is_organellar)
  kept    <- df_annotated %>% filter(!is_organellar)
  list(kept = kept, removed = removed)
}

# Remove duplicate assemblies from the table, preferring RefSeq over Genbank when both exist
remove_duplicates <- function(df_annotated) {
  picked <- df_annotated %>%
    group_by(pair_key) %>%
    arrange(desc(is_refseq), assembly_accession) %>%  # stable tie-break
    mutate(.rank = row_number()) %>%
    ungroup()

  kept <- picked %>% filter(.rank == 1)
  removed <- picked %>% filter(.rank > 1)

  list(kept = kept, removed = removed)
}

# Generate a summary of what was removed 
cleanup_summary <- function(df_input,
                            kept_final,
                            removed_duplicates = NULL,
                            removed_organellar = NULL,
                            out = stdout()) {
  n_input <- nrow(df_input)
  n_org   <- if (is.null(removed_organellar)) 0L else nrow(removed_organellar)
  n_dup   <- if (is.null(removed_duplicates)) 0L else nrow(removed_duplicates)
  n_out   <- nrow(kept_final)

  cat(
    "\nClean-up summary:\n",
    "  Input rows:               ", n_input, "\n",
    "  Organellar rows removed:  ", n_org, "\n",
    "  Duplicate rows removed:   ", n_dup, "\n",
    "  Output rows:              ", n_out, "\n\n",
    sep = ""
  )

  if (!is.null(removed_organellar) && n_org > 0) {
    cat("Organellar rows removed:\n")
    removed_organellar %>%
      transmute(
        removed_accession = assembly_accession,
        organism_name,
        organism_taxonomic_id,
        assembly_name
      ) %>%
      write_tsv(out)
    cat("\n")
  }

  if (!is.null(removed_duplicates) && n_dup > 0) {
    cat("Duplicate rows removed (showing kept accession for same pair_key):\n")
    removed_report <- removed_duplicates %>%
      select(pair_key,
             removed_accession = assembly_accession,
             name = organism_name) %>%
      left_join(
        kept_final %>% select(pair_key, kept_accession = assembly_accession),
        by = "pair_key"
      ) %>%
      select(-pair_key)

    write_tsv(removed_report, out)
    cat("\n")
  }

  invisible(NULL)
}

# End-to-end cleanup, outputs a nice formatted final table
clean_assembly_list <- function(infile,
                                warn_organellar = TRUE,
                                remove_organellar = FALSE,
                                report = TRUE,
                                out = stdout()) {
  df <- read_and_validate(infile)
  df_annotated <- find_issues(df)

  # optional: warn/list organellar hits (without removing)
  if (warn_organellar) {
    find_organellar(df_annotated, out = out)
  }

  # optional: remove organellar before de-dup
  removed_organellar <- NULL
  df_for_dedup <- df_annotated
  if (remove_organellar) {
    org_step <- remove_organellar_genomes(df_annotated)
    df_for_dedup <- org_step$kept
    removed_organellar <- org_step$removed
  }

  # de-dup genbank/refseq pairs
  dedup <- remove_duplicates(df_for_dedup)

  # report what happened
  if (report) {
    cleanup_summary(
      df_input = df_annotated,
      kept_final = dedup$kept,
      removed_duplicates = dedup$removed,
      removed_organellar = removed_organellar,
      out = out
    )
  }

  # final output for downstream steps
  dedup$kept %>%
    transmute(
      accession = assembly_accession,
      name      = organism_name,
      taxid     = organism_taxonomic_id,
      filename  = paste0(assembly_accession, "_", assembly_name, "_genomic.fna")
    )

}
