library(tidyverse)
library(stringr)

# read in files
setwd("CMSC701_final")
report_dir <- "kraken2_reports"
key_file <- "key.txt"

# read in key file
key <- read_tsv(
  key_file,
  col_types = cols(.default = "c")) %>%
  select(Accession, Taxid) %>%
  rename(
    accession = Accession,
    true_taxid = Taxid)

# grab all report files
report_files <- list.files(
  report_dir,
  pattern = "\\.out$",
  full.names = TRUE)

# report processing fxn
process_report <- function(report_file, key) {
  
  db_name <- tools::file_path_sans_ext(
    basename(report_file))
  
  cat("Processing:", db_name, "\n")
  
  kraken <- read_tsv(
    report_file,
    col_names = FALSE,
    col_types = cols(.default = "c"))
  
  colnames(kraken) <- c(
    "status",
    "read_id",
    "classification",
    "read_length",
    "kmer_hits")
  
  kraken <- kraken %>%
    mutate(
      accession = str_extract(
        read_id,
        "GCF_\\d+\\.\\d+"))
  
  # add column for the true taxid
  kraken <- kraken %>%
    left_join(key, by = "accession")

  # grab predicted taxid
  kraken <- kraken %>%
    mutate(
      predicted_taxid = str_extract(
        classification,
        "(?<=taxid )\\d+"))
  
  # deal with unclassified reads
  kraken <- kraken %>%
    mutate(
      predicted_taxid = ifelse(
        status == "U",
        NA,
        predicted_taxid))
  
  # check if correct
  kraken <- kraken %>%
    mutate(
      correct = predicted_taxid == true_taxid)
  
  # Calculate TP/FP/FN
  TP <- sum(kraken$correct, na.rm = TRUE)
  
  FP <- sum(
    !is.na(kraken$predicted_taxid) &
      kraken$predicted_taxid != kraken$true_taxid,
    na.rm = TRUE)
  
  FN <- sum(
    is.na(kraken$predicted_taxid) |
      kraken$predicted_taxid != kraken$true_taxid,
    na.rm = TRUE)
  
  # calculate precision and recall
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  f1 <- 2 * precision * recall / (precision + recall)
  
  # make a summary table
  tibble(
    Database = db_name,
    TP = TP,
    FP = FP,
    FN = FN,
    Precision = precision,
    Recall = recall,
    F1 = f1
  )
}

# run everything
all_metrics <- map_dfr(
  report_files,
  process_report,
  key = key)

print(all_metrics)
write_tsv(
  all_metrics,
  "all_database_metrics.tsv")

# PLOT
plot_data <- all_metrics %>%
  pivot_longer(
    cols = c(Precision, Recall, F1),
    names_to = "Metric",
    values_to = "Score")

ggplot(
  plot_data,
  aes(
    x = Database,
    y = Score,
    fill = Metric)
  ) +
  
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
    ) +
  
  geom_text(
    aes(label = sprintf("%.2f", Score)),
    position = position_dodge(width = 0.8),
    vjust = -0.4,
    size = 4
    ) +
  
  scale_y_continuous(
    limits = c(0, 1)
    ) +
  
  scale_fill_manual(
    values = c(
      "Precision" = "steelblue",
      "Recall" = "darkorange",
      "F1" = "forestgreen")
    ) +
  
  labs(
    title = "Classification Performance",
    x = "Database",
    y = "Score",
    fill = "Metric") +
  
  theme_bw(base_size = 14) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1)
  )
