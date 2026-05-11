library(tidyverse)

# Load genbank report
clean_report <- read_tsv("/Users/elizabeth.hunter/Library/CloudStorage/OneDrive-FDA/Documents/PhD_coursework/CMSC-701_Algorithms/enterobacter_clean.tsv")
full_report <- read_tsv("/Users/elizabeth.hunter/Library/CloudStorage/OneDrive-FDA/Documents/PhD_coursework/CMSC-701_Algorithms/enterobacter.tsv")

# Fix headers
colnames(full_report) <- c("assembly_name", "accession", "paired_accession", "name", "taxid", "seq_length", "contig_N50", "scaffold_N50", "GC", "checkM_complete", "checkM_contam")

# Merge
merged_report <- full_report %>%
  right_join(clean_report, by = c("accession", "taxid", "name"))

# Replace NA checkM contam with 0 
merged_report <- merged_report %>%
  mutate(checkM_contam = if_else(is.na(checkM_contam), 0, checkM_contam))

# Filter by quality for the best assembly per taxid
best_reps <- merged_report %>%
  mutate(checkm_score = checkM_complete - 5 * checkM_contam) %>%
  group_by(taxid) %>%
  arrange(
    desc(checkm_score),
    desc(checkM_complete),
    checkM_contam,
    desc(scaffold_N50),
    desc(seq_length),
    accession
  ) %>%
  slice(1) %>%
  ungroup()
