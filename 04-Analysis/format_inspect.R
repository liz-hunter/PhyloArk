library(tidyverse)

# read in files
setwd("/Users/lizhunter/Documents/school/classes/CMSC701_final")
inspect_file <- "results/inspect/refseq_inspect.txt"
key_file <- "key.txt"

# import the key file
key <- read_tsv(
  key_file,
  col_types = cols(.default = "c")
) %>%
  select(Accession, Taxid) %>%
  rename(
    accession = Accession,
    taxid = Taxid
  )

# read in the kraken inspect file
inspect <- read_tsv(
  inspect_file,
  col_names = FALSE,
  col_types = cols(.default = "c")
)

# format the kraken file
colnames(inspect) <- c(
  "perc_comp",
  "min_incl",
  "min_excl",
  "level",
  "taxid",
  "name"
)

