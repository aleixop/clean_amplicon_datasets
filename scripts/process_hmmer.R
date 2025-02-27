#!/usr/bin/env Rscript

library(argparser, quietly = T)

parser <- arg_parser(
  name = "process_hmmer.R",
  description = "Remove ASVs not belonging to 18S from a seqtab based on a hmmsearch tblout table",
  hide.opts = TRUE
)

parser <- add_argument(parser,
  arg = "--seqtab",
  help = "Path to the seqtab file."
)

parser <- add_argument(parser,
  arg = "--hmmer",
  help = "Path to the hmmsearch file."
)

parser <- add_argument(parser,
  arg = "--min_evalue",
  default = 1e-5,
  help = "Minimum e-value to filter hmmsearch results."
)

parser <- add_argument(parser,
  arg = "--out_seqtab",
  default = "filtered_seqtab.rds",
  help = "Path where filtered seqtab file should be written."
)

parser <- add_argument(parser,
  arg = "--out_removed_seqs",
  default = "removed_seqs.txt",
  help = "Path where file with removed ASVs should be written."
)

# Load libraries ----------------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))


# Aux functions -----------------------------------------------------------

read_hmmer <- function(hmmer) {
  hmmer_colnames <-
    c(
      "domain_name", "domain_accession", "query_name", "query_accession", "sequence_evalue", "sequence_score",
      "sequence_bias", "best_domain_evalue", "best_domain_score", "best_domain_bis", "domain_number_exp", "domain_number_reg",
      "domain_number_clu", "domain_number_ov", "domain_number_env", "domain_number_dom", "domain_number_rep", "domain_number_inc", "description"
    )

  read_table(hmmer, col_names = hmmer_colnames, comment = "#", show_col_types = F)
}

# Read args ---------------------------------------------------------------

args <- parse_args(parser)

seqtab <- readRDS(args$seqtab)
hmmer <- read_hmmer(args$hmmer)
min_evalue <- args$min_evalue
out_seqtab <- args$out_seqtab
out_removed_seqs <- args$out_removed_seqs

# Filter hmmer ------------------------------------------------------------

hmmer_filtered <-
  hmmer |>
  filter(sequence_evalue <= min_evalue)

keep_asvs <-
  colnames(seqtab) %in% hmmer_filtered$domain_name

discard_asvs <-
  !keep_asvs

# Filter seqtab -----------------------------------------------------------

filtered_seqtab <-
  seqtab[, keep_asvs]

removed_asvs <-
  colSums(seqtab[, discard_asvs]) |>
  as_tibble(rownames = "ASV") |>
  dplyr::rename(reads = value) |>
  mutate(step_removed = "3.HMMER search 18S")

# Write output ------------------------------------------------------------

saveRDS(filtered_seqtab, out_seqtab)
write_tsv(removed_asvs, out_removed_seqs)
