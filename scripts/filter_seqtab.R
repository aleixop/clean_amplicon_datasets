#!/usr/bin/env Rscript

library(argparser, quietly = T)

parser <- arg_parser(
  name = "filter_seqtab.R",
  description = "Filter a dada2 seqtab file by minimum abundance, occurrence and length.",
  hide.opts = TRUE
)

parser <- add_argument(parser,
  arg = "--seqtab",
  help = "Path to the seqtab file."
)

parser <- add_argument(parser,
  arg = "--min_abundance",
  help = "ASVs with less than this number of overall reads will be removed.",
  default = 1,
  type = "integer"
)

parser <- add_argument(parser,
  arg = "--min_occurrence",
  help = "ASVs present in less than this number of samples will be removed.",
  default = 1,
  type = "integer"
)

parser <- add_argument(parser,
  arg = "--min_bp",
  help = "ASVs shorter than this number of base pairs will be removed.",
  default = 32,
  type = "integer"
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

# Read args ---------------------------------------------------------------

args <- parse_args(parser)

seqtab <- readRDS(args$seqtab)
min_abundance <- args$min_abundance
min_occurrence <- args$min_occurrence
min_bp <- args$min_bp
out_seqtab <- args$out_seqtab
out_removed_seqs <- args$out_removed_seqs

# Keep ASVs that meet the filtering criteria ------------------------------

keep_asvs <-
  colSums(seqtab) >= min_abundance &
    colSums(seqtab > 0) >= min_occurrence &
    nchar(colnames(seqtab)) >= min_bp

discard_asvs <-
  !keep_asvs

# Filter seqtab -----------------------------------------------------------

filtered_seqtab <-
  seqtab[, keep_asvs]

removed_asvs <-
  colSums(seqtab[, discard_asvs]) |>
  as_tibble(rownames = "ASV") |>
  dplyr::rename(reads = value) |>
  mutate(
    step_removed =
      case_when(
        nchar(ASV) < min_bp ~ "Filtering - Too short",
        reads < min_abundance ~ "Filtering - Low abundance",
        TRUE ~ "Filtering - Low occurrence"
      )
  )

# Write output ------------------------------------------------------------

saveRDS(filtered_seqtab, out_seqtab)
write_tsv(removed_asvs, out_removed_seqs)