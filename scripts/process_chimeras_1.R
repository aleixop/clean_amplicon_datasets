#!/usr/bin/env Rscript

library(argparser, quietly = T)

parser <- arg_parser(
  name = "process_chimeras_1.R",
  description = "Take a blast of seqs with trimmed tails against the full sequences and output putative chimeras.",
  hide.opts = TRUE
)

parser <- add_argument(parser,
  arg = "--seqtab",
  help = "Path to the seqtab file."
)

parser <- add_argument(parser,
  arg = "--blast",
  help = "Path to the blast file."
)

parser <- add_argument(parser,
  arg = "--fasta",
  help = "Path to the fasta file used for the blast search."
)

parser <- add_argument(parser,
  arg = "--out_list",
  default = "out_list.txt",
  help = "Path where ASVs to consider in the second chimera step should be written."
)


# Load libraries ----------------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))

# Aux functions -----------------------------------------------------------

read_blast <- function(blast) {
  read_tsv(blast,
    col_names = c(
      "qseqid", "sseqid", "pident", "length", "mismatch",
      "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
      "qlen", "slen"
    ),
    show_col_types = F
  ) |>
    filter(qseqid != sseqid) # remove self hits
}

# Read args ---------------------------------------------------------------

args <- parse_args(parser)

seqtab <- readRDS(args$seqtab)
blast_tab <- read_blast(args$blast)
fasta <- readDNAStringSet(args$fasta)
out_list <- args$out_list

# Process fasta -----------------------------------------------------------

fasta_df <-
  tibble(
    ASVid = names(fasta),
    ASV = as.character(fasta)
  )

# Process seqtab ----------------------------------------------------------

asvs_by_rank <-
  seqtab |>
  colSums() |>
  as_tibble(rownames = "ASV") |>
  arrange(-value) |>
  dplyr::rename(total_counts = value) |>
  mutate(rank = row_number()) |>
  left_join(fasta_df, by = "ASV")

# Process blast -----------------------------------------------------------

## Only keep comparisons with max difference between ranks. We want the query to have a low rank and subject a high rank

blast_processed <-
  blast_tab |>
  left_join(asvs_by_rank |> select(qseqid = ASVid, rank_qseqid = rank), by = "qseqid") |>
  left_join(asvs_by_rank |> select(sseqid = ASVid, rank_sseqid = rank), by = "sseqid") |>
  mutate(diff_rank = rank_qseqid / rank_sseqid) |>
  filter(diff_rank > 1) |> # force query to always have a lower rank than subject
  group_by(qseqid) |>
  slice_max(order_by = diff_rank, n = 1) # for each query, keep the comparison with highest rank difference

asv_for_step2 <- # those that will be mapped again using their tails
  blast_processed |>
  pull(qseqid)

write_lines(asv_for_step2, out_list)
