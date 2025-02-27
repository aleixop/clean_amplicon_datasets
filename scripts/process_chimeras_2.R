#!/usr/bin/env Rscript

library(argparser, quietly = T)

parser <- arg_parser(
  name = "process_chimeras_2.R",
  description = "Take 2 blast of seqs with (1) trimmed tails and (2) trimmed heads against the full sequences. Detect chimeras and filter seqtab.",
  hide.opts = TRUE
)

parser <- add_argument(parser,
  arg = "--seqtab",
  help = "Path to the seqtab file."
)

parser <- add_argument(parser,
  arg = "--blast1",
  help = "Path to the blast file with trimmed tails."
)

parser <- add_argument(parser,
  arg = "--blast2",
  help = "Path to the blast file with trimmed heads."
)

parser <- add_argument(parser,
  arg = "--fasta",
  help = "Path to the fasta file used for the blast search."
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

parser <- add_argument(parser,
  arg = "--out_chimeras",
  default = "chimeras.txt",
  help = "Path where chimeras report should be written."
)

parser <- add_argument(parser,
  arg = "--max_pident_chimera_parent",
  default = 99,
  help = "Maximum percentage of sequence identity between a chimera and both of their parents."
)

parser <- add_argument(parser,
  arg = "--max_pident_parents",
  default = 95,
  help = "Maximum percentage of sequence identity between the parents of a chimera."
)

parser <- add_argument(parser,
  arg = "--max_chimera_occurrence",
  default = 5,
  help = "Maximum occurrence (in percentage) that a chimera can have."
)

# Load libraries ----------------------------------------------------------

suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))
suppressMessages(library(dada2))

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

pident_2seqs <- function(query, subject) {
  length_query <- nchar(query)
  alignment <- nwhamming(query, subject)
  pident <- 100 * (length_query - alignment) / length_query

  return(pident)
}

# Read args ---------------------------------------------------------------

args <- parse_args(parser)

seqtab <- readRDS(args$seqtab)
blast_tab <- read_blast(args$blast1)
blast_tab2 <- read_blast(args$blast2)
fasta <- readDNAStringSet(args$fasta)
out_seqtab <- args$out_seqtab
out_removed_seqs <- args$out_removed_seqs
out_chimeras <- args$out_chimeras
max_pident_query_parent <- args$max_pident_chimera_parent
max_pident_parents <- args$max_pident_parents
max_occurrence <- args$max_chimera_occurrence

# Process fasta -----------------------------------------------------------

fasta_df <-
  tibble(
    ASVid = names(fasta),
    ASV = as.character(fasta)
  )

# Process seqtab ----------------------------------------------------------

occurrence <-
  colSums(seqtab > 0) |>
  as_tibble(rownames = "ASV") |>
  dplyr::rename(occurrence = value) |>
  mutate(occurrence_perc = 100 * occurrence / nrow(seqtab))

asvs_by_rank <-
  seqtab |>
  colSums() |>
  as_tibble(rownames = "ASV") |>
  arrange(-value) |>
  dplyr::rename(total_counts = value) |>
  mutate(rank = row_number()) |>
  left_join(fasta_df, by = "ASV") |>
  left_join(occurrence, by = "ASV")

# Process blast files -----------------------------------------------------

blast_processed <-
  blast_tab |>
  left_join(asvs_by_rank |> select(qseqid = ASVid, rank_qseqid = rank), by = "qseqid") |>
  left_join(asvs_by_rank |> select(sseqid = ASVid, rank_sseqid = rank), by = "sseqid") |>
  mutate(diff_rank = rank_qseqid / rank_sseqid) |>
  filter(diff_rank > 1) |> # force query to always have a lower rank than subject
  group_by(qseqid) |>
  slice_max(order_by = diff_rank, n = 1) # for each query, keep the comparison with highest rank difference

blast2_processed <-
  blast_tab2 |>
  left_join(asvs_by_rank |> select(qseqid = ASVid, rank_qseqid = rank), by = "qseqid") |>
  left_join(asvs_by_rank |> select(sseqid = ASVid, rank_sseqid = rank), by = "sseqid") |>
  mutate(diff_rank = rank_qseqid / rank_sseqid) |>
  filter(diff_rank > 1) |> # force query to always have a lower rank than subject
  group_by(qseqid) |>
  filter(pident == max(pident)) |> # keep highest pident per query
  slice_max(order_by = diff_rank, n = 1) # for each query, keep the comparison with highest rank difference

# Parents file raw --------------------------------------------------------

blast_step1 <-
  blast_processed |>
  dplyr::rename_with(~ paste0(.x, "_1"), -qseqid)

blast_step2 <-
  blast2_processed |>
  dplyr::rename_with(~ paste0(.x, "_2"), -qseqid)

parents_file_0 <-
  blast_step1 |>
  filter(qseqid %in% blast_step2$qseqid) |>
  arrange(rank_qseqid_1) |>
  left_join(blast_step2, by = "qseqid") |>
  select(-contains("rank"))

# Compare parents ---------------------------------------------------------

parents_comparison <-
  parents_file_0 |>
  select(qseqid, sseqid_1, sseqid_2) |>
  left_join(fasta_df |> select(qseqid = ASVid, seq_qseqid = ASV), by = "qseqid") |>
  left_join(fasta_df |> select(sseqid_1 = ASVid, seq_1 = ASV), by = "sseqid_1") |>
  left_join(fasta_df |> select(sseqid_2 = ASVid, seq_2 = ASV), by = "sseqid_2") |>
  mutate(
    qseqid_seq1 = pident_2seqs(seq_qseqid, seq_1),
    qseqid_seq2 = pident_2seqs(seq_qseqid, seq_2),
    seq1_seq2 = pident_2seqs(seq_1, seq_2)
  )

# Final parents file ------------------------------------------------------

parents_file <-
  parents_file_0 |>
  filter(sseqid_1 != sseqid_2) |> # parents should be different seqs
  left_join(parents_comparison, by = c("qseqid", "sseqid_1", "sseqid_2")) |>
  select(qseqid, sseqid_1, sseqid_2, qseqid_seq1, qseqid_seq2, seq1_seq2, starts_with("seq_"))

# Detect chimeras ---------------------------------------------------------

chimeras <-
  parents_file |>
  filter(
    if_all(c(qseqid_seq1, qseqid_seq2), ~ .x < max_pident_query_parent), # chimera should be different from parents
    seq1_seq2 < max_pident_parents
  ) |> # parents should be relatively distant
  left_join(asvs_by_rank |> select(qseqid = ASVid, total_counts, occurrence_perc), by = "qseqid") |>
  filter(occurrence_perc < max_occurrence) |> # chimera should have low occurrence
  select(qseqid, total_counts, occurrence_perc, everything())

# Filter fasta ------------------------------------------------------------

chimeras_list <-
  chimeras |>
  pull(qseqid)

fasta_filtered <-
  fasta[!names(fasta) %in% chimeras_list]

keep_asvs <-
  colnames(seqtab) %in% as.character(fasta_filtered)

discard_asvs <-
  !keep_asvs

# Filter seqtab -----------------------------------------------------------

filtered_seqtab <-
  seqtab[, keep_asvs]

removed_asvs <-
  colSums(seqtab[, discard_asvs]) |>
  as_tibble(rownames = "ASV") |>
  dplyr::rename(reads = value) |>
  mutate(step_removed = "5.Chimera")

# Write output ------------------------------------------------------------

saveRDS(filtered_seqtab, out_seqtab)
write_tsv(removed_asvs, out_removed_seqs)
write_tsv(chimeras, out_chimeras)