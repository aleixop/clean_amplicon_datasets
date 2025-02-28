#!/usr/bin/env Rscript

library(argparser, quietly = T)

parser <- arg_parser(
  name = "process_blast_internal_gaps.R",
  description = "Remove ASVs with internal gaps from a seqtab based on a blast search",
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
  arg = "--min_query_coverage",
  default = 0.99,
  help = "Minimum fraction of the query sequence that must be aligned for removal consideration."
)

parser <- add_argument(parser,
  arg = "--min_gaps_subject",
  default = 15,
  help = "Minimum length of internal gaps in the subject sequence for query removal."
)

parser <- add_argument(parser,
  arg = "--diff_ranks_log",
  default = -0.4,
  help = "Maximum difference between log-transformed ranks [log(rank_query) - log(rank_subject)] to remove the subject instead of the query. This ensures that query instead of subject is kept when the abundance rank of query is considerably lower (i.e., query is more abundant than subject)."
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
suppressMessages(library(IRanges))

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

calculate_aligned_positions <- function(start, end) {
  IRanges(start, end) |>
    reduce() |>
    width() |>
    sum()
}

print_aligned_positions <- function(start, end) {
  iranges_object <-
    IRanges(start, end) |> IRanges::reduce()

  positions <-
    paste0(start(iranges_object), "-", end(iranges_object), collapse = ",")

  return(positions)
}

# Read args ---------------------------------------------------------------

args <- parse_args(parser)

seqtab <- readRDS(args$seqtab)
blast_tab <- read_blast(args$blast)
fasta <- readDNAStringSet(args$fasta)
out_seqtab <- args$out_seqtab
out_removed_seqs <- args$out_removed_seqs
query_covered_threshold <- args$min_query_coverage
min_gaps_subject <- args$min_gaps_subject
diff_ranks_log <- args$diff_ranks_log

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

# Detect sequences with inverted fragments --------------------------------
## these need to be removed directly as they may cause problems with IRanges

inverted_signals <-
  blast_tab |>
  filter(sstart - send > 0) |>
  rowwise() |>
  mutate(pair = paste(sort(c(qseqid, sseqid)), collapse = "-")) |>
  group_by(pair) |>
  slice_head(n = 1)

inverted_seqs <-
  tibble(remove = c(inverted_signals$qseqid, inverted_signals$sseqid)) |>
  group_by(remove) |>
  tally() |>
  filter(n > 1) |>
  select(remove)

# Filter blast to remove query-subject pairs that only appear once --------

blast_tab_pairs <-
  blast_tab |>
  group_by(qseqid, sseqid) |>
  filter(n() > 1) |>
  left_join(asvs_by_rank |> select(qseqid = ASVid, rank_qseqid = rank), by = "qseqid") |>
  left_join(asvs_by_rank |> select(sseqid = ASVid, rank_sseqid = rank), by = "sseqid") |>
  filter(
    !qseqid %in% inverted_seqs$remove,
    !sseqid %in% inverted_seqs$remove
  )

# Is the query covered? ---------------------------------------------------

querys_coverage <-
  blast_tab_pairs |>
  group_by(qseqid) |>
  filter(rank_sseqid == min(rank_sseqid)) |> # take only the comparison against the most abundant ASV from all pairs
  group_by(qseqid, sseqid, qlen, slen, rank_qseqid, rank_sseqid) |>
  summarise(
    total_aln_length_query = calculate_aligned_positions(qstart, qend),
    total_aln_length_subject = calculate_aligned_positions(sstart, send),
    aligned_positions = print_aligned_positions(qstart, qend),
    .groups = "drop_last"
  ) |>
  mutate(
    pair = paste(sort(c(qseqid, sseqid)), collapse = "-"),
    gap_length_subject = slen - total_aln_length_subject
  ) |>
  group_by(pair) |>
  slice_min(order_by = qlen, with_ties = F)

covered_querys <-
  querys_coverage |>
  filter(
    total_aln_length_query >= qlen * query_covered_threshold, # most of the query should be covered
    gap_length_subject >= min_gaps_subject # a portion of the subject should be uncovered
  ) |>
  mutate(remove = case_when(
    log(rank_qseqid) - log(rank_sseqid) <= diff_ranks_log ~ sseqid,
    TRUE ~ qseqid
  ))

# Sequences to remove -----------------------------------------------------

seqs_with_gaps <- # to remove
  covered_querys |>
  ungroup() |>
  select(remove) |>
  unique() |>
  bind_rows(inverted_seqs) |>
  pull(remove)

inverted_seqs_sequences <-
  fasta[names(fasta) %in% inverted_seqs$remove] |>
  as.character()

# Filter fasta ------------------------------------------------------------

fasta_filtered <-
  fasta[!names(fasta) %in% seqs_with_gaps]

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
  mutate(step_removed = case_when(
    ASV %in% inverted_seqs_sequences ~ "3.Internal gaps - Inverted sequence",
    TRUE ~ "3.Internal gaps - Too many gaps"
  ))

# Write output ------------------------------------------------------------

saveRDS(filtered_seqtab, out_seqtab)
write_tsv(removed_asvs, out_removed_seqs)
