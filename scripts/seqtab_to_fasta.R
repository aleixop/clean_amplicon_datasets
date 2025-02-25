#!/usr/bin/env Rscript

library(argparser, quietly = T)

parser <- arg_parser(
  name = "seqtab_to_fasta.R",
  description = "Export a fasta file from a dada2 seqtab file.",
  hide.opts = TRUE
)

parser <- add_argument(parser,
  arg = "--seqtab",
  help = "Path to the seqtab file."
)

parser <- add_argument(parser,
  arg = "--out_fasta",
  default = "seqtab.fasta",
  help = "Path where fasta file should be written."
)

parser <- add_argument(parser,
  arg = "--add_sizes",
  flag = TRUE,
  help = "Append total size at the end of each header. This will be added with the ';size=<nreads>' suffix"
)

parser <- add_argument(parser,
  arg = "--add_names",
  flag = TRUE,
  help = "Create new names. These will be added with the form 'ASV_n' and ordered by decreasing total reads. By default, sequences will be used as headers."
)


# Read args ---------------------------------------------------------------

args <- parse_args(parser)

seqtab <- readRDS(args$seqtab)
out_fasta <- args$out_fasta
add_sizes <- args$add_sizes
add_names <- args$add_names

suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))


# Create out fasta --------------------------------------------------------

seqs <-
  tibble(
    seq = colnames(seqtab),
    size = colSums(seqtab)
  )

if (add_names) {
  seqs_name <-
    seqs |>
    arrange(-size) |>
    mutate(name = paste0("ASV_", row_number()))
} else {
  seqs_name <-
    seqs |>
    mutate(name = seq)
}

if (add_sizes) {
  fasta <-
    seqs_name |>
    mutate(name = paste0(name, ";size=", size))
} else {
  fasta <-
    seqs_name
}

fasta <-
  fasta |>
  select(name, seq) |>
  deframe() |>
  DNAStringSet()

writeXStringSet(fasta, out_fasta, width = 10000)
