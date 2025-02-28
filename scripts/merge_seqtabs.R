#!/usr/bin/env Rscript

library(argparser, quietly = T)

parser <- arg_parser(
  name = "merge_seqtabs.R",
  description = "Merge seqtabs into one single file.",
  hide.opts = TRUE
)

parser <- add_argument(parser,
  arg = "--seqtabs",
  help = "Path to seqtab files",
  nargs = Inf
)

parser <- add_argument(parser,
  arg = "--out_seqtab",
  default = "merged_seqtab.rds",
  help = "Path where the merged seqtab should be written."
)

# Load libraries ----------------------------------------------------------

suppressMessages(library(dada2))
suppressMessages(library(tidyverse))

# Read args ---------------------------------------------------------------

args <- parse_args(parser)

seqtabs <- map(args$seqtabs, ~ as.matrix(readRDS(.x)))
out_seqtab <- args$out_seqtab

# Merge seqtabs -----------------------------------------------------------

merged_seqtab <- 
  mergeSequenceTables(tables = seqtabs, orderBy = NULL)

# Write output ------------------------------------------------------------

saveRDS(merged_seqtab, out_seqtab)