#!/usr/bin/env Rscript

library(argparser, quietly = T)

parser <- arg_parser(
  name = "final_report.R",
  description = "Create a report of removed sequences from the amplicon cleaning pipeline.",
  hide.opts = TRUE
)

parser <- add_argument(parser,
  arg = "--original_seqtab",
  help = "Path to the original seqtab file."
)

parser <- add_argument(parser,
  arg = "--final_seqtab",
  help = "Path to the final seqtab file."
)

parser <- add_argument(parser,
  arg = "--removed_seqs",
  help = "Path to the files containing sequences removed during all the cleaning process. A total of 5 files must be supplied.",
  nargs = Inf
)

parser <- add_argument(parser,
  arg = "--out_report",
  default = "final_report.txt",
  help = "Path where the final report should be written."
)

parser <- add_argument(parser,
  arg = "--out_removed_seqs",
  default = "removed_seqs.txt",
  help = "Path where file with all removed ASVs should be written."
)

# Load libraries ----------------------------------------------------------

suppressMessages(library(tidyverse))

# Read args ---------------------------------------------------------------

args <- parse_args(parser)

original_seqtab <- readRDS(args$original_seqtab)
final_seqtab <- readRDS(args$final_seqtab)
removed_seqs <- args$removed_seqs
out_report <- args$out_report
out_removed_seqs <- args$out_removed_seqs

if (length(removed_seqs) != 5) {
  len <- length(removed_seqs)
  stop(paste0("Script expects 5 report files and you provided ", len, ". Check this and come back later."))
}


# Process reports ---------------------------------------------------------

removed_seqs_df <-
  map_df(removed_seqs, ~ read_tsv(.x, col_types = cols("c", "n", "c")))

# Create final report -----------------------------------------------------

all_steps <-
  c(
    "Original seqtab",
    "1.Filtering - Too short",
    "1.Filtering - Low abundance",
    "1.Filtering - Low occurrence",
    "2.Clustering",
    "3.HMMER search 18S",
    "4.Internal gaps - Too many gaps",
    "4.Internal gaps - Inverted sequence",
    "5.Chimera",
    "Total removed",
    "Final seqtab"
  )

seqtab_stats <-
  tibble(
    Step = c("Original seqtab", "Final seqtab", "Total removed"),
    asvs = c(ncol(original_seqtab), ncol(final_seqtab), nrow(removed_seqs_df)),
    reads = c(sum(original_seqtab), sum(final_seqtab), sum(removed_seqs_df$reads))
  )

report <-
  removed_seqs_df |>
  group_by(Step = step_removed) |>
  summarise(
    `asvs` = n(),
    `reads` = sum(reads),
    .groups = "drop"
  ) |>
  bind_rows(seqtab_stats) |>
  complete(Step = all_steps, fill = list(asvs = 0, reads = 0)) |>
  mutate(Step = factor(Step, levels = all_steps)) |>
  arrange(Step)

# Write output ------------------------------------------------------------

write_tsv(removed_seqs_df, out_removed_seqs)
write_tsv(report, out_report)

# Print final report ------------------------------------------------------

nasvs <- ncol(original_seqtab) - ncol(final_seqtab)
nreads <- sum(original_seqtab) - sum(final_seqtab)
perc_nasvs <- round(100 * nasvs / ncol(original_seqtab), 2)
perc_nreads <- round(100 * nreads / sum(original_seqtab), 2)

cat(paste0("A total of ", nasvs, " ASVs (", perc_nasvs, "%) were removed. These represent ",nreads," reads (",perc_nreads,'%).\n'))
cat(paste0("Here you have a more complete report. You can also find it in '",out_report,"'."))
knitr::kable(report)
