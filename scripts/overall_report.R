#!/usr/bin/env Rscript

library(argparser, quietly = T)

parser <- arg_parser(
  name = "overall_report.R",
  description = "Create an overall report of initial and final reads/asvs per dataset.",
  hide.opts = TRUE
)

parser <- add_argument(parser,
  arg = "--merged_seqtab",
  help = "Path to the merged seqtab file."
)

parser <- add_argument(parser,
  arg = "--clust_seqtab",
  help = "Path to the clustered merged seqtab file."
)

parser <- add_argument(parser,
  arg = "--reports",
  help = "Paths to the report file of each dataset.",
  nargs = Inf
)

parser <- add_argument(parser,
  arg = "--out_report",
  default = "overall_report.txt",
  help = "Path where the final report should be written."
)

# Load libraries ----------------------------------------------------------

suppressMessages(library(tidyverse))

# Read args ---------------------------------------------------------------

args <- parse_args(parser)

merged_seqtab <- readRDS(args$merged_seqtab)
clust_seqtab <- readRDS(args$clust_seqtab)
reports <- args$reports
out_report <- args$out_report


# Stats for merged tables -------------------------------------------------

merged_stats <-
  tibble(
    dataset = rep("final_merged", 2),
    step = c("original", "final"),
    asvs = c(ncol(merged_seqtab), ncol(clust_seqtab)),
    reads = c(sum(merged_seqtab), sum(clust_seqtab))
  )

# Process reports ---------------------------------------------------------

report <-
  map_df(
    reports,
    ~ read_tsv(.x,
      show_col_types = F
    ) |>
      filter(str_detect(step, "seqtab")) |>
      mutate(
        dataset = str_remove(basename(.x), "_report.tsv$"),
        step = str_to_lower(str_remove(step, " seqtab$"))
      )
  ) |> 
  bind_rows(merged_stats) |> 
  pivot_wider(names_from = step,
              values_from = c(asvs, reads))

# Write output ------------------------------------------------------------

write_tsv(report, out_report)