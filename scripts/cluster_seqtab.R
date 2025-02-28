#!/usr/bin/env Rscript

library(argparser, quietly = T)

parser <- arg_parser(
  name = "cluster_seqtab.R",
  description = "Cluster a dada2 seqtab file base on percentage identity and minimum coverage.",
  hide.opts = TRUE
)

parser <- add_argument(parser,
  arg = "--seqtab",
  help = "Path to the seqtab file."
)

parser <- add_argument(parser,
  arg = "--perc_identity",
  help = "Clustering identity in a 0-100 scale.",
  default = 100,
  type = "integer"
)

parser <- add_argument(parser,
  arg = "--min_coverage",
  help = "Minimum fraction of sequence positions that must be overlapping for a sequence to be clustered with the cluster representative.",
  default = 0.9,
  type = "integer"
)

parser <- add_argument(parser,
  arg = "--representative_method",
  help = "Choose whether representatives should be chosen based on 'abundance' (most abundant ASV as representative) or 'length' (longest ASV as representative).",
  default = "abundance"
)

parser <- add_argument(parser,
  arg = "--out_seqtab",
  default = "clustered_seqtab.rds",
  help = "Path where filtered seqtab file should be written."
)

parser <- add_argument(parser,
  arg = "--out_clusters",
  default = "clusters.tsv",
  help = "Path where correspondence table between ASVs and cluster representatives should be written."
)

parser <- add_argument(parser,
  arg = "--out_removed_seqs",
  default = "removed_seqs.txt",
  help = "Path where file with removed ASVs should be written."
)

# Read args ---------------------------------------------------------------

args <- parse_args(parser)

seqtab <- readRDS(args$seqtab)
perc_identity <- args$perc_identity
representative_method <- args$representative_method
min_coverage <- args$min_coverage
out_seqtab <- args$out_seqtab
out_clusters <- args$out_clusters
out_removed_seqs <- args$out_removed_seqs

suppressMessages(library(tidyverse))
suppressMessages(library(DECIPHER))
suppressMessages(library(Biostrings))
suppressMessages(library(data.table))

# Create ASV df -----------------------------------------------------------

asv_df <-
  tibble(
    seq = colnames(seqtab),
    size = colSums(seqtab)
  ) |>
  mutate(
    length = nchar(seq),
    name = paste0("ASV_", row_number())
  ) # add names to avoid duplicates


# Create DNAStringSet for clustering --------------------------------------

dna <-
  asv_df |>
  select(name, seq) |>
  deframe() |>
  DNAStringSet()


# Find clusters of ASVs ---------------------------------------------------

cutoff <- (100 - perc_identity) / 100

clusters <-
  Clusterize(
    myXStringSet = dna,
    cutoff = cutoff,
    minCoverage = min_coverage,
    processors = NULL
  ) |>
  as_tibble(rownames = "name")

asv_df_clusters <-
  asv_df |>
  left_join(clusters, by = "name")

# Choose representatives --------------------------------------------------

if (representative_method == "abundance") {
  clusters_out <-
    asv_df_clusters |>
    group_by(cluster) |>
    mutate(representative = seq[size == max(size)][1]) |>
    select(-name)
} else if (representative_method == "length") {
  clusters_out <-
    asv_df_clusters |>
    group_by(cluster) |>
    mutate(representative = seq[length == max(length)][1]) |>
    select(-name)
} else {
  stop("Representative method selection must be one of the following: 'abundance', 'length'")
}

representatives <-
  clusters_out |>
  ungroup() |>
  select(seq, representative)

export_clusters <-
  clusters_out |>
  group_by(cluster) |>
  filter(n() > 1) |>
  arrange(cluster, -size) |>
  mutate(type = case_when(
    seq != representative ~ "hit",
    TRUE ~ "representative"
  ))

cat(paste0("A total of ", sum(export_clusters$type == "hit"), " sequences were clustered.\n"))

write_tsv(export_clusters, out_clusters)

# Export list of "removed" seqs -------------------------------------------

removed_asvs <- # not really removed, but clustered, that is why reads = 0, this is needed for final reporting
  clusters_out |>
  ungroup(cluster) |>
  filter(seq != representative) |>
  select(ASV = seq) |>
  mutate(
    reads = 0,
    step_removed = "Clustering"
  )

write_tsv(removed_asvs, out_removed_seqs)

# Create clustered seqtab -------------------------------------------------

seqtab_with_representatives <-
  seqtab |>
  as.data.table(keep.rownames = "sample") |> # transform to data.table
  melt(
    id.vars = "sample",
    variable.name = "seq",
    value.name = "abundance"
  ) |> # put in long format
  merge(representatives, by = "seq") # join representatives

merged_seqtab <-
  seqtab_with_representatives[, .(abundance = sum(abundance)), keyby = .(representative, sample)] |> # collapse countsby representatives
  dcast(sample ~ representative,
    value.var = "abundance",
    fill = 0
  ) |> # put in wide format
  as_tibble() |> # transform to table
  column_to_rownames("sample") |> # put sample as rownames to keep original format
  setcolorder(unique(representatives$representative)) # keep original ASV order as much as possible

saveRDS(merged_seqtab, out_seqtab)
