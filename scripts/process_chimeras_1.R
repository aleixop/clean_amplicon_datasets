suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))

args <- commandArgs(trailingOnly = T)
fasta_file <- args[1]
seqtab_file <- args[2]
blast_file <- args[3]
out_list <- args[4]

# Read blast --------------------------------------------------------------

blast_tab <- 
  read_tsv(blast_file, 
           col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", 
                         "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                         "qlen", "slen")) |> 
  filter(qseqid != sseqid) # remove self hits

# Read and process fasta --------------------------------------------------

fasta <- 
  readDNAStringSet(fasta_file)

fasta_df <- 
  tibble(ASVid = names(fasta),
         ASV = as.character(fasta))

# Read and process counts -------------------------------------------------

counts <- 
  readRDS(seqtab_file)[,as.character(fasta)]

asvs_by_rank <-
  counts |> 
  colSums() |> 
  as_tibble(rownames = 'ASV') |> 
  arrange(-value) |> 
  dplyr::rename(total_counts = value) |> 
  mutate(rank = row_number()) |> 
  left_join(fasta_df)
 
# Process blast -----------------------------------------------------------

## Only keep comparisons with max difference between ranks. We want the query to have a low rank and subject a high rank

blast_processed <- 
  blast_tab |> 
  left_join(asvs_by_rank |> select(qseqid = ASVid, rank_qseqid = rank)) |> 
  left_join(asvs_by_rank |> select(sseqid = ASVid, rank_sseqid = rank)) |> 
  mutate(diff_rank = rank_qseqid/rank_sseqid) |> 
  filter(diff_rank > 1) |> # force query to always have a lower rank than subject
  group_by(qseqid) |> 
  slice_max(order_by = diff_rank, n = 1) # for each query, keep the comparison with highest rank difference

asv_for_step2 <- # those that will be mapped again using their tails
  blast_processed |> 
  select(qseqid)

write_tsv(asv_for_step2, out_list, col_names = F)