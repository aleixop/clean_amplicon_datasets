suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))
suppressMessages(library(dada2))

args <- commandArgs(trailingOnly = T)
fasta_file <- args[1]
seqtab_file <- args[2]
blast_file1 <- args[3]
blast_file2 <- args[4]
out_fasta <- args[5]
out_seqtab <- args[6]

# Read blast --------------------------------------------------------------

blast_tab <- 
  read_tsv(blast_file1, 
           col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", 
                         "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                         "qlen", "slen")) |> 
  filter(qseqid != sseqid) # remove self hits

blast_tab2 <- 
  read_tsv(blast_file2, 
           col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", 
                         "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                         "qlen", "slen")) |> 
  filter(qseqid != sseqid)

# Read and process fasta --------------------------------------------------

fasta <- 
  readDNAStringSet(fasta_file)

fasta_df <- 
  tibble(ASVid = names(fasta),
         ASV = as.character(fasta))

# Read and process counts -------------------------------------------------

counts <- 
  readRDS(seqtab_file)[,as.character(fasta)]

occurrence <- 
  colSums(counts > 0) |> 
  as_tibble(rownames = 'ASV') |> 
  dplyr::rename(occurrence = value) |> 
  mutate(occurrence_perc = 100*occurrence/nrow(counts))

asvs_by_rank <-
  counts |> 
  colSums() |> 
  as_tibble(rownames = 'ASV') |> 
  arrange(-value) |> 
  dplyr::rename(total_counts = value) |> 
  mutate(rank = row_number()) |> 
  left_join(fasta_df) |> 
  left_join(occurrence)
 
# Process blast files -----------------------------------------------------

blast_processed <- 
  blast_tab |> 
  left_join(asvs_by_rank |> select(qseqid = ASVid, rank_qseqid = rank)) |> 
  left_join(asvs_by_rank |> select(sseqid = ASVid, rank_sseqid = rank)) |> 
  mutate(diff_rank = rank_qseqid/rank_sseqid) |> 
  filter(diff_rank > 1) |> # force query to always have a lower rank than subject
  group_by(qseqid) |> 
  slice_max(order_by = diff_rank, n = 1) # for each query, keep the comparison with highest rank difference

blast2_processed <- 
  blast_tab2 |> 
  left_join(asvs_by_rank |> select(qseqid = ASVid, rank_qseqid = rank)) |> 
  left_join(asvs_by_rank |> select(sseqid = ASVid, rank_sseqid = rank)) |> 
  mutate(diff_rank = rank_qseqid/rank_sseqid) |> 
  filter(diff_rank > 1) |> # force query to always have a lower rank than subject
  group_by(qseqid) |> 
  filter(pident == max(pident)) |> # keep highest pident per query
  slice_max(order_by = diff_rank, n = 1) # for each query, keep the comparison with highest rank difference

# Parents file raw --------------------------------------------------------

blast_step1 <- 
  blast_processed |> 
  dplyr::rename_with(~ paste0(.x,'_1'), -qseqid)
  
blast_step2 <- 
  blast2_processed |> 
  dplyr::rename_with(~ paste0(.x,'_2'), -qseqid)

parents_file_0 <- 
  blast_step1 |> 
  filter(qseqid %in% blast_step2$qseqid) |>
  arrange(rank_qseqid_1) |> 
  left_join(blast_step2) |> 
  select(-contains('rank'))

# Compare parents ---------------------------------------------------------

pident_2seqs <- function(query, subject){
  
  require(dada2)
  
  length_query = nchar(query)
  alignment = nwhamming(query, subject)
  pident = 100*(length_query-alignment)/length_query
  
  return(pident)
}

parents_comparison <- 
  parents_file_0 |> 
  select(qseqid, sseqid_1, sseqid_2) |> 
  left_join(fasta_df |> select(qseqid = ASVid, seq_qseqid = ASV)) |> 
  left_join(fasta_df |> select(sseqid_1 = ASVid, seq_1 = ASV)) |> 
  left_join(fasta_df |> select(sseqid_2 = ASVid, seq_2 = ASV)) |> 
  mutate(qseqid_seq1 = pident_2seqs(seq_qseqid, seq_1),
         qseqid_seq2 = pident_2seqs(seq_qseqid, seq_2),
         seq1_seq2 = pident_2seqs(seq_1, seq_2))

# Final parents file ------------------------------------------------------

parents_file <- 
  parents_file_0 |> 
  filter(sseqid_1 != sseqid_2) |> 
  left_join(parents_comparison) |> 
  select(qseqid, sseqid_1, sseqid_2, qseqid_seq1, qseqid_seq2, seq1_seq2, everything())

# Detect chimeras ---------------------------------------------------------

max_pident_query_parent = 99
max_pident_parents = 95
max_occurrence = 5

chimeras <- 
  parents_file |> 
  filter(if_all(c(qseqid_seq1, qseqid_seq2), ~ .x < max_pident_query_parent),
         seq1_seq2 < max_pident_parents) |> 
  left_join(asvs_by_rank |> select(qseqid = ASVid, total_counts, occurrence_perc)) |>
  filter(occurrence_perc < max_occurrence) |> 
  select(qseqid, total_counts, occurrence_perc, everything())

out_dir <- str_match(out_fasta, '(.*)\\/')[,2]
if (is.na(out_dir)){ out_dir <- './'}

write_tsv(chimeras, paste0(out_dir,'chimeras_parents_file.tsv'))

# Filter fasta ------------------------------------------------------------

chimeras_list <-
  chimeras |> 
  pull(qseqid)

fasta_filtered <- 
  fasta[!names(fasta) %in% chimeras_list]

writeXStringSet(fasta_filtered, out_fasta)

# Filter counts -----------------------------------------------------------

counts_filtered <- 
  counts[,as.character(fasta_filtered)]

saveRDS(counts_filtered, out_seqtab)

# Report ------------------------------------------------------------------

perc_asvs_removed <- paste0(' (',round(100*(length(fasta)-length(fasta_filtered))/length(fasta),2),'%)')
perc_reads_removed <- paste0(' (',round(100*(sum(counts)-sum(counts_filtered))/sum(counts),3),'%)')

cat('\n### Remove chimeric ASVs ###')
cat(paste0('\n# ASVs removed: ',length(fasta)-length(fasta_filtered), perc_asvs_removed))
cat(paste0('\n# Reads removed: ',sum(counts)-sum(counts_filtered), perc_reads_removed))
cat(paste0('\n# ASVs kept: ', length(fasta_filtered)))
cat(paste0('\n# Reads kept: ', sum(counts_filtered)))

chimeras_stats <- 
  asvs_by_rank |> 
  filter(ASVid %in% chimeras_list)

cat('\n# Summary of occurrence (%) of removed ASVS:\n')
summary(chimeras_stats$occurrence_perc)
cat('# Summary of total reads of removed ASVS:\n')
summary(chimeras_stats$total_counts)
cat('\n')