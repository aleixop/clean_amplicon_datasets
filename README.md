# Snakemake workflow to clean amplicon datasets

This is a pipeline to clean short eukaryotic amplicon datasets processed with dada2. It has been fine-tuned for **sequences of the V4 region of the 18S ribosomal RNA gene**.

## Summary

- [Overview](#overview)
- [How to run this workflow](#how-to-run-this-workflow)
  - [Step 1: Clone this repository](#step-1-clone-this-repository)
  - [Step 2: Install required software](#step-2-install-required-software)
    - [Manual installation](#manual-installation)
    - [Conda installation](#conda-installation)
  - [Step 3: Prepare your input files](#step-3-prepare-your-input-files)
  - [Step 4: Modify parameters](#step-4-modify-parameters)
  - [Step 5: Run the pipeline](#step-5-run-the-pipeline)
    - [With manual installation](#with-manual-installation)
    - [With conda installation](#with-conda-installation)
- [Pipeline steps in detail](#pipeline-steps-in-detail)
  - [Filtering](#filtering)
  - [Clustering](#clustering)
  - [HMM search to remove non-ribosomal sequences](#hmm-search-to-remove-non-ribosomal-sequences)
  - [Remove internal gaps](#remove-internal-gaps)
  - [Chimera removal](#chimera-removal)
  - [Merge datasets and clustering](#merge-datasets-and-clustering)
- [Output explained](#output-explained)
- [Additional help](#additional-help)

## Overview

The cleaning process is based in five main steps and an additional one when processing more than one dataset:

1. **Filtering**: Removes sequences based on minimum abundance, occurrence, and length thresholds.
2. **Clustering**: Clusters identical sequences (if a 100% identity is used) that may have longer ends or terminal deletions.
3. **HMM search**: Removes sequences that do not belong to the 18S (e.g. functional genes).
4. **Internal gaps**: Removes sequences that contain large internal gaps.
5. **Chimeras**: Removes chimeric sequences that derive from two different sequences in the dataset.

Additional step when cleaning more than one dataset:

6. **Merging and clustering**: Merges the different datasets into one single seqtab and clusters identical sequences (if a 100% identity is used) that may have longer ends or terminal deletions.

## How to run this workflow

### Step 1: clone this repository

Run this command:

```
git clone https://github.com/aleixop/clean_amplicon_datasets.git
```

### Step 2: install required software

Required software is pretty common in bioinformatic analyses, so your cluster may already have them all installed. If not, software can be installed either manually or through conda.

#### Manual installation

The required software for this pipeline is the following:

- [`snakemake`](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [`BLAST`](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- [`R`](https://cran.r-project.org/)
- [`seqkit`](https://bioinf.shenwei.me/seqkit/download/)
- [`mafft`](https://mafft.cbrc.jp/alignment/software/)
- [`hmmer`](http://hmmer.org/download.html)

And the following R packages (you can just copy the following chunk into the R terminal):

```
install.packages(c("tidyverse","argparser","knitr"), repos = "https://cloud.r-project.org")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings","IRanges","dada2","DECIPHER"))

```

#### Conda installation

First of all, if you don't have mamba installed, follow the steps explained [here](https://snakemake.readthedocs.io/en/stable/tutorial/setup.html#step-1-installing-miniforge) to do so. 

Alternatively, you can install mamba like this:

```
conda activate base
conda install -n base -c conda-forge mamba
```

Then, activate the conda base environment, enter the `clean_amplicon_datasets` directory and install all required software into an isolated Conda environment with the name `clean_amplicon_datasets` via:

```
conda activate base
cd clean_amplicon_datasets
mamba env create --name clean_amplicon_datasets --file environment.yaml
```

### Step 3: prepare your input files

In case you want to test the pipeline, this repository contains files for testing in `data/input/`. To run the pipeline on your own samples just remove these files and add your seqtab (or seqtabs) to `data/input/`. These should follow this naming:

```
data/input/<dataset1>.rds
data/input/<dataset2>.rds
...
```

### Step 4: modify parameters

You can modify the workflow parameters by editing the config file `config.yaml`, located in the root of the project. Remember that the default parameters are fine-tuned for eukaryotic datasets of the V4 region of the 18S rRNA gene. 

### Step 5: run the pipeline

#### With manual installation

Load all the [required software](#manual-installation) or make sure that paths for software are exported and run this code from the root of the project (where the `Snakefile` is located). You can write the number of threads you want to use with `--cores`:

```
snakemake --cores <threads>
```

If using an HPC with SLURM, you can find a template to run this pipeline [here](scripts/clean_amplicon_datasets.sh). Just modify it for your case and run it with:

```
sbatch clean_amplicon_datasets.sh
```

#### With conda installation

Activate the environment you created in [Step 2](#conda-installation):

```
conda activate clean_amplicon_datasets
```

And run this code from the root of the project (where the `Snakefile` is located). You can write the number of threads you want to use with `--cores`:

```
snakemake --cores <threads>
```

## Pipeline steps in detail 

### Filtering

This step takes the initial seqtab and filters it by thresholds of abundance (default is to remove singletons), occurrence and minimum bp length. It uses the script in `filter_seqtab.R`:

```
usage: filter_seqtab.R [--] [--help] [--seqtab SEQTAB] [--min_abundance
       MIN_ABUNDANCE] [--min_occurrence MIN_OCCURRENCE] [--min_bp
       MIN_BP] [--out_seqtab OUT_SEQTAB] [--out_removed_seqs
       OUT_REMOVED_SEQS]

Filter a dada2 seqtab file by minimum abundance, occurrence and length.

flags:
  -h, --help           show this help message and exit

optional arguments:
  -s, --seqtab         Path to the seqtab file.
  -m, --min_abundance  ASVs with less than this number of overall reads
                       will be removed. [default: 1]
  --min_occurrence     ASVs present in less than this number of samples
                       will be removed. [default: 1]
  --min_bp             ASVs shorter than this number of base pairs will
                       be removed. [default: 32]
  -o, --out_seqtab     Path where filtered seqtab file should be
                       written. [default: filtered_seqtab.rds]
  --out_removed_seqs   Path where file with removed ASVs should be
                       written. [default: removed_seqs.txt]
```

### Clustering

If the input seqtab was created by merging different datasets, it may contain identical sequences with longer ends or terminal deletions. The idea here is to cluster the seqtab at 100% identity to merge those sequences. By default, the most abundant sequence is chosen as representative. Clustering is performed by the script `cluster_seqtab.R`, using DECIPHER's function `Clusterize`:

```
usage: cluster_seqtab.R [--] [--help] [--seqtab SEQTAB]
       [--perc_identity PERC_IDENTITY] [--min_coverage MIN_COVERAGE]
       [--representative_method REPRESENTATIVE_METHOD] [--out_seqtab
       OUT_SEQTAB] [--out_clusters OUT_CLUSTERS] [--out_removed_seqs
       OUT_REMOVED_SEQS]

Cluster a dada2 seqtab file base on percentage identity and minimum
coverage.

flags:
  -h, --help                   show this help message and exit

optional arguments:
  -s, --seqtab                 Path to the seqtab file.
  -p, --perc_identity          Clustering identity in a 0-100 scale.
                               [default: 100]
  -m, --min_coverage           Minimum fraction of sequence positions
                               that must be overlapping for a sequence
                               to be clustered with the cluster
                               representative. [default: 0.9]
  -r, --representative_method  Choose whether representatives should be
                               chosen based on 'abundance' (most
                               abundant ASV as representative) or
                               'length' (longest ASV as
                               representative). [default: abundance]
  -o, --out_seqtab             Path where filtered seqtab file should
                               be written. [default:
                               clustered_seqtab.rds]
  --out_clusters               Path where correspondence table between
                               ASVs and cluster representatives should
                               be written. [default: clusters.tsv]
  --out_removed_seqs           Path where file with removed ASVs should
                               be written. [default: removed_seqs.txt]
```

### HMM search to remove non-ribosomal sequences

Some amplicon datasets may contain sequences that do not belong to the 18S (e.g., sequences derived from functional genes). In this step, a fasta is created from the seqtab file obtained from previous step. This fasta is then aligned with `mafft` (auto mode) and a HMM profile is built from the obtained aligment (with HMMER's `hmmbuild`). Then, `hmmsearch` is performed using the fasta file against the HMM profile. Sequences that do no map or have a high evalue are removed by the script `process_hmmer.R`. The idea behind this step is that the majority of sequences in the dataset should belong to the 18S. Thus the HMM profile should capture the 18S structure and sequences not belonging to the 18S should not properly align to it.

```
usage: process_hmmer.R [--] [--help] [--seqtab SEQTAB] [--hmmer HMMER]
       [--min_evalue MIN_EVALUE] [--out_seqtab OUT_SEQTAB]
       [--out_removed_seqs OUT_REMOVED_SEQS]

Remove ASVs not belonging to 18S from a seqtab based on a hmmsearch
tblout table

flags:
  -h, --help          show this help message and exit

optional arguments:
  -s, --seqtab        Path to the seqtab file.
  --hmmer             Path to the hmmsearch file.
  -m, --min_evalue    Minimum e-value to filter hmmsearch results.
                      [default: 1e-05]
  -o, --out_seqtab    Path where filtered seqtab file should be
                      written. [default: filtered_seqtab.rds]
  --out_removed_seqs  Path where file with removed ASVs should be
                      written. [default: removed_seqs.txt]
```

### Remove internal gaps

Some sequences may contain internal gaps. In this step, an all versus all BLAST search is performed to look for sequences that have splitted hits to one or more sequences in the dataset. The obtained blast is processed with the script `process_blast_internal_gaps.R`.

```
usage: process_blast_internal_gaps.R [--] [--help] [--seqtab SEQTAB]
       [--blast BLAST] [--fasta FASTA] [--min_query_coverage
       MIN_QUERY_COVERAGE] [--min_gaps_subject MIN_GAPS_SUBJECT]
       [--diff_ranks_log DIFF_RANKS_LOG] [--out_seqtab OUT_SEQTAB]
       [--out_removed_seqs OUT_REMOVED_SEQS]

Remove ASVs with internal gaps from a seqtab based on a blast search

flags:
  -h, --help                show this help message and exit

optional arguments:
  -s, --seqtab              Path to the seqtab file.
  -b, --blast               Path to the blast file.
  -f, --fasta               Path to the fasta file used for the blast
                            search.
  -m, --min_query_coverage  Minimum fraction of the query sequence that
                            must be aligned for removal consideration.
                            [default: 0.99]
  --min_gaps_subject        Minimum length of internal gaps in the
                            subject sequence for query removal.
                            [default: 15]
  -d, --diff_ranks_log      Maximum difference between log-transformed
                            ranks [log(rank_query) - log(rank_subject)]
                            to remove the subject instead of the query.
                            This ensures that query instead of subject
                            is kept when the abundance rank of query is
                            considerably lower (i.e., query is more
                            abundant than subject). [default: -0.4]
  -o, --out_seqtab          Path where filtered seqtab file should be
                            written. [default: filtered_seqtab.rds]
  --out_removed_seqs        Path where file with removed ASVs should be
                            written. [default: removed_seqs.txt]
```

### Chimera removal

This step detects chimeras (bimeras) by looking for the 2 parents that originated them. It consists of 2 rounds of blast. The first one maps the first 180 bp of each sequence against the complete sequences and looks for exact matches. For each query, it keeps the subject with highest rank. This is done by script `process_chimeras_1.R`:

```
usage: process_chimeras_1.R [--] [--help] [--seqtab SEQTAB] [--blast
       BLAST] [--fasta FASTA] [--out_list OUT_LIST]

Take a blast of seqs with trimmed tails against the full sequences and
output putative chimeras.

flags:
  -h, --help      show this help message and exit

optional arguments:
  -s, --seqtab    Path to the seqtab file.
  -b, --blast     Path to the blast file.
  -f, --fasta     Path to the fasta file used for the blast search.
  -o, --out_list  Path where ASVs to consider in the second chimera
                  step should be written. [default: out_list.txt]
```

Then, for the sequences that mapped, a new fasta is created by keeping only their tails (from 241 to end). This sequences are then blasted against the full sequences again. Then, the identity between the putative chimera and each of their parents and the identity between parents are computed. The idea is that a chimera should be different from parents and that parents should be relatively distant. Also, occurrence is considered, as chimeras should be present in a low number of samples. All this processing is performed by script `process_chimeras_2.R`:

```
usage: process_chimeras_2.R [--] [--help] [--seqtab SEQTAB] [--blast1
       BLAST1] [--blast2 BLAST2] [--fasta FASTA] [--out_seqtab
       OUT_SEQTAB] [--out_removed_seqs OUT_REMOVED_SEQS]
       [--out_chimeras OUT_CHIMERAS] [--max_pident_chimera_parent
       MAX_PIDENT_CHIMERA_PARENT] [--max_pident_parents
       MAX_PIDENT_PARENTS] [--max_chimera_occurrence
       MAX_CHIMERA_OCCURRENCE]

Take 2 blast of seqs with (1) trimmed tails and (2) trimmed heads
against the full sequences. Detect chimeras and filter seqtab.

flags:
  -h, --help                       show this help message and exit

optional arguments:
  -s, --seqtab                     Path to the seqtab file.
  -b, --blast1                     Path to the blast file with trimmed
                                   tails.
  --blast2                         Path to the blast file with trimmed
                                   heads.
  -f, --fasta                      Path to the fasta file used for the
                                   blast search.
  -o, --out_seqtab                 Path where filtered seqtab file
                                   should be written. [default:
                                   filtered_seqtab.rds]
  --out_removed_seqs               Path where file with removed ASVs
                                   should be written. [default:
                                   removed_seqs.txt]
  --out_chimeras                   Path where chimeras report should be
                                   written. [default: chimeras.txt]
  -m, --max_pident_chimera_parent  Maximum percentage of sequence
                                   identity between a chimera and both
                                   of their parents. [default: 99]
  --max_pident_parents             Maximum percentage of sequence
                                   identity between the parents of a
                                   chimera. [default: 95]
  --max_chimera_occurrence         Maximum occurrence (in percentage)
                                   that a chimera can have. [default:
                                   5]

```

### Merge datasets and clustering

The script automatically detects if there is more than one dataset in the input. If so, it takes the final seqtabs obtained for each datasets and merges them with the script `merge_seqtabs.R`:

```
usage: merge_seqtabs.R [--] [--help] [--seqtabs SEQTABS] [--out_seqtab
       OUT_SEQTAB]

Merge seqtabs into one single file.

flags:
  -h, --help        show this help message and exit

optional arguments:
  -s, --seqtabs     Path to seqtab files
  -o, --out_seqtab  Path where the merged seqtab should be written.
                    [default: merged_seqtab.rds]
```

Then, the obtained seqtab is clustered again to collapse identical sequences by the script `cluster_seqtab.R` and the final seqtab is obtained.

## Output explained

Here is an example of the output generated with the test data provided in `data/input/`:

```
results/
│
├── datasets/
│   ├── dataset1/
│   │   ├── 1-filtering/
│   │   │   ├── dataset1_filtering_seqtab.rds
│   │   │   ├── dataset1_filtering_removed_seqs.txt
│   │   ├── 2-clustering/
│   │   │   ├── dataset1_clustering_seqtab.rds
│   │   │   ├── dataset1_clusters.tsv
│   │   │   ├── dataset1_clustering_removed_seqs.txt
│   │   ├── 3-hmmsearch/
│   │   │   ├── dataset1_hmmsearch_seqtab.rds
│   │   │   ├── dataset1_hmmsearch_removed_seqs.txt
│   │   ├── 4-internal_gaps/
│   │   │   ├── dataset1_internal-gaps_seqtab.rds
│   │   │   ├── dataset1_internal-gaps_removed_seqs.txt
│   │   ├── 5-chimeras/
│   │   │   ├── dataset1_chimeras_seqtab.rds
│   │   │   ├── dataset1_removed_seqs.txt
│   │   │   ├── dataset1_chimeras.txt
│   │   ├── dataset1_final_seqtab.rds  (Symlink to seqtab in 5-chimeras)
│   │   ├── dataset1_removed_seqs.txt
│   │   ├── dataset1_report.tsv
│   │
│   ├── dataset2/
│   │   ├── (Same structure as dataset1)
│
└── final/
    ├── clusters.tsv
    ├── final_seqtab.rds
    ├── overall_report.tsv
```

All results are stored in `results/`. For each dataset, a subdirectory is created in `results/datasets/`. Inside the dataset directory, directories for each cleaning step are also created. Inside each step directory, the filtered seqtab from that specific step and a list of the removed sequences are written. Additionally, for the clustering step a file with clusters is created; and for the chimera removal step, a file with chimeras is also written. Inside the root directory of each dataset, the final filtered seqtab, all removed sequences in all steps and a report of the the process (ASVs and reads removed in each step) are given.

When working with more than one dataset, a directory `results/final` is created. There, the final clustered seqtab, the clusters file and an overall report are written. 

## Additional help

If you need some help to obtain the seqtabs that are used as input in this pipeline, you can take a look at our [dada2_guidelines](https://github.com/adriaaula/dada2_guidelines).