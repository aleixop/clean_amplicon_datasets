# Snakemake workflow to clean amplicon datasets

This is a pipeline to clean short eukaryotic amplicon datasets processed with dada2. It has been fine-tuned for sequences of the V4 region of the 18S ribosomal RNA gene. If you need some help using dada2 you can take a look at our [dada2_guidelines](https://github.com/adriaaula/dada2_guidelines).

## Summary

- [Overview](#overview)
- [How to run this workflow](#how-to-run-this-workflow)
  - [Step 1: clone this repository](#step-1-clone-this-repository)
  - [Step 2: install required software](#step-2-install-required-software)
    - [Manual installation](#manual-installation)
    - [Conda installation](#conda-installation)
  - [Step 3: prepare your input files](#step-3-prepare-your-input-files)
  - [Step 4: run the pipeline](#step-4-run-the-pipeline)
    - [With manual installation](#with-manual-installation)
    - [With conda installation](#with-conda-installation)
- [Pipeline steps in detail](#pipeline-steps-in-detail)

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

### Step 4: run the pipeline

#### With manual installation

Load all the [required software](#manual-installation) or make sure that paths for software are exported and run this code from the root of the project (where the `Snakefile` is located). You must write the number of threads you want to use:

```
snakemake --cores <threads>
```

If using an HPC with SLURM, you can find a template to run this pipeline [here](scripts/clean_amplicon_datasets.sh)

#### With conda installation

Activate the environment you created in [Step 2](#conda-installation):

```
conda activate clean_amplicon_datasets
```

And run this code from the root of the project (where the `Snakefile` is located). You must write the number of threads you want to use:

```
snakemake --cores <threads>
```

## Pipeline steps in detail