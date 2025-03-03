#!/bin/bash

#SBATCH --job-name=clean-amplicons
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=clean-amplicons_%J.out
#SBATCH --error=clean-amplicons_%J.err

## Load modules if you did not use conda for installation

module load mafft
module load hmmer
module load blast
module load seqkit
module load R
module load snakemake

## Activate the environment if you used conda for installation

# module load miniconda # you may need to change this for your cluster
# conda activate clean_amplicon_datasets

## Run the pipeline

snakemake --cores ${SLURM_CPUS_PER_TASK}