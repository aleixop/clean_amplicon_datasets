#!/bin/bash

#SBATCH --job-name=clean-amplicons
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=clean-amplicons_%J.out
#SBATCH --error=clean-amplicons_%J.err

module load mafft
module load hmmer
module load blast
module load seqkit
module load R
module load snakemake

snakemake --cores ${SLURM_CPUS_PER_TASK}