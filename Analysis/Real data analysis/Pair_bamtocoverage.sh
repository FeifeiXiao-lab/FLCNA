#!/bin/sh

#SBATCH --job-name=Pair_bamtocoverage
#SBATCH --output=job%j.out
#SBATCH --error job%j.err
#SBATCH -p gcai-lab
###Number of Cores Max 20
#SBATCH -n 2
###Number of Nodes
#SBATCH -N 1

module load python3/anaconda/2020.02
export R_LIBS_USER=/home/fqin/Packages/R/4.1.0
source activate /work/fqin/R_env
Rscript /home/fqin/SRA001/fastq/Pair_bamtocoverage.R

