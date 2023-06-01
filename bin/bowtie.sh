#!/bin/bash -l
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=600GB
#SBATCH --tmp=300GB
#SBATCH --partition ag2tb
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

####################### PREREQUISITES #####################################

####################### SCRIPT PURPOSE #####################################

####################### HOW TO RUN SCRIPT #####################################

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################


####################### OUTPUT FILES #####################################

####################### LOAD SOFTWARE #####################################
module load bowtie2
mkdir Phex_SEP22_index
bowtie2-build --threads 16 Phex_SEP22.fa Phex_SEP22_index/Phex_SEP22.fa
