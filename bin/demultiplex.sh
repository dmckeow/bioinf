#!/bin/bash -l
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=240GB
#SBATCH --partition ag2tb
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

####################### PREREQUISITES #####################################
## nanopore seq reads data

####################### SCRIPT PURPOSE #####################################

####################### HOW TO RUN SCRIPT #####################################

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################
### $1 = directory with fastq.gz
### $2 = name
####################### OUTPUT FILES #####################################

####################### LOAD SOFTWARE #####################################
### LOAD software available via shared environment on server:
module load porechop
module load pigz

zcat ${1}/*fastq.gz > ${2}.fastq

pigz -p $SLURM_CPUS_PER_TASK ${2}.fastq

porechop -t $SLURM_CPUS_PER_TASK -i ${2}.fastq.gz -b ./demultiplex_${2}_unclassified_STRICT --barcode_threshold 85 --require_two_barcodes

porechop -t $SLURM_CPUS_PER_TASK -i ${2}.fastq.gz -b ./demultiplex_${2}_unclassified_LENIENT --barcode_threshold 60 --barcode_diff 1
