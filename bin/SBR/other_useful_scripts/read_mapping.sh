#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition long
#SBATCH --cpus-per-task 16
#SBATCH --mem 64GB

###### MANUAL STEPS AND PREREQUISITES ######
## this is to realign reads to a contig/assembly/etc to check possible sequence errors

####################### INPUTS #####################################
###### load software ######
module load minimap2/2.18; module load samtools/1.9

##### command line arguments required #####
## $1 = reference fasta query (for reads to be aligned against)
## IF YOU HAVE MULTIPLE ILLUMINA RUNS FOR SAME GENOME, CAT THEM FIRST BY FORWARD AND REVERSE
## $2 = forward (1) reads OR unpaired/nanopore/etc reads
## $3 = reverse (2) reads

##### variables #####
n1=$(basename "$1" | cut -d "." -f 1)
n2=$(basename "$2" | cut -d "." -f 1)
n3=$(basename "$3" | cut -d "." -f 1)

####################### SCRIPT #####################################

################ STEP 001 ##############
### Illumina paired end alignment (variables c2 and c3 to merge all forward runs togther and all reverse runs together)
#minimap2 -ax sr $1 $2 $3 > minimap2_"$n1"_"$n2"_"$n3".sam
#samtools view -bS -o minimap2_"$n1"_"$n2"_"$n3".bam minimap2_"$n1"_"$n2"_"$n3".sam
#samtools sort minimap2_"$n1"_"$n2"_"$n3".bam -o minimap2_"$n1"_"$n2"_"$n3"_sorted.bam
#samtools index minimap2_"$n1"_"$n2"_"$n3"_sorted.bam

## Oxford nanopore alignment
#minimap2 -ax map-ont $1 $2 > minimap2_"$n1"_"$n2".sam
#samtools view -bS -o minimap2_"$n1"_"$n2".bam minimap2_"$n1"_"$n2".sam
#samtools sort minimap2_"$n1"_"$n2".bam -o minimap2_"$n1"_"$n2"_sorted.bam
#samtools index minimap2_"$n1"_"$n2"_sorted.bam


samtools index $1 > ./$(basename "$1").bai
