#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 8
#SBATCH --mem 32GB

###### MANUAL STEPS AND PREREQUISITES ######
## run from project folder (which contains archives input script tmp finalresult)
## if you have configured virsorter; run STEP 000 (either run interactively or delete leading # and run via sbatch)
##### command line arguments required #####
## $1 = genome assembly for analysis (.fa.gz)

################### FINAL OUTPUTS ##########################

################### INPUTS ##########################
###### load software ######
module load virsorter/2.2; module load seqkit/0.14.0

##### input files for analysis #####


##### NCVOG defintion files #####

############################### SCRIPT ###############################################

############### STEP 000 ############################
####### SETUP DATABASE (only needs to be done once)
#virsorter setup -d tmp/virsorter_db -j 4


############### STEP 001 ############################

###### run virsorter
##### virsorter is not optimised for eukaryotic genomes - even for artificial Ectocarpus contig with EsV-1 genome inserted in it -
##### viral boundaries are ~ 100 kb or more off actual boundaries
### testing trimmed/untrimmed contigs and --viral-gene-required ON/OFF had no effect on output
cd finalresult/IVEX001.5
### untrimmed, --viral-gene-required OFF
virsorter run --rm-tmpdir --min-length 1500 --include-groups NCLDV -j 8 -w $(basename $1 .fa.gz).out -i $1 all

sed '1,1d' $(basename $1 .fa.gz).out/final-viral-boundary.tsv | sort -t $'\t' -k 1,1V | cut -f 1,4,5 | awk -F"\t" 'BEGIN{OFS="\t"};{print $1,"virsorter2","virus_region",$2,$3,".\t.\t.\tID=VS2_"$1"_"$2"-"$3";Name=size_"$3-$2}' > IVEX001.5_final_$(basename $1 .fa.gz).gff
