#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 8
#SBATCH --mem 64GB

######### databases
RVDB="/projet/fr2424/sib/dmckeown/db/virus/C-RVDBv22_0"
NR="/db/nt/current/blast/nt"
NRP="/scratch2/fr2424/sib/dmckeown/db/nr.dmnd"

####################### INPUTS #####################################
###### load software ######
module load blast/2.9.0; module load diamond/2.0.9

#### for finding repeats
#blastn -query $1 -db $RVDB -out $(basename "$1" .fa).blast -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore salltitles" -task dc-megablast -evalue 0.001 -max_target_seqs 10 -num_threads 6 -dust no -soft_masking no


##### for blastp, excluding viruses --taxonlist to restrict to taxa
diamond blastp -b8 -c1 -k 50 -e 1e-3 --unal 1 --taxonlist 33634 -d $NRP -q $1 -o $(basename $1 .faa).NR.blast -f 6 qseqid salltitles evalue bitscore staxids;
