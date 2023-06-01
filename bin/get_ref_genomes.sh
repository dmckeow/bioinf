#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --ntasks=8
#SBATCH --mem=72GB
#SBATCH --tmp=72GB
#SBATCH --partition amdlarge
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err


####################### PREREQUISITES #####################################
## nanopore seq reads data
## rm -f slurm* tmp* ## cleanup between tests
####################### SCRIPT PURPOSE #####################################

####################### HOW TO RUN SCRIPT #####################################

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################
## $1 = input file directory (this will also be location of genomes)
## $2 = taxid to search for (can be any level)
####################### OUTPUT FILES #####################################
## in placed in directory flye at location of input file

####################### LOAD SOFTWARE #####################################
module load minimap2/2.17
module load samtools
CANU="/home/dcschroe/dcschroe/dmckeow/canu-2.2/bin/canu"
KRAKEN="/panfs/roc/msisoft/kraken/2.0.7beta/kraken2"
KDB="/panfs/roc/msisoft/kraken/kraken_db"
SEQKIT="/home/dcschroe/dcschroe/dmckeow/seqkit"
####################### SET VARIABLES/ARGUMENTS #####################################

################################################################################################
####################### RUNNING THE SCRIPT #####################################################
################################################################################################
cd $1
rm -f new_taxdump.tar.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar -zxf new_taxdump.tar.gz
SEARCH="$2";
grep -w "$SEARCH" taxidlineage.dmp | cut -f 1 | sort -Vu | grep -w -f - rankedlineage.dmp > house_ref_genomes1

#### manually check house_ref_genomes1 and remove all entries you dont want, and save as house_ref_genomes2

cut -f 1 house_ref_genomes2 | sort -Vu | sed -E 's/^|$/_____/g' >  house_ref_genomes3


#### get fastas for matches to reads from NCBI
rm -f assembly_summary_refseq.txt*; wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
sed -i -e 's/\t/\t_____/5' -e 's/\t/_____\t/6' assembly_summary_refseq.txt
awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/[^A-Za-z0-9_]/,"_",$8)}1' assembly_summary_refseq.txt | awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/_+/,"_",$8)}1' | awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/^_|_$/,"_",$8)}1' > tmp && mv tmp assembly_summary_refseq.txt

grep -wF -f house_ref_genomes3 assembly_summary_refseq.txt | awk -F "\t" '!a[$6]++' | cut -f 20 | sed -E 's/^https(.*)(\/.*)/wget ftp\1\2\2_genomic.fna.gz/g' > tmp_get_ncbi_genomes
rm -f *.fna.gz*; bash tmp_get_ncbi_genomes

grep -wF -f house_ref_genomes3 assembly_summary_refseq.txt | awk -F "\t" '!a[$6]++' | awk -F"\t" '{print "mv "$1"_"$16"_genomic.fna.gz "$6"."$8".ref1.fna.gz"}' | sed 's/_____//g' > ref1_namefix
bash ref1_namefix
