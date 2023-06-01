#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=8
#SBATCH --mem=72GB
#SBATCH --tmp=72GB
#SBATCH --partition agsmall
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err


####################### PREREQUISITES #####################################
## nanopore seq reads data
## rm -f slurm* tmp* ## cleanup between tests
####################### SCRIPT PURPOSE #####################################

####################### HOW TO RUN SCRIPT #####################################
#### make a directory for all results

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################
## $1 = input file directory (this will also be location of results)
####################### OUTPUT FILES #####################################
## in placed in directory flye at location of input file

####################### LOAD SOFTWARE #####################################
module load minimap2/2.17
module load samtools
module load porechop/0.2.4
CANU="/home/dcschroe/dcschroe/dmckeow/canu-2.2/bin/canu"
KRAKEN="/panfs/roc/msisoft/kraken/2.0.7beta/kraken2"
KDB="/panfs/roc/msisoft/kraken/kraken_db"
SEQKIT="/home/dcschroe/dcschroe/dmckeow/seqkit"
####################### SET VARIABLES/ARGUMENTS #####################################
HOUSE="$HOME/data/house_ref_genomes/*.fna.gz" ## house reference genomes location
OUTD="$1/quick_ONT_checkmap"
################################################################################################
####################### RUNNING THE SCRIPT #####################################################
################################################################################################

mkdir $OUTD
rm -f $1/tmp_allreads.fastq
cat $1/*.fastq > $1/tmp_allreads.fastq
porechop -i $1/tmp_allreads.fastq -o $OUTD/tmp_allreads_trimmed.fastq

$KRAKEN --db $KDB $OUTD/tmp_allreads_trimmed.fastq --use-mpa-style --threads 8 --use-names --output $OUTD/"$(basename $1)".kraken2 --report $OUTD/"$(basename $1)".kraken2.report --unclassified-out $OUTD/"$(basename $1)".kraken2.unclassified --classified-out $OUTD/"$(basename $1)".kraken2.classified

#### get fastas for matches to reads from NCBI
rm -f $OUTD/assembly_summary_refseq.txt*; wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -P $OUTD
sed -i -e 's/\t/\t_____/5' -e 's/\t/_____\t/6' $OUTD/assembly_summary_refseq.txt
awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/[^A-Za-z0-9_]/,"_",$8)}1' $OUTD/assembly_summary_refseq.txt | awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/_+/,"_",$8)}1' | awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/^_|_$/,"_",$8)}1' > $OUTD/tmp && mv $OUTD/tmp $OUTD/assembly_summary_refseq.txt

cut -f 3 $OUTD/"$(basename $1)".kraken2 | sed -E 's/.+taxid ([0-9]+)\)/_____\1_____/g' | sort -Vu | grep -wF -f - $OUTD/assembly_summary_refseq.txt | awk -F "\t" '!a[$6]++' | cut -f 20 | sed -E 's/^https(.*)(\/.*)/wget ftp\1\2\2_genomic.fna.gz -P $OUTD/g' > $OUTD/tmp_get_ncbi_genomes
rm -f $OUTD/*.fna.gz*; bash $OUTD/tmp_get_ncbi_genomes

a="$(basename $1)"; cut -f 3 $OUTD/"$(basename $1)".kraken2 | sed -E 's/.+taxid ([0-9]+)\)/_____\1_____/g' | sort -Vu | grep -wF -f - $OUTD/assembly_summary_refseq.txt | awk -F "\t" '!a[$6]++' | awk -F"\t" -v a=$a -v b=$OUTD '{print "mv "b$1"_"$16"_genomic.fna.gz "b"tmp_"a"."$6"."$8".ref1.fna.gz"}' | sed 's/_____//g' > $OUTD/ref1_namefix
bash $OUTD/ref1_namefix
