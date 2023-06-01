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
#### make job submission file for convenience
#find . -name "barcode*" | awk '/fastq_pass\/barcode[0-9]+$/' > tmp; for f in $(cat tmp); do realpath $f; done | sed -E 's/(.+)/sbatch \/home\/dcschroe\/dmckeow\/projects\/DWV\/script\/viral_pipeline1.sh \1/g' > run_viral_pipeline1
## to check if reads were processed
#find . -name "barcode*" | awk '/fastq_pass\/barcode[0-9]+$/' > tmp; for f in $(cat tmp); do realpath $f; done | sed -E 's/(.+)/ls \1\/DENOVO*/g' | bash - | grep ".trimmedReads.fasta.gz" > tmp_trimread_success
####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################
## $1 = input file directory (this will also be location of results)

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
HOUSE="$HOME/data/house_ref_genomes/*.fna.gz" ## house reference genomes location
################################################################################################
####################### RUNNING THE SCRIPT #####################################################
################################################################################################
######### classify the DENOVO contigs
sed -i 's/ /___/g' $1/DENOVO_"$(basename $1)"/"$(basename $1)".contigs.fasta
sed -i 's/ /___/g' $1/DENOVO_"$(basename $1)"/"$(basename $1)".unassembled.fasta

cat $1/DENOVO_"$(basename $1)"/"$(basename $1)".contigs.fasta $1/DENOVO_"$(basename $1)"/"$(basename $1)".unassembled.fasta | sed 's/ /___/g' > $1/DENOVO_"$(basename $1)"/"$(basename $1)".contigsall.fasta

$KRAKEN --db $KDB $1/DENOVO_"$(basename $1)"/"$(basename $1)".contigsall.fasta --use-mpa-style --threads 8 --use-names --output $1/DENOVO_"$(basename $1)"/"$(basename $1)".contigs.kraken2 --report $1/DENOVO_"$(basename $1)"/"$(basename $1)".contigs.report.kraken2 --unclassified-out $1/DENOVO_"$(basename $1)"/"$(basename $1)".contigs.unclassified.kraken2 --classified-out $1/DENOVO_"$(basename $1)"/"$(basename $1)".contigs.classified.kraken2

#################### STEP 001 #################################
##### REFERENCE ASSEMBLY
#### run kraken to classify reads
rm -fr REF_"$(basename $1)"; mkdir REF_"$(basename $1)"
cd $1/REF_"$(basename $1)"/

$KRAKEN --db $KDB $1/DENOVO_"$(basename $1)"/"$(basename $1)".trimmedReads.fasta.gz --use-mpa-style --threads 8 --use-names --output $1/REF_"$(basename $1)"/"$(basename $1)".trimmedReads.kraken2 --report $1/REF_"$(basename $1)"/"$(basename $1)".trimmedReads.report.kraken2 --unclassified-out $1/REF_"$(basename $1)"/"$(basename $1)".trimmedReads.unclassified.kraken2 --classified-out $1/REF_"$(basename $1)"/"$(basename $1)".trimmedReads.classified.kraken2

#### get fastas for matches to reads from NCBI
rm -f $1/assembly_summary_refseq.txt*; wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt -P $1
sed -i -e 's/\t/\t_____/5' -e 's/\t/_____\t/6' $1/assembly_summary_refseq.txt
awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/[^A-Za-z0-9_]/,"_",$8)}1' $1/assembly_summary_refseq.txt | awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/_+/,"_",$8)}1' | awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/^_|_$/,"_",$8)}1' > $1/tmp && mv $1/tmp $1/assembly_summary_refseq.txt

cut -f 3 $1/REF_"$(basename $1)"/"$(basename $1)".trimmedReads.kraken2 | sed -E 's/.+taxid ([0-9]+)\)/_____\1_____/g' | sort -Vu | grep -wF -f - $1/assembly_summary_refseq.txt | awk -F "\t" '!a[$6]++' | cut -f 20 | sed -E 's/^https(.*)(\/.*)/wget ftp\1\2\2_genomic.fna.gz/g' > $1/tmp_get_ncbi_genomes
rm -f $1/*.fna.gz*; bash $1/tmp_get_ncbi_genomes

a="$(basename $1)"; b="$1/"; cut -f 3 "$(basename $1)".trimmedReads.kraken2 | sed -E 's/.+taxid ([0-9]+)\)/_____\1_____/g' | sort -Vu | grep -wF -f - assembly_summary_refseq.txt | awk -F "\t" '!a[$6]++' | awk -F"\t" -v a=$a -v b=$b '{print "mv "b$1"_"$16"_genomic.fna.gz "b"tmp_"a"."$6"."$8".ref1.fna.gz"}' | sed 's/_____//g' > $1/ref1_namefix
bash $1/ref1_namefix

####### copy house genomes to working directory (these are genomes you have previously chosen to always be used as reference genomes)
for f in $HOUSE; do cp $f $1/tmp_"$(basename $1)"."$(basename $f)"; done

#### map ALL reads against reference genomes identified by Kraken AND keep only aligned reads
for f in $1/tmp_*.ref1.fna.gz; do a=$1/DENOVO_"$(basename $1)"/"$(basename $1)".trimmedReads.fasta.gz; minimap2 -t 8 -ax map-ont $f $a | samtools sort -O BAM - > $1/REF_"$(basename $1)"/"$(basename $f .ref1.fna.gz)".bam; done
for f in $1/REF_"$(basename $1)"/tmp_*.ref1.fna.gz; do samtools view -b -F 4 $1/REF_"$(basename $1)"/"$(basename $f .ref1.fna.gz)".bam > $1/REF_"$(basename $1)"/tmp && mv $1/REF_"$(basename $1)"/tmp $1/REF_"$(basename $1)"/"$(basename $f .ref1.fna.gz)".bam; done

######## fix file names
realpath $1/REF_"$(basename $1)"/tmp_*.bam | sed -E 's/(.+)/mv \1 \1/g' | sed 's/\/tmp_/\//2' > $1/REF_"$(basename $1)"/tmp_rename
realpath $1/REF_"$(basename $1)"/tmp_*.ref1.fna.gz | sed -E 's/(.+)/ mv \1 \1/g' | sed 's/\/tmp_/\//2' | paste -d ";" $1/REF_"$(basename $1)"/tmp_rename - > $1/REF_"$(basename $1)"/tmp_manual_rename

############## MANUAL STEP #################################
## manually remove the tmp_ prefix from the .bam files you want to keep, then do rm -f tmp*
## step 002 will not run if any tmp*.bam files exist in $1
