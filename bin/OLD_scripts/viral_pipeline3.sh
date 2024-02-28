#!/bin/bash -l
#SBATCH --time=24:00:00
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
#### make job submission file for convenience
#find . -name "barcode*" | awk '/fastq_pass\/barcode[0-9]+$/' > tmp; for f in $(cat tmp); do realpath $f; done | sed -E 's/(.+)/sbatch \/home\/dcschroe\/dmckeow\/projects\/DWV\/script\/viral_pipeline3.sh \1/g' > run_viral_pipeline3

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################
## $1 = input file directory (this will also be location of results)

####################### OUTPUT FILES #####################################
## in placed in directory flye at location of input file

####################### LOAD SOFTWARE #####################################
module load minimap2/2.17
module load samtools
module load racon/1.4.20
CANU="/home/dcschroe/dcschroe/dmckeow/canu-2.2/bin/canu"
KRAKEN="/panfs/roc/msisoft/kraken/2.0.7beta/kraken2"
KDB="/panfs/roc/msisoft/kraken/kraken_db"
SEQKIT="/home/dcschroe/dcschroe/dmckeow/seqkit"
module load quast
####################### SET VARIABLES/ARGUMENTS #####################################

################################################################################################
####################### RUNNING THE SCRIPT #####################################################
################################################################################################
cd $1
cd REF_"$(basename $1)"

for f in "$(basename $1)"*.mapped.fastq; do
  sed -i 's/ /___/g' "$(basename $f .mapped.fastq)"/"$(basename $f .mapped.fastq)".contigs.fasta
  sed -i 's/ /___/g' "$(basename $f .mapped.fastq)"/"$(basename $f .mapped.fastq)".unassembled.fasta
done

#################### STEP 003 #################################
for f in "$(basename $1)"*.mapped.fastq; do
  cat "$(basename $f .mapped.fastq)"/"$(basename $f .mapped.fastq)".contigs.fasta "$(basename $f .mapped.fastq)"/"$(basename $f .mapped.fastq)".unassembled.fasta | sed 's/ /___/g' > "$(basename $f .mapped.fastq)"/"$(basename $f .mapped.fastq)".contigsall.fasta

  minimap2 -t 8 -ax map-ont "$(basename $f .mapped.fastq)"/"$(basename $f .mapped.fastq)".contigsall.fasta ../tmp_allreads.fastq > "$(basename $f .mapped.fastq)"/"$(basename $f .ref1.fna.gz)".racon.sam

  racon -u -m 8 -x -6 -g -8 -w 500 -t 8 ../tmp_allreads.fastq "$(basename $f .mapped.fastq)"/"$(basename $f .ref1.fna.gz)".racon.sam "$(basename $f .mapped.fastq)"/"$(basename $f .mapped.fastq)".contigsall.fasta > "$(basename $f .mapped.fastq)"/"$(basename $f .mapped.fastq)".contigsall.racon.fasta
done

cd ../DENOVO_"$(basename $1)"

#################### STEP 004 #################################
minimap2 -t 8 -ax map-ont "$(basename $1)".contigsall.fasta ../tmp_allreads.fastq > "$(basename $1)".racon.sam

racon -u -m 8 -x -6 -g -8 -w 500 -t 8 ../tmp_allreads.fastq "$(basename $1)".racon.sam "$(basename $1)".contigsall.fasta > "$(basename $1)".contigsall.racon.fasta


#################### STEP 005 #################################
#### assess the assembly qualities
#find REF_"$(basename $1)" -name "*.fasta" > tmp_quast_assemblies
#find DENOVO_"$(basename $1)" -name "*.fasta" >> tmp_quast_assemblies
#sed 's/$/ /g' tmp_quast_assemblies | tr -d '\n' > tmp && mv tmp tmp_quast_assemblies

##############################
#### get fastas for matches to reads from NCBI
#rm -f assembly_summary_refseq.txt*; wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
#sed -i -e 's/\t/\t_____/5' -e 's/\t/_____\t/6' assembly_summary_refseq.txt
#awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/[^A-Za-z0-9_]/,"_",$8)}1' assembly_summary_refseq.txt | awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/_+/,"_",$8)}1' | awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/^_|_$/,"_",$8)}1' > tmp && mv tmp assembly_summary_refseq.txt

#cut -f 3 REF_"$(basename $1)"/"$(basename $1)".trimmedReads.kraken2 | sed -E 's/.+taxid ([0-9]+)\)/_____\1_____/g' | sort -Vu | grep -wF -f - assembly_summary_refseq.txt | awk -F "\t" '!a[$6]++' | cut -f 20 | sed -E 's/^https(.*)(\/.*)/wget ftp\1\2\2_genomic.fna.gz/g' > tmp_get_ncbi_genomes
#rm -f *.fna.gz*; bash tmp_get_ncbi_genomes

#a="$(basename $1)"; cut -f 3 REF_"$(basename $1)"/"$(basename $1)".trimmedReads.kraken2 | sed -E 's/.+taxid ([0-9]+)\)/_____\1_____/g' | sort -Vu | grep -wF -f - assembly_summary_refseq.txt | awk -F "\t" '!a[$6]++' | awk -F"\t" -v a=$a '{print "mv "$1"_"$16"_genomic.fna.gz tmp_"a"."$6"."$8".ref1.fna.gz"}' | sed 's/_____//g' > ref1_namefix
#bash ref1_namefix
#################################

#ls *.fna.gz | sed 's/$/ /g' | tr -d '\n' > tmp && mv tmp tmp_quast_refs

#metaquast $(cat tmp_quast_assemblies) -R $(cat tmp_quast_refs)
