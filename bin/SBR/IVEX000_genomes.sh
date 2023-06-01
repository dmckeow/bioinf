#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition long
#SBATCH --cpus-per-task 8
#SBATCH --mem 64GB

###### EXAMPLE WORKFLOW ######
## This script is example of how to FIX inputs for IVEX000.sh
## This is only needed if running on genomes from outside phaeoexplorer (such as PUBLIC_ genomes)

########### input requirements: #########
##### .gff (gff3 format) - secondary gene predictions to merge with those of GeneMarkS requirements:
## tab-delimited fields 1-9, then semicolon delimited beyond, any number of fields, as long as 1-9 are present as in example
## NO lines which are empty or start with #
## [contig_name_matching_assembly_defline_exactly] in $1
## [mRNA] in $3
## [ID=gene_name_matching_protein_fasta_deflines_exactly] in $9
## any non-mRNA features can be of any format or content. The mRNA features can have any other fields beyond 1-9, but must be semicolon delimited
## gff example:
## S-promiscuus_M_contig6  Gmove   mRNA    57545   58921   324     -       .       ID=mRNA_S-promiscuus_M_contig6.15943.1
##### assembly fasta defline example:
## >S-promiscuus_M_contig6
##### protein fasta defline example:
## >mRNA_S-promiscuus_M_contig6.15943.1
#### for ease of use (simple paths to specify as arguments), first place all gff/assembly/protein fastas in single directory (even if they won't all need fixed)
#place fixed genomes in: input/IVEX000_genomes

## Here are some example lines for some public genomes which needed fixed:
## Any public gff/fastas will neeed to inspected manually, and fixed in various ways
## gzip the assembly fasta copies

##### all genomes
## get all genome assemblies
#cp /shared/projects/phaeoexplorer/01__FINAL_GENOME_ASSEMBLIES/*.fa .
## get all protein fastas and fix Phaeoexplorer deflines (public genomes still need fixed, but this still copies them)
#for f in /shared/projects/phaeoexplorer/04__FINAL_PROTEOMES/*_proteins.fa; do sed -E 's/>prot_(.+) assembled CDS.*/>mRNA_\1/g' $f > $(basename $f); done
## get all gffs and remove any empty or # lines
#for f in /shared/projects/phaeoexplorer/02__FINAL_ANNOTATIONS/*.gff; do sed -e '/^#/d' -e '/^$/d' $f > $(basename $f); done
######## specific fixes for genomes not matching format
## Undaria pinnatifida
#sed -i -Ee 's/\|/;;;/1' -Ee 's/\|/;;;/1' -Ee 's/>.+;;;(.+);;;.+/>TR.\1/g' PUBLIC_Undaria-pinnatifida_proteins.fa
## Sargassum fusiforme
#sed -i -e 's/\ttranscript\t/\tmRNA\t/g' -e 's/\t/\tID=/8' PUBLIC_Sargassum-fusiforme.gff
## Saccharina_japonica
#sed -i -e 's/\tmRNA\t/\t;;;;;\t/g' -e 's/\tgene\t/\tmRNA\t/g' -e 's/\t;;;;;\t/\tgene\t/g' PUBLIC_Saccharina-japonica.gff
#sed -i -e 's/ /;;;/1' -Ee 's/(>.+);;;.+/\1/g' PUBLIC_Saccharina-japonica_proteins.fa
## Ectocarpus_sp7
#sed -i -e 's/ /;;;/1' -Ee 's/(>.+);;;.+/\1/g' PUBLIC_Ectocarpus-sp7_proteins.fa
## Nemacystus decipiens
#sed -i -e 's/\ttranscript\t/\tmRNA\t/g' PUBLIC_Nemacystus-decipiens.gff
## Cladosiphon okamuranus
#sed -i -e 's/\ttranscript\t/\tmRNA\t/g' PUBLIC_Cladosiphon-okamuranus.gff
#sed -i -E 's/>Cok_S_s[0-9]+_([0-9]+\..+)/>g\1/g' PUBLIC_Cladosiphon-okamuranus_proteins.fa


################ GET READ fastqs from the web
### $2 = https to SRA archive
### $3 = path to SRA archive file
### $4 = output directory name

if [[ "$1" =~ "A1" ]]; then

module load sra-tools/2.11.0
cd $2
wget $3
rm -fr $5
mkdir $5
fasterq-dump --threads $SLURM_CPUS_PER_TASK --skip-technical --split-files --outdir $5 $4
pigz -p $SLURM_CPUS_PER_TASK $5/*.fastq
fi

############ MAP reads vs contig_size
### $2 = contigs.fasta
### $3 = illumina reads .1
### $4 = illumina reads .2
### $5 = output name (unique for each read pair)
### $6 = output directory (BIG storage NEEDED)

if [[ "$1" =~ "index" ]]; then
#bwa index -p $2 $2
module load bowtie2
bowtie2-build --threads $SLURM_CPUS_PER_TASK $2 $2
fi


if [[ "$1" =~ "A2" ]]; then

#module load bwa-mem2
module load samtools
module load bowtie2

#bwa mem -M -t 8 $2 $3 $4 > ${5}.sam
#samtools fixmate -O bam ${5}.sam ${5}.bam
#samtools sort -T tmp.${5}.prefix -O bam -@ 8 -o ${5}.sorted.bam ${5}.bam

bowtie2 --threads $SLURM_CPUS_PER_TASK -x $2 -1 $3 -2 $4 --no-unal -S ${6}/${5}.sam
samtools view -@ $SLURM_CPUS_PER_TASK -F 4 -bS ${6}/${5}.sam > ${6}/${5}.RAW.bam
rm -f ${6}/${5}.sam
samtools sort -T ${6}/tmp.${5}.RAW.bam -O bam -@ $SLURM_CPUS_PER_TASK -o ${6}/${5}.bam ${6}/${5}.RAW.bam
rm -f ${6}/${5}.RAW.bam

fi

if [[ "$1" =~ "RESTART" ]]; then
rm -fr /scratch2/fr2424/sib/dmckeown/phex_mapping
mkdir /scratch2/fr2424/sib/dmckeown/phex_mapping
fi
