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
## $1 = input file directory (this will also be location of results) .../barcodeXX
## $2 = ONT sequencing run name
## $3 = name of specific region of interest in mapping

####################### OUTPUT FILES #####################################
## in placed in directory flye at location of input file

####################### LOAD SOFTWARE #####################################
module load minimap2/2.17
module load samtools
KRAKEN="/panfs/roc/msisoft/kraken/2.0.7beta/kraken2"
KDB="/panfs/roc/msisoft/kraken/kraken_db"
SEQKIT="/home/dcschroe/dcschroe/dmckeow/seqkit"
####################### SET VARIABLES/ARGUMENTS #####################################
###### set reference genomes, in ls order
REF2="$HOME/data/house_ref_genomes/1690666.Bombyx_mori_iflavirus.ref1.fna.gz"
REF3="$HOME/data/house_ref_genomes/1906245.Moku_virus.ref1.fna.gz"
REF4="$HOME/data/house_ref_genomes/198112.Deformed_wing_virus_genome_assembly_Devon_Type_C_scaffold_Type_C_consensus_final.ref1.fna.gz"
REF5="$HOME/data/house_ref_genomes/198112.Deformed_wing_virus_isolate_Egypt_1977_polyprotein_gene.ref1.fna.gz"
REF6="$HOME/data/house_ref_genomes/198112.Deformed_wing_virus.ref1.fna.gz"
REF7="$HOME/data/house_ref_genomes/2201278.Darwin_bee_virus_3_isolate_NT_12_polyprotein_gene.ref1.fna.gz"
REF8="$HOME/data/house_ref_genomes/232800.Varroa_destructor_virus_1.ref1.fna.gz"
REF9="$HOME/data/house_ref_genomes/458132.Slow_bee_paralysis_virus.ref1.fna.gz"
REF10="$HOME/data/house_ref_genomes/89463.Sacbrood_virus.ref1.fna.gz"

###### set region of interest in genomes
SREG2="8419-9889" ## good
SREG3="8736-9447" ## good
SREG4="8762-9623" ## good
SREG5="8770-9610" ## good
SREG6="8791-9631" ## good
SREG7="8466-9288" ## good
SREG8="8770-9605" ## good
SREG9="8155-9028" ## good
SREG10="7681-8557" ## good

################################################################################################
####################### RUNNING THE SCRIPT #####################################################
################################################################################################

OUTNAME="${2}_$(basename $1)"

#### map ALL reads against reference genomes identified by Kraken AND keep only aligned reads
rm -f $1/tmp_allreads.fastq; cat $1/*.fastq > $1/tmp_allreads.fastq
cat $REF2 $REF3 $REF4 $REF5 $REF6 $REF7 $REF8 $REF9 $REF10 > tmp_${OUTNAME}_multiple_reference_genomes.fna.gz
minimap2 -t 8 -ax map-ont tmp_${OUTNAME}_multiple_reference_genomes.fna.gz $1/tmp_allreads.fastq | samtools sort -O BAM - > ${OUTNAME}_multiple_reference_genomes.bam

$SEQKIT stat -T $1/tmp_allreads.fastq | cut -f 4 | sed '1,1d' > tmp_${OUTNAME}.all.read.count

########## label reads aligned to both  DWV A amd B with ___PRC___
##samtools sort -O SAM ${OUTNAME}_multiple_reference_genomes.bam -o ${OUTNAME}_multiple_reference_genomes.sam

##samtools view -F 4 ${OUTNAME}_multiple_reference_genomes.sam | awk '/NC_006494.1|NC_004830.2/' | cut -f 1-3 | sort -V -k 1,1 > tmp_${OUTNAME}_multiple_reference_genomes

##cut -f 1 tmp_${OUTNAME}_multiple_reference_genomes | uniq -d | grep -f - tmp_${OUTNAME}_multiple_reference_genomes | awk -F "\t" -v n=${OUTNAME}_multiple_reference_genomes.sam '{print "sed -i -E \"s/"$1"/___PRC___"$1"/g\" "n}' | sort -Vu | bash -

##samtools sort -O ${OUTNAME}_multiple_reference_genomes.sam -o ${OUTNAME}_multiple_reference_genomes.bam


############ count reads whole genome

ACC=$(zcat $REF2 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b -F 4 ${OUTNAME}_multiple_reference_genomes.bam | samtools sort -O SAM | grep -c "$ACC" > tmp_${OUTNAME}_$(basename $REF2 .ref1.fna.gz).mapped.read.count

ACC=$(zcat $REF3 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b -F 4 ${OUTNAME}_multiple_reference_genomes.bam | samtools sort -O SAM | grep -c "$ACC" > tmp_${OUTNAME}_$(basename $REF3 .ref1.fna.gz).mapped.read.count

ACC=$(zcat $REF4 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b -F 4 ${OUTNAME}_multiple_reference_genomes.bam | samtools sort -O SAM | grep -c "$ACC" > tmp_${OUTNAME}_$(basename $REF4 .ref1.fna.gz).mapped.read.count

ACC=$(zcat $REF5 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b -F 4 ${OUTNAME}_multiple_reference_genomes.bam | samtools sort -O SAM | grep -c "$ACC" > tmp_${OUTNAME}_$(basename $REF5 .ref1.fna.gz).mapped.read.count

ACC=$(zcat $REF6 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b -F 4 ${OUTNAME}_multiple_reference_genomes.bam | samtools sort -O SAM | grep -c "$ACC" > tmp_${OUTNAME}_$(basename $REF6 .ref1.fna.gz).mapped.read.count

ACC=$(zcat $REF7 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b -F 4 ${OUTNAME}_multiple_reference_genomes.bam | samtools sort -O SAM | grep -c "$ACC" > tmp_${OUTNAME}_$(basename $REF7 .ref1.fna.gz).mapped.read.count

ACC=$(zcat $REF8 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b -F 4 ${OUTNAME}_multiple_reference_genomes.bam | samtools sort -O SAM | grep -c "$ACC" > tmp_${OUTNAME}_$(basename $REF8 .ref1.fna.gz).mapped.read.count

ACC=$(zcat $REF9 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b -F 4 ${OUTNAME}_multiple_reference_genomes.bam | samtools sort -O SAM | grep -c "$ACC" > tmp_${OUTNAME}_$(basename $REF9 .ref1.fna.gz).mapped.read.count

ACC=$(zcat $REF10 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b -F 4 ${OUTNAME}_multiple_reference_genomes.bam | samtools sort -O SAM | grep -c "$ACC" > tmp_${OUTNAME}_$(basename $REF10 .ref1.fna.gz).mapped.read.count

######### index
samtools index ${OUTNAME}_multiple_reference_genomes.bam

########### count reads specific region

ACC=$(zcat $REF2 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b ${OUTNAME}_multiple_reference_genomes.bam ${ACC}:${SREG2} | samtools view -b -c - > tmp_${OUTNAME}_$(basename $REF2 .ref1.fna.gz).mapped.read.${3}.count

ACC=$(zcat $REF3 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b ${OUTNAME}_multiple_reference_genomes.bam ${ACC}:${SREG3} | samtools view -b -c - > tmp_${OUTNAME}_$(basename $REF3 .ref1.fna.gz).mapped.read.${3}.count

ACC=$(zcat $REF4 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b ${OUTNAME}_multiple_reference_genomes.bam ${ACC}:${SREG4} | samtools view -b -c - > tmp_${OUTNAME}_$(basename $REF4 .ref1.fna.gz).mapped.read.${3}.count

ACC=$(zcat $REF5 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b ${OUTNAME}_multiple_reference_genomes.bam ${ACC}:${SREG5} | samtools view -b -c - > tmp_${OUTNAME}_$(basename $REF5 .ref1.fna.gz).mapped.read.${3}.count

ACC=$(zcat $REF6 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b ${OUTNAME}_multiple_reference_genomes.bam ${ACC}:${SREG6} | samtools view -b -c - > tmp_${OUTNAME}_$(basename $REF6 .ref1.fna.gz).mapped.read.${3}.count

ACC=$(zcat $REF7 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b ${OUTNAME}_multiple_reference_genomes.bam ${ACC}:${SREG7} | samtools view -b -c - > tmp_${OUTNAME}_$(basename $REF7 .ref1.fna.gz).mapped.read.${3}.count

ACC=$(zcat $REF8 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b ${OUTNAME}_multiple_reference_genomes.bam ${ACC}:${SREG8} | samtools view -b -c - > tmp_${OUTNAME}_$(basename $REF8 .ref1.fna.gz).mapped.read.${3}.count

ACC=$(zcat $REF9 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b ${OUTNAME}_multiple_reference_genomes.bam ${ACC}:${SREG9} | samtools view -b -c - > tmp_${OUTNAME}_$(basename $REF9 .ref1.fna.gz).mapped.read.${3}.count

ACC=$(zcat $REF10 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); samtools view -b ${OUTNAME}_multiple_reference_genomes.bam ${ACC}:${SREG10} | samtools view -b -c - > tmp_${OUTNAME}_$(basename $REF10 .ref1.fna.gz).mapped.read.${3}.count



################# merge count data

ACC=$(zcat $REF2 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); paste -d "\t" tmp_${OUTNAME}.all.read.count tmp_${OUTNAME}_$(basename $REF2 .ref1.fna.gz).mapped.read.count tmp_${OUTNAME}_$(basename $REF2 .ref1.fna.gz).mapped.read.${3}.count | awk -F"\t" -v aa=$OUTNAME -v bb=$(basename $REF2 .ref1.fna.gz) -v cc=${ACC} -v dd=${3} -v ee=${SREG2} 'BEGIN{OFS="\t"};{print aa,$0,100*$2/$1,100*$3/$1,bb,cc,dd,ee}' | sed 's/_readmap_/\t/1' > tmp_${OUTNAME}_$(basename $REF2 .ref1.fna.gz).read.count.summary

ACC=$(zcat $REF3 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); paste -d "\t" tmp_${OUTNAME}.all.read.count tmp_${OUTNAME}_$(basename $REF3 .ref1.fna.gz).mapped.read.count tmp_${OUTNAME}_$(basename $REF3 .ref1.fna.gz).mapped.read.${3}.count | awk -F"\t" -v aa=$OUTNAME -v bb=$(basename $REF3 .ref1.fna.gz) -v cc=${ACC} -v dd=${3} -v ee=${SREG3} 'BEGIN{OFS="\t"};{print aa,$0,100*$2/$1,100*$3/$1,bb,cc,dd,ee}' | sed 's/_readmap_/\t/1' > tmp_${OUTNAME}_$(basename $REF3 .ref1.fna.gz).read.count.summary

ACC=$(zcat $REF4 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); paste -d "\t" tmp_${OUTNAME}.all.read.count tmp_${OUTNAME}_$(basename $REF4 .ref1.fna.gz).mapped.read.count tmp_${OUTNAME}_$(basename $REF4 .ref1.fna.gz).mapped.read.${3}.count | awk -F"\t" -v aa=$OUTNAME -v bb=$(basename $REF4 .ref1.fna.gz) -v cc=${ACC} -v dd=${3} -v ee=${SREG4} 'BEGIN{OFS="\t"};{print aa,$0,100*$2/$1,100*$3/$1,bb,cc,dd,ee}' | sed 's/_readmap_/\t/1' > tmp_${OUTNAME}_$(basename $REF4 .ref1.fna.gz).read.count.summary

ACC=$(zcat $REF5 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); paste -d "\t" tmp_${OUTNAME}.all.read.count tmp_${OUTNAME}_$(basename $REF5 .ref1.fna.gz).mapped.read.count tmp_${OUTNAME}_$(basename $REF5 .ref1.fna.gz).mapped.read.${3}.count | awk -F"\t" -v aa=$OUTNAME -v bb=$(basename $REF5 .ref1.fna.gz) -v cc=${ACC} -v dd=${3} -v ee=${SREG5} 'BEGIN{OFS="\t"};{print aa,$0,100*$2/$1,100*$3/$1,bb,cc,dd,ee}' | sed 's/_readmap_/\t/1' > tmp_${OUTNAME}_$(basename $REF5 .ref1.fna.gz).read.count.summary

ACC=$(zcat $REF6 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); paste -d "\t" tmp_${OUTNAME}.all.read.count tmp_${OUTNAME}_$(basename $REF6 .ref1.fna.gz).mapped.read.count tmp_${OUTNAME}_$(basename $REF6 .ref1.fna.gz).mapped.read.${3}.count | awk -F"\t" -v aa=$OUTNAME -v bb=$(basename $REF6 .ref1.fna.gz) -v cc=${ACC} -v dd=${3} -v ee=${SREG6} 'BEGIN{OFS="\t"};{print aa,$0,100*$2/$1,100*$3/$1,bb,cc,dd,ee}' | sed 's/_readmap_/\t/1' > tmp_${OUTNAME}_$(basename $REF6 .ref1.fna.gz).read.count.summary

ACC=$(zcat $REF7 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); paste -d "\t" tmp_${OUTNAME}.all.read.count tmp_${OUTNAME}_$(basename $REF7 .ref1.fna.gz).mapped.read.count tmp_${OUTNAME}_$(basename $REF7 .ref1.fna.gz).mapped.read.${3}.count | awk -F"\t" -v aa=$OUTNAME -v bb=$(basename $REF7 .ref1.fna.gz) -v cc=${ACC} -v dd=${3} -v ee=${SREG7} 'BEGIN{OFS="\t"};{print aa,$0,100*$2/$1,100*$3/$1,bb,cc,dd,ee}' | sed 's/_readmap_/\t/1' > tmp_${OUTNAME}_$(basename $REF7 .ref1.fna.gz).read.count.summary

ACC=$(zcat $REF8 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); paste -d "\t" tmp_${OUTNAME}.all.read.count tmp_${OUTNAME}_$(basename $REF8 .ref1.fna.gz).mapped.read.count tmp_${OUTNAME}_$(basename $REF8 .ref1.fna.gz).mapped.read.${3}.count | awk -F"\t" -v aa=$OUTNAME -v bb=$(basename $REF8 .ref1.fna.gz) -v cc=${ACC} -v dd=${3} -v ee=${SREG8} 'BEGIN{OFS="\t"};{print aa,$0,100*$2/$1,100*$3/$1,bb,cc,dd,ee}' | sed 's/_readmap_/\t/1' > tmp_${OUTNAME}_$(basename $REF8 .ref1.fna.gz).read.count.summary

ACC=$(zcat $REF9 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); paste -d "\t" tmp_${OUTNAME}.all.read.count tmp_${OUTNAME}_$(basename $REF9 .ref1.fna.gz).mapped.read.count tmp_${OUTNAME}_$(basename $REF9 .ref1.fna.gz).mapped.read.${3}.count | awk -F"\t" -v aa=$OUTNAME -v bb=$(basename $REF9 .ref1.fna.gz) -v cc=${ACC} -v dd=${3} -v ee=${SREG9} 'BEGIN{OFS="\t"};{print aa,$0,100*$2/$1,100*$3/$1,bb,cc,dd,ee}' | sed 's/_readmap_/\t/1' > tmp_${OUTNAME}_$(basename $REF9 .ref1.fna.gz).read.count.summary

ACC=$(zcat $REF10 | grep ">" | head -1 | cut -d " " -f 1 | sed 's/>//g'); paste -d "\t" tmp_${OUTNAME}.all.read.count tmp_${OUTNAME}_$(basename $REF10 .ref1.fna.gz).mapped.read.count tmp_${OUTNAME}_$(basename $REF10 .ref1.fna.gz).mapped.read.${3}.count | awk -F"\t" -v aa=$OUTNAME -v bb=$(basename $REF10 .ref1.fna.gz) -v cc=${ACC} -v dd=${3} -v ee=${SREG10} 'BEGIN{OFS="\t"};{print aa,$0,100*$2/$1,100*$3/$1,bb,cc,dd,ee}' | sed 's/_readmap_/\t/1' > tmp_${OUTNAME}_$(basename $REF10 .ref1.fna.gz).read.count.summary


#### when all done do:
#find . -name "tmp_*.read.count.summary" > tmp
#for f in $(cat tmp); do cat $f ; done > allbarcodes.read.count.summary
#sed -i -z 's/^/ONT_run\tONT_barcode\ttotal_reads_all\ttotal_reads_mapped_to_genome\ttotal_reads_mapped_to_specific_region\tpercent_reads_mapped_to_genome\tpercent_reads_mapped_to_specific_region\treference_genome\treference_accession\tspecific_region\tspecific_region_coords\n/g' allbarcodes.read.count.summary

#for f in *genomes.bam; do samtools view -F 4 $f | awk '/NC_006494.1|NC_004830.2/' | cut -f 1-3 | sort -V -k 1,1 > tmp_$(basename $f .bam); samtools sort -O SAM $f -o $(basename $f .bam).sam; cut -f 1 tmp_$(basename $f .bam) | uniq -d | grep -f - tmp_$(basename $f .bam) | awk -F "\t" -v n=$(basename $f .bam).sam '{print "sed -i -E \"s/"$1"/___PRC___"$1"/g\" "n}' | sort -Vu | bash - ; samtools sort -O BAM $(basename $f .bam).sam -o $f; done
