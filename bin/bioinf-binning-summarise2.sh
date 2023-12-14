#!/bin/bash

umask 002
green='\033[32m';red='\033[31m';cyan='\033[36m';nocolor='\033[m'
echo -e "${green}RUNNING bioinf-binning-summarise with the following parameters: $@ ${nocolor}"

Help()
{
echo -e "########### This script is for generating various outputs for checking contigs within a bin. You must have ran anvio and anvi-summarise for this script to work. This script is designed to run on a single bin per job submission and within the folder for each bin. RUN this script within the BINNING_ output folder of your binning run ###########\n"
echo -e "-r -reads [REQUIRED]\tFile listing the full paths to read fastq files in single column, to map against contigs"
echo -e "-c -contigs [REQUIRED]\tFile listing the full paths to contig fasta files in single column, keep separate by sample"
echo -e "-s -step\tWhich script steps to run"
echo -e "-t -threads\tNumber of parallel processes"

}

### set euclidean options to default setting here

##

while getopts r:c:i:s:t:h option
do 
    case "${option}" in
    r)reads=${OPTARG};;
    c)contigs=${OPTARG};;
    i)input=${OPTARG};;
    s)step=${OPTARG};;
    t)threads=${OPTARG};;
    h)Help; exit;;
    esac

done

if [[ -z "${reads}" ]]; then echo "-r, -reads REQUIRED"; Help; exit; fi
if [[ -z "${contigs}" ]]; then echo "-c, -contigs REQUIRED"; Help; exit; fi


### THREADS

if [[ ! -z "${SLURM_CPUS_PER_TASK}" ]]; then THREADS="${SLURM_CPUS_PER_TASK}"; echo -e "${green}Using threads set by SLURM_CPUS_PER_TASK:${nocolor}"; fi

if [[ ! -z "${threads}" ]]; then THREADS="${threads}"; echo -e "${green}Using threads set by -threads flag:${nocolor}"; fi

if [[ -z "${threads}" ]] && [[ -z "${SLURM_CPUS_PER_TASK}" ]]; then THREADS="1"; echo -e "${green}No SLURM_CPUS_PER_TASK or -threads set, using default threads:${nocolor}"; fi
echo -e "\t${THREADS} threads"


####################### PREREQUISITES #####################################

####################### SCRIPT PURPOSE #####################################

####################### HOW TO RUN SCRIPT #####################################

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################

####################### OUTPUT FILES #####################################

####################### SOFTWARE, DATABSES, ETC #####################################
### LOAD software available via shared environment on server:
eval "$(conda shell.bash hook)"
conda activate bioinftools

########################################################################
########################################################################

tmp="tmp.bioinf_binning_summarise"
PROJECT_NAME=$(basename $contigs)
REFFASTA="/panfs/jay/groups/27/dcschroe/dmckeow/data/house_ref_genomes/Spillover_FINAL_refs_formap.fa"

####################### THE SCRIPT #####################################

#################################################################

########################################## STEP A1 - get bin info for anvio run (repeat for each run you want to include frmo within main folder - concatenate them into a single "bin_list" and run rest of script using that)
if [[ -z "${step}" ]] || [[ "$step" =~ "A1" ]]; then
###########################################
rm -f bin_list; touch bin_list
if [[ -z "${input}" ]]; then echo "-i, -input REQUIRED"; Help; exit; fi


for f in $(find */bin_by_bin/* -name "*-contigs.fa"); do
    B=$(echo $(basename $PWD)"___"$(basename $f) | sed -e 's/-contigs\.fa//g' -e 's/-SUMMARY___/___/g')
    grep ">" $f | sed "s/>/$B\t/g"
done >> bin_list

sed -e 's/;\|\t/,/g' -e 's/\//;/g' ${input} | awk -F"," '{print "s/^"$3"$/"$0"/g"}' > tmp.fixer.input
cut -f 2 bin_list | sed -E 's/contig_(.+)_[0-9]+$/\1/g' | sed -f tmp.fixer.input - | paste bin_list - | sed -e 's/,/\t/g' -e 's/;/\//g' > tmp && mv tmp bin_list

sed -i 's/;/,/g' bin_list

awk -F "\t" '{print "s/^"$1"$/"$2"/g"}' *.deflinekey > tmp.fixer

cut -f 2 bin_list | sed -f tmp.fixer - | paste bin_list - > tmp && mv tmp bin_list

sort -t $'\t' -k 5,5V bin_list -o bin_list


###########################################
fi
###########################################

########################################## STEP A2 - prepare mapping plots for the contigs within the bin
if [[ -z "${step}" ]] || [[ "$step" =~ "A2" ]]; then
###########################################

#### dereplicate contigs before mapping
rm -fr tmp.derep; mkdir tmp.derep
for f in $(cat ${contigs}); do
    seqkit seq --min-len 2000 $f > tmp.derep/$(basename $f).0
    cd-hit-est -c 0.95 -s 0.8 -n 4 -i tmp.derep/$(basename $f).0 -o tmp.derep/$(basename $f)
    rm -f tmp.derep/$(basename $f).0
done

rm -f derep-contigs.fa; touch derep-contigs.fa
rm -f tmp.derep/*.clstr

for f in tmp.derep/*; do
    N=$(basename $f)
    sed -i "s/>/>${N}______/g" $f
done

cat tmp.derep/* >> derep-contigs.fa

#### ALL contigs vs ALL reads for each sample in the bin and count them per contig
rm -fr tmp.mapping; mkdir tmp.mapping

    #### prepare file to run mapping from
    for f in $(cat ${reads}); do
        Q="sbatch --time=12:00:00 --cpus-per-task=4 --mem=48GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval \"\$(conda shell.bash hook)\"; conda activate bioinftools;"
        R="'"
        echo "$Q seqkit seq --min-len 200 ${f} | minimap2 -t \$SLURM_CPUS_PER_TASK -ax map-ont --secondary=no derep-contigs.fa - | samtools view -h -F 0x900 | samtools sort -O BAM - > tmp.mapping/$(basename $f).bam; samtools view -F 4 tmp.mapping/$(basename $f).bam | cut -f 1,3 | sort -Vu | cut -f 2 | sort -V | uniq -c | sed -E \"s/ +([0-9]+) (.+)/\1\t\2/g\" > tmp.mapping/$(basename $f).bam.reads_mapped$R"

    done > run_mapping

bash run_mapping

###########################################
fi
###########################################


################################################################################################
if [[ -z "${step}" ]] || [[ "$step" =~ "A3" ]]; then
################################################################################################
######### run kaiju on the bctrimmed reads (runs separately from any other steps, but must have completed Step A1 )

for f in "$bioinfdb"/KAIJU/*; do
    kaiju -t ${f}/nodes.dmp -f ${f}/*.fmi -i derep-contigs.fa -o derep-contigs.fa.kaiju.$(basename $f) -z $THREADS -v
    sort -t $'\t' -V -k 2,2 derep-contigs.fa.kaiju.$(basename $f) -o derep-contigs.fa.kaiju.$(basename $f)
    kaiju-addTaxonNames -t ${f}/nodes.dmp -n ${f}/names.dmp -i derep-contigs.fa.kaiju.$(basename $f) -o derep-contigs.fa.kaiju.$(basename $f).names -r superkingdom,phylum,order,class,family,genus,species
done

awk '{print $0"\t0\tNA\tNA\tNA\tNA"}' derep-contigs.fa.kaiju.*.names | cut -f 1-8 | sort -t $'\t' -k 2,2V -k 4,4nr | awk -F "\t" '!a[$2]++' > derep-contigs.fa.kaiju.merge

################################################################################################
fi
################################################################################################

########################################## STEP A4 
if [[ -z "${step}" ]] || [[ "$step" =~ "A4" ]]; then
###########################################
   
grep ">" derep-contigs.fa | sed -e 's/>//g' -e 's/ .*//g' > tmp.list
    
    for f in tmp.mapping/*.bam.reads_mapped; do
        N=$(echo $f | sed -e 's/\.bam\.reads_mapped//g' -e 's/.*\///g')
        cut -f 2 $f | grep -vwF -f - tmp.list | sed 's/^/0\t/g' > ${f}_2
        cat $f >> ${f}_2
        sed -i "s/^/$N\t/g" ${f}_2
        sort -t $'\t' -k 3,3V ${f}_2 -o ${f}_2
    done

rm -f ALL.bam.reads_mapped.ALL; touch ALL.bam.reads_mapped.ALL
cat tmp.mapping/*.bam.reads_mapped_2 >> ALL.bam.reads_mapped.ALL

#### get reads mapped total per contig and by sample (that reads came from)

awk -F'\t' '{C = $3; n=$2; sum[C] += n;} END{for (C in sum)print C"\t"sum[C];}' ALL.bam.reads_mapped.ALL | awk -F "\t" '{print "s/^"$1"$/"$1"\\t"$2"/g"}' > reads_mapped_bycontig
cut -f 3 ALL.bam.reads_mapped.ALL | sed -f reads_mapped_bycontig - | cut -f 2 | paste ALL.bam.reads_mapped.ALL - > tmp && mv tmp ALL.bam.reads_mapped.ALL

awk -F'\t' '{C = $1; n=$2; sum[C] += n;} END{for (C in sum)print C"\t"sum[C];}' ALL.bam.reads_mapped.ALL | awk -F "\t" '{print "s/^"$1"$/"$1"\\t"$2"/g"}' > reads_mapped_bysample
cut -f 1 ALL.bam.reads_mapped.ALL | sed -f reads_mapped_bysample - | cut -f 2 | paste ALL.bam.reads_mapped.ALL - > tmp && mv tmp ALL.bam.reads_mapped.ALL

### get total of reads in each sample
rm -f read_total_per_sample; touch read_total_per_sample

for f in $(cat ${reads}); do
    seqkit stat -j $SLURM_CPUS_PER_TASK -T $f | cut -f 1,4 | sed "/num_seqs/d" >> read_total_per_sample
done

sed 's/.*\///g' read_total_per_sample | awk -F "\t" '{print "s/^"$1"$/"$1"\\t"$2"/g"}' > tmp && mv tmp read_total_per_sample

cut -f 1 ALL.bam.reads_mapped.ALL | sed -f read_total_per_sample - | cut -f 2 | paste ALL.bam.reads_mapped.ALL - > tmp && mv tmp ALL.bam.reads_mapped.ALL

#### get contig lengths

samtools faidx derep-contigs.fa
cut -f 1,2 derep-contigs.fa.fai | awk -F "\t" '{print "s/^"$1"$/"$1"\\t"$2"/g"}' > contig_length

cut -f 3 ALL.bam.reads_mapped.ALL | sed -f contig_length - | cut -f 2 | paste ALL.bam.reads_mapped.ALL - > tmp && mv tmp ALL.bam.reads_mapped.ALL

###########################################
fi
###########################################

########################################## STEP A5 
if [[ -z "${step}" ]] || [[ "$step" =~ "A5" ]]; then
###########################################

########## add kaiju info
awk -F "\t" '{print "s/^"$2"$/"$0"/g"}' derep-contigs.fa.kaiju.merge | sed 's/\t/\\t/g' > kaiju_fixer

cut -f 3 ALL.bam.reads_mapped.ALL | sed -f kaiju_fixer - | paste ALL.bam.reads_mapped.ALL - > tmp && mv tmp ALL.bam.reads_mapped.ALL

##### add bin_list if it exists
if [ -f bin_list ]; then
  awk -F "\t" '{print $3"______"$6"\t"$1}' bin_list | sed -e 's/ .*\t/\t/g' -e 's/.*\///g' | awk -F "\t" '{print "s/^"$1"$/"$2"/g"}' > tmp.binfixer

  ##### set to replace contig names with no bin as NA
  cut -f 3 ALL.bam.reads_mapped.ALL | sort -Vu > tmp.list
  awk -F "\t" '{print $3"______"$6"\t"$1}' bin_list | sed -e 's/ .*\t/\t/g' -e 's/.*\///g' | cut -f 1 | sort -Vu | cat - tmp.list | sort -V | uniq -u | awk -F "\t" '{print "s/^"$1"$/NA/g"}' >> tmp.binfixer

  cut -f 3 ALL.bam.reads_mapped.ALL | sed -f tmp.binfixer - | paste ALL.bam.reads_mapped.ALL - > tmp && mv tmp ALL.bam.reads_mapped.ALL
fi

#### if no binning info, jsut enter NA
if ! [ -f bin_list ]; then
    sed -i 's/$/\tNA/g' ALL.bam.reads_mapped.ALL
fi

#### headers
sed -i -z "s/^/\
specific_read_source\t\
num_specific_read_source_reads_mapped_to_specific_contig\t\
specific_contig\t\
total_num_reads_mapped_to_specific_contig_from_all_read_sources\t\
total_num_reads_mapped_to_all_contigs_from_specific_read_source\t\
num_reads_in_specific_read_source\t\
specific_contig_length\t\
kaiju_C_or_U\t\
kaiju_specific_contig\t\
kaiju_taxid\t\
kaiju_score\t\
kaiju_taxids_all\t\
kaiju_accessions\t\
kaiju_aligned_seq\t\
kaiju_taxonomy\t\
bin_of_specific_contig\t\
\n/g" ALL.bam.reads_mapped.ALL


#####################################
fi
#####################################

########################################## STEP A6
if [[ -z "${step}" ]] || [[ "$step" =~ "A6" ]]; then
###########################################

########### breadth of coverage

rm -f ALL.coverage.stats.ALL; touch ALL.coverage.stats.ALL

for f in tmp.mapping/*.bam; do
N=$(basename $f)
samtools index $f

## get contig lengths
samtools idxstats $f | cut -f 1,2 | sed '/^\*/d' > ${f}.idxstats

## total coverage per base for each contig
samtools depth $f | awk -F "\t" -v N=$N '{print N"\t"$0}' > ${f}.samtoolsdepth

done

## contig length fix
rm -f tmp.ALLcontig; touch tmp.ALLcontig
cat tmp.mapping/*.idxstats >> tmp.ALLcontig
sort -Vu tmp.ALLcontig -o tmp.ALLcontig

## get total number of bases covered at MIN_COVERAGE_DEPTH or higher PER contig
rm -f tmp.ALLstats.list.count0 tmp.ALLstats.list.count1 tmp.ALLstats.list.count2 tmp.ALLstats.list.count3; touch tmp.ALLstats.list.count0 tmp.ALLstats.list.count1 tmp.ALLstats.list.count2 tmp.ALLstats.list.count3
for f in tmp.mapping/*.samtoolsdepth; do
    cut -f 1,2 $f | sort -Vu >> tmp.ALLstats.list.count0
    awk '$4 >= 5' $f | cut -f 1,2 | sort -V | uniq -c | sed -E 's/ +(.+) (.+)/\2\t\1/g' >> tmp.ALLstats.list.count1
    awk '$4 >= 10' $f | cut -f 1,2 | sort -V | uniq -c | sed -E 's/ +(.+) (.+)/\2\t\1/g' >> tmp.ALLstats.list.count2
    awk '$4 >= 100' $f | cut -f 1,2 | sort -V | uniq -c | sed -E 's/ +(.+) (.+)/\2\t\1/g' >> tmp.ALLstats.list.count3
done

### 1X coverage
rm -f tmp.ALLstats.list.count0.5; touch tmp.ALLstats.list.count0.5
for f in tmp.mapping/*.samtoolsdepth; do
    cut -f 1,2 $f | sort -V | uniq -c | sed -E 's/ +(.+) (.+)/\2\t\1/g' >> tmp.ALLstats.list.count0.5
done

## add headers
sed -i 's/\.bam\t/\t/g' tmp.ALLcontig
sed -i 's/\.bam\t/\t/g' tmp.ALLstats.list.count*

sed -i -z 's/^/contig\tcontig_length\n/1' tmp.ALLcontig
sed -i -z 's/^/reads\tcontig\n/1' tmp.ALLstats.list.count0
sed -i -z 's/^/reads\tcontig\ttotal_bases_over_1_coverage\n/1' tmp.ALLstats.list.count0.5
sed -i -z 's/^/reads\tcontig\ttotal_bases_over_5_coverage\n/1' tmp.ALLstats.list.count1
sed -i -z 's/^/reads\tcontig\ttotal_bases_over_10_coverage\n/1' tmp.ALLstats.list.count2
sed -i -z 's/^/reads\tcontig\ttotal_bases_over_100_coverage\n/1' tmp.ALLstats.list.count3


#####################################
fi
#####################################



########################################## STEP A7
if [[ -z "${step}" ]] || [[ "$step" =~ "A7" ]]; then
###########################################

## map reads of each sample against reference genomes , to later identify the coverage of regions of interest - i.e. the RDRP primer target

#### ALL contigs vs ALL reads for each sample in the bin and count them per contig
rm -fr tmp.mapping.ref; mkdir tmp.mapping.ref

    #### prepare file to run mapping from
    for f in $(cat ${reads}); do
        Q="sbatch --time=12:00:00 --cpus-per-task=4 --mem=48GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval \"\$(conda shell.bash hook)\"; conda activate bioinftools;"
        R="'"
        echo "$Q seqkit seq --min-len 200 ${f} | minimap2 -t \$SLURM_CPUS_PER_TASK -ax map-ont --secondary=no $REFFASTA - | samtools view -h -F 0x900 | samtools sort -O BAM - > tmp.mapping.ref/$(basename $f).bam$R"

    done > run_mapping_vs_ref

bash run_mapping_vs_ref

#####################################
fi
#####################################

########################################## STEP A8
if [[ -z "${step}" ]] || [[ "$step" =~ "A8" ]]; then
###########################################

################ get the intervals by window size for sequences that are targets for mapping
samtools faidx $REFFASTA
cut -f 1,2 ${REFFASTA}.fai > tmp.reffasta
WINDOW=100

while IFS=$'\t' read -r sequence length; do
    for ((start=1; start<=length; start+=WINDOW)); do
        end=$((start+WINDOW-1))
        if [ "$end" -gt "$length" ]; then
            end="$length"
        fi
        echo -e "$sequence\t$start\t$end"
    done
done < tmp.reffasta > $(basename ${REFFASTA}).bed



rm -f tmp.mapping.ref/*.meanwindowdepth
rm -f tmp.mapping.ref/*.idxstats
##rm -f tmp.mapping.ref/*.samtoolsdepth

for f in tmp.mapping.ref/*.bam; do
    N=$(basename $f)
    bedtools coverage -mean -a $(basename ${REFFASTA}).bed -b $f | awk -F "\t" -v N=$N '{print N"\t"$0}' > ${f}.meanwindowdepth

    ##### get the contig length and total reads mapped for each reference sequence 

    samtools index $f
    samtools idxstats $f | sed '/^\*/d' | cut -f 1,2 | awk -F "\t" -v N=$N '{print N"\t"$0}' > ${f}.idxstats
    NR=$(samtools view -c $f)
    sed -i "s/$/\t$NR/g" ${f}.idxstats ## add the true num reads in library
    ###### samtools depth replaced by bedtools, but might be useful for checking things
    ###samtools depth -a $f | awk -F "\t" -v N=$N '{print N"\t"$0}' > ${f}.samtoolsdepth


done

rm -f ALL.mapping.ref.idxstats; touch ALL.mapping.ref.idxstats

for f in tmp.mapping.ref/*.idxstats; do
    cat $f >> ALL.mapping.ref.idxstats
done

rm -f ALL.mapping.ref.meanwindowdepth; touch ALL.mapping.ref.meanwindowdepth

for f in tmp.mapping.ref/*.meanwindowdepth; do
    cat $f >> ALL.mapping.ref.meanwindowdepth
done

### add headers
sed -i -z 's/^/specific_read_source\treference\tref_length\tNumReadsInSpecificReadSource\n/g' ALL.mapping.ref.idxstats
sed -i -z 's/^/specific_read_source\treference\tref_start\tref_end\tavg_coverage_by_window\n/g' ALL.mapping.ref.meanwindowdepth

##### get the contig length and total reads mapped for each reference sequence 




#### get total number of bases covered at MIN_COVERAGE_DEPTH or higher PER contig
###rm -f tmp.ALLstats.list.count0 tmp.ALLstats.list.count1 tmp.ALLstats.list.count2 tmp.ALLstats.list.count3; touch tmp.ALLstats.list.count0 tmp.ALLstats.list.count1 tmp.ALLstats.list.count2 tmp.ALLstats.list.count3
###for f in tmp.mapping/*.samtoolsdepth; do
## #   cut -f 1,2 $f | sort -Vu >> tmp.ALLstats.list.count0
##  #  awk '$4 >= 5' $f | cut -f 1,2 | sort -V | uniq -c | sed -E 's/ +(.+) (.+)/\2\t\1/g' >> tmp.ALLstats.list.count1
##   # awk '$4 >= 10' $f | cut -f 1,2 | sort -V | uniq -c | sed -E 's/ +(.+) (.+)/\2\t\1/g' >> tmp.ALLstats.list.count2
##    #awk '$4 >= 100' $f | cut -f 1,2 | sort -V | uniq -c | sed -E 's/ +(.+) (.+)/\2\t\1/g' >> tmp.ALLstats.list.count3
###done
##
#### add headers
###sed -i 's/\.bam\t/\t/g' tmp.ALLcontig
###sed -i 's/\.bam\t/\t/g' tmp.ALLstats.list.count*
##
###sed -i -z 's/^/contig\tcontig_length\n/1' tmp.ALLcontig
###sed -i -z 's/^/reads\tcontig\n/1' tmp.ALLstats.list.count0
###sed -i -z 's/^/reads\tcontig\ttotal_bases_over_1_coverage\n/1' tmp.ALLstats.list.count0.5
###sed -i -z 's/^/reads\tcontig\ttotal_bases_over_5_coverage\n/1' tmp.ALLstats.list.count1
###sed -i -z 's/^/reads\tcontig\ttotal_bases_over_10_coverage\n/1' tmp.ALLstats.list.count2
###sed -i -z 's/^/reads\tcontig\ttotal_bases_over_100_coverage\n/1' tmp.ALLstats.list.count3


#####################################
fi
#####################################

########################################## STEP A9
if [[ -z "${step}" ]] || [[ "$step" =~ "A9" ]]; then
###########################################

####MAXTARGETS=$(seqkit stat -T ${project}-SPLITS.fa | sed '1,1d' | awk -F "\t" '{printf "%0.0f\n", $8/100}')
##diamond blastx --max-target-seqs $MAXTARGETS --evalue 1e-20 --sensitive -p $THREADS -d $f -q ${project}-SPLITS.fa -f 6 -o $(basename $f .dmnd).evalue_1e-20.dmnd.blastx

DmndDb1="/panfs/jay/groups/27/dcschroe/shared/bioinfdb/DMND/vog.dmnd"

diamond blastx --range-culling --top 20 -F 15 --evalue 1e-20 --sensitive -p $THREADS -d $DmndDb1 -q $REFFASTA -f 6 -o $(basename $REFFASTA).dmnd.blastx


### best hit per split subject combo
awk -F "\t" '!a[$1,$2]++' $(basename $REFFASTA).dmnd.blastx > $(basename $REFFASTA).dmnd.blastx.best


####################### BLAST #############################
### IF any blast database files are found in bioinfdb/BLAST, then also do a BLASTx search against those

BlastDb1="/common/bioref/blast/latest/nr"
blastx -evalue 1e-20 -num_threads $THREADS -db $BlastDb1 -query $REFFASTA -outfmt 6 -out $(basename $REFFASTA).ncbi.blastx

### best hit per split subject combo
 awk -F "\t" '!a[$1,$2]++' $(basename $REFFASTA).ncbi.blastx > $(basename $REFFASTA).ncbi.blastx.best

##############################################################
fi
##############################################################
