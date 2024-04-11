#!/bin/bash

umask 002
green='\033[32m';red='\033[31m';cyan='\033[36m';nocolor='\033[m'
echo -e "${green}RUNNING bioinf-binning-summarise with the following parameters: $@ ${nocolor}"

Help()
{
echo -e "########### This script is for generating various outputs for checking contigs within a bin. You must have ran anvio and anvi-summarise for this script to work. This script is designed to run on a single bin per job submission and within the folder for each bin. RUN this script within the BINNING_ output folder of your binning run ###########\n"
echo -e "-i -input\tsame input file used for running anvio binning script - three columns, ; delimited, full/path/to/contig;/full/path/to/reads;sample_name"
echo -e "-s -step [default is run steps A, B, C]\tOptionally, you can specify which script steps to run. This script has multiple steps: A, B, C, and D which must complete before the next step is run. Step D si completely optional and requires a file to be manually edited before running step D. "
echo -e "-r -reference \tOptional - path to a nucleotide fasta to map reads against. Useful for external reference genomes, etc. Can provide multiple paths separated by commas. Sequences within a single file will be subject to competitive mapping, whereas each separate reference fasta will be a separate mapping "
echo -e "-t -threads\tOptional. Number of parallel processes"
}

### set euclidean options to default setting here

##

while getopts r:i:s:t:h option
do 
    case "${option}" in
    i)input=${OPTARG};;
    s)step=${OPTARG};;
    t)threads=${OPTARG};;
    r)reference=${OPTARG};;
    h)Help; exit;;
    esac

done

if [[ -z "${input}" ]]; then echo "-i, -input REQUIRED"; Help; exit; fi


### THREADS

if [[ ! -z "${SLURM_CPUS_PER_TASK}" ]]; then THREADS="${SLURM_CPUS_PER_TASK}"; echo -e "${green}Using threads set by SLURM_CPUS_PER_TASK:${nocolor}"; fi

if [[ ! -z "${threads}" ]]; then THREADS="${threads}"; echo -e "${green}Using threads set by -threads flag:${nocolor}"; fi

if [[ -z "${threads}" ]] && [[ -z "${SLURM_CPUS_PER_TASK}" ]]; then THREADS="1"; echo -e "${green}No SLURM_CPUS_PER_TASK or -threads set, using default threads:${nocolor}"; fi
echo -e "\t${THREADS} threads"

if [[ -z "${step}" ]]; then step="A_B_C"; fi

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
PROJECT_NAME=$(basename $input)

####################### THE SCRIPT #####################################

#################################################################




########################################## STEP A2 - prepare mapping plots for the contigs
if [[ "${step}" =~ "A" ]] || [[ "$step" =~ "a2" ]]; then
###########################################

#### dereplicate contigs before mapping
rm -fr tmp.derep; mkdir tmp.derep
rm -f all-contigs.fa; touch all-contigs.fa

for f in $(cat ${input}); do
    CONTIGS=$(echo $f | cut -d ";" -f 1)
    seqkit seq --min-len 1000 $CONTIGS > tmp.derep/$(basename $CONTIGS).0
    cd-hit-est -c 0.9 -s 0.8 -n 4 -i tmp.derep/$(basename $CONTIGS).0 -o tmp.derep/$(basename $CONTIGS)
    ####### get all the contigs without any filtering for all other tools
    cat $CONTIGS > tmp.derep/$(basename $CONTIGS).0
done

rm -f derep-contigs.fa; touch derep-contigs.fa

for f in tmp.derep/*; do
    N=$(basename $f)
    sed -i "s/>/>${N}______/g" $f
done

for f in tmp.derep/*.0; do
    N=$(basename $f)
    sed -i "s/\.0______/______/g" $f
done

cat tmp.derep/*.0 > all-contigs.fa

#################################################
#########################################################

if [ -f input_merge_bc ]; then

#### dereplicate further by specific merged barcodes, using last contig in each row as merged file name
rm -fr tmp.derep.merge; mkdir tmp.derep.merge
sed -E -e 's/^(.*)\t(.*)$/\1\t\2.contigs.fasta > tmp \&\& mv tmp tmp.derep.merge\/\2/g' -e 's/^/cat tmp.derep\//g' -e 's/$/.contigs.fasta/g' -e 's/\t/.contigs.fasta tmp.derep\//g' input_merge_bc > tmp.run_bc_merge
bash tmp.run_bc_merge

for f in tmp.derep.merge/*.contigs.fasta; do
    cd-hit-est -c 0.95 -s 0.9 -n 4 -i $f -o ${f}.derep
    mv ${f}.derep $f
done

sed -E 's/^cat (.*) >.*/rm -f \1/g' tmp.run_bc_merge > tmp.run_bc_merge.rm
bash tmp.run_bc_merge.rm
for f in tmp.derep.merge/*.contigs.fasta; do mv $f tmp.derep/ ; done

fi

#################################################
#################################################

cat tmp.derep/*.fasta >> derep-contigs.fa

#### ALL contigs vs ALL reads for each sample in the bin and count them per contig
rm -fr tmp.mapping; mkdir tmp.mapping

    #### prepare file to run mapping from
    for f in $(cat ${input}); do
        READS=$(echo $f | cut -d ";" -f 2)
        Q="sbatch --time=12:00:00 --cpus-per-task=4 --mem=48GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval \"\$(conda shell.bash hook)\"; conda activate bioinftools;"
        R="'"
        echo "$Q seqkit seq --min-len 200 ${READS} | minimap2 -t \$SLURM_CPUS_PER_TASK -ax map-ont --secondary=no derep-contigs.fa - | samtools view -h -F 2308 | samtools sort -O BAM - > tmp.mapping/$(basename $READS).bam; samtools view -F 4 tmp.mapping/$(basename $READS).bam | cut -f 1,3 | sort -Vu | cut -f 2 | sort -V | uniq -c | sed -E \"s/ +([0-9]+) (.+)/\1\t\2/g\" > tmp.mapping/$(basename $READS).bam.reads_mapped$R"

    done > run_mapping

bash run_mapping

##### record what contigs removed by dereplication
grep ">" derep-contigs.fa | sed -e 's/ .*//g' -e 's/>//g' > tmp.deflines.derep.contigs
grep ">" all-contigs.fa | sed -e 's/ .*//g' -e 's/>//g' | grep -vwF -f tmp.deflines.derep.contigs - > tmp.deflines.all-contigs

sed -i 's/$/\tNO/g' tmp.deflines.derep.contigs
sed -i 's/$/\tYES/g' tmp.deflines.all-contigs
cat tmp.deflines.derep.contigs tmp.deflines.all-contigs | sort -V -t $'\t' -k 1,1 > tmp.deflines.all-contigs.1

seqkit seq --max-len 999 all-contigs.fa | grep ">" | sed -e 's/ .*//g' -e 's/$/\tYES/g' -e 's/>//g' > tmp.deflines.all-contigs.2
seqkit seq --min-len 1000 all-contigs.fa | grep ">" | sed -e 's/ .*//g' -e 's/$/\tNO/g' -e 's/>//g' >> tmp.deflines.all-contigs.2
sort -V -t $'\t' -k 1,1 -o tmp.deflines.all-contigs.2 tmp.deflines.all-contigs.2
paste tmp.deflines.all-contigs.1 tmp.deflines.all-contigs.2 | awk -F "\t" '{if($4 == "YES") print $1"\tNO\t"$4; else print $1"\t"$2"\t"$4}' | sed -z 's/^/contig\tWasContigRemovedByDerep\tWasContigRemovedByLengthCutoff\n/g' > SeqsRemovedByClustering


## mmseqs cluster the derep contigs

############# cluster at 90 % sequence ID over 80 % of length (of the shorter sequence in pairwise comparison) [Zayed et al. 2022, Science 376, 156-162]
mmseqs easy-cluster all-contigs.fa binned_contigs-mmseqs tmp.binned_contigs-mmseqs --min-seq-id 0.9 -c 0.8 --cov-mode 5

### give mmseqs clusters group names
awk -F "\t" '{if($1 == $2) print $1}' binned_contigs-mmseqs_cluster.tsv > bin_mmseq
counter=1
while IFS=$'\t' read -r line; do
    group=$(printf '%03d' $counter)
    echo "/^$line\\t/ s/$/\t$group/g"
    ((counter++))
done < bin_mmseq | sed -f - binned_contigs-mmseqs_cluster.tsv > tmp && mv tmp bin_mmseq

cut -f 3 bin_mmseq | sort -V | uniq -c | sed -E 's/^ +(.+) (.+)/\/\\t\2\$\/ s\/$\/\\t\1\/g/g' | sed -f - bin_mmseq | cut -f 2-4 > tmp_bin_mmseq && mv tmp_bin_mmseq bin_mmseq

sed -i -z 's/^/Contig\tMmseqCluster\tNumInMmseqCluster\n/g' bin_mmseq

for f in ${bioinfdb}/DMND/*; do
    diamond blastx --range-culling --max-target-seqs 1 -F 15 --evalue 1e-20 --sensitive -p $THREADS -d $f -q all-contigs.fa -f 6 -o $(basename $f .dmnd)_AllContigs.dmnd.blastx
done


for f in ${bioinfdb}/DMND/*; do
    awk -F "\t" '!a[$1]++' $(basename $f .dmnd)_AllContigs.dmnd.blastx > tmp.$(basename $f .dmnd)_AllContigs.dmnd.blastx && mv tmp.$(basename $f .dmnd)_AllContigs.dmnd.blastx $(basename $f .dmnd)_AllContigs.dmnd.blastx
    sed -i -z 's/^/qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n/g' $(basename $f .dmnd)_AllContigs.dmnd.blastx
done



samtools faidx all-contigs.fa; cut -f 1,2 all-contigs.fa.fai > AllContigLengths.tsv
sed -i -z 's/^/Contig\tLength\n/g' AllContigLengths.tsv

###########################################
fi
###########################################

################################################################################################ B1
if [[ "${step}" =~ "B" ]] || [[ "$step" =~ "b1" ]]; then
################################################################################################
######### run kaiju on the bctrimmed reads (runs separately from any other steps, but must have completed Step A1 )

for f in "$bioinfdb"/KAIJU/*; do
    kaiju -t ${f}/nodes.dmp -f ${f}/*.fmi -i all-contigs.fa -o all-contigs.fa.kaiju.$(basename $f) -z $THREADS -v
    sort -t $'\t' -V -k 2,2 all-contigs.fa.kaiju.$(basename $f) -o all-contigs.fa.kaiju.$(basename $f)
    kaiju-addTaxonNames -t ${f}/nodes.dmp -n ${f}/names.dmp -i all-contigs.fa.kaiju.$(basename $f) -o all-contigs.fa.kaiju.$(basename $f).names -r superkingdom,phylum,order,class,family,genus,species
    awk '{print $0"\t0\tNA\tNA\tNA\tNA"}' all-contigs.fa.kaiju.$(basename $f).names | cut -f 1-8 | sort -t $'\t' -k 2,2V -k 4,4nr | awk -F "\t" '!a[$2]++' | sed 's/0\tNA\tNA\tNA\tNA//g' > tmp.all-contigs.fa.kaiju.$(basename $f).names && mv tmp.all-contigs.fa.kaiju.$(basename $f).names all-contigs.fa.kaiju.$(basename $f).names
    sed -i -z 's/^/kaiju_C_or_U\tkaiju_specific_contig\tkaiju_taxid\tkaiju_score\tkaiju_taxids_all\tkaiju_accessions\tkaiju_aligned_seq\tkaiju_taxonomy\n/g' all-contigs.fa.kaiju.$(basename $f).names
done



#ISJOBRUNNING=$(squeue --me | grep -s "wrap" | wc -l); 
#until [[ $ISJOBRUNNING -le 0 ]]; do 
 #   sleep 1m; ISJOBRUNNING=$(squeue --me | grep -s "wrap" | wc -l);
#done;

################################################################################################
fi
################################################################################################

########################################## STEP B2
if [[ "${step}" =~ "B" ]] || [[ "$step" =~ "b2" ]]; then
###########################################
       
### get total of reads in each sample
rm -f read_total_per_sample; touch read_total_per_sample

for f in $(cat ${input}); do
    READS=$(echo $f | cut -d ";" -f 2)
    seqkit stat -j $SLURM_CPUS_PER_TASK -T $READS | cut -f 1,4 | sed "/num_seqs/d" >> read_total_per_sample
done

sed 's/.*\///g' read_total_per_sample | awk -F "\t" '{print $1"\t"$2}' > tmp && mv tmp read_total_per_sample

sed -i -z 's/^/specific_read_source\tnum_reads_in_specific_read_source\n/g' read_total_per_sample

############ cat the idxstats
rm -f ALLSamples.idxstats; touch ALLSamples.idxstats
for f in tmp.mapping/*.bam; do
        N=$(echo $f | sed -e 's/\.bctrimmedreads\.fastq\.gz\.bam//g' -e 's/.*\///g')
        samtools idxstats $f | sed "s/^/$N\t/g" >> ALLSamples.idxstats
done

awk -F "\t" '{if($4 >= 1) print $0}' ALLSamples.idxstats | sed -z 's/^/ReadSource\tContig\tContigLength\tNumReadsMapped\tNumUnmappedReads\n/g' > tmp.ALLSamples.idxstats && mv tmp.ALLSamples.idxstats ALLSamples.idxstats

###########################################
fi
###########################################

CoverageDepth () {
for f in ${CoverageDepthInput}/*.bam; do
    N=$(basename $f)
    samtools index $f

    ## get contig lengths
    samtools idxstats $f | cut -f 1,2 | sed '/^\*/d' > ${f}.contiglengths

    ## total coverage per base for each contig
    samtools depth $f | awk -F "\t" -v N=$N '{print N"\t"$0}' > ${f}.samtoolsdepth

done

## contig length fix
rm -f tmp.ALLcontig; touch tmp.ALLcontig
cat tmp.mapping/*.contiglengths >> tmp.ALLcontig
sort -Vu tmp.ALLcontig -o tmp.ALLcontig

find CoverageDepthInput -type f -empty -print -delete

}


##########################################
if [[ "${step}" =~ "B" ]] || [[ "$step" =~ "b4" ]]; then
###########################################

########### breadth of coverage

CoverageDepthInput="tmp.mapping"
CoverageDepth

cd $CoverageDepthInput
sbatch --time=24:00:00 --cpus-per-task=12 --mem=240GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; Rscript /panfs/jay/groups/27/dcschroe/dmckeow/bioinf/bin/spillover-CovDepth.r'
cd ..
##### OBSOLETE - now do in R using normalised reads counts
## get total number of bases covered at MIN_COVERAGE_DEPTH or higher PER contig
##rm -f tmp.ALLstats.list.count0 tmp.ALLstats.list.count1 tmp.ALLstats.list.count2 tmp.ALLstats.list.count3; touch tmp.ALLstats.list.count0 tmp.ALLstats.list.count1 tmp.ALLstats.list.count2 tmp.ALLstats.list.count3
##rm -f tmp.ALLstats.list.count0.5; touch tmp.ALLstats.list.count0.5

##for f in tmp.mapping/*.samtoolsdepth; do
  ##  cut -f 1,2 $f | sort -Vu >> tmp.ALLstats.list.count0
   ## awk '$4 >= 5' $f | cut -f 1,2 >> tmp.ALLstats.list.count1
    ##awk '$4 >= 10' $f | cut -f 1,2 >> tmp.ALLstats.list.count2
    ##awk '$4 >= 100' $f | cut -f 1,2 >> tmp.ALLstats.list.count3
##done

##sort -V tmp.ALLstats.list.count0 | uniq -c | sed -E 's/ +(.+) (.+)/\2\t\1/g' > tmp.ALLstats.list.count0.5
##sort -V tmp.ALLstats.list.count1 | uniq -c | sed -E 's/ +(.+) (.+)/\2\t\1/g' > tmp.ALLstats && mv tmp.ALLstats tmp.ALLstats.list.count1
##sort -V tmp.ALLstats.list.count2 | uniq -c | sed -E 's/ +(.+) (.+)/\2\t\1/g' > tmp.ALLstats && mv tmp.ALLstats tmp.ALLstats.list.count2
##sort -V tmp.ALLstats.list.count3 | uniq -c | sed -E 's/ +(.+) (.+)/\2\t\1/g' > tmp.ALLstats && mv tmp.ALLstats tmp.ALLstats.list.count3


## add headers
sed -i 's/\.bam\t/\t/g' tmp.ALLcontig
sed -i 's/\.bam\t/\t/g' tmp.ALLstats.list.count*

sed -i -z 's/^/contig\tcontig_length\n/1' tmp.ALLcontig
##sed -i -z 's/^/reads\tcontig\n/1' tmp.ALLstats.list.count0
##sed -i -z 's/^/reads\tcontig\ttotal_bases_over_5_coverage\n/1' tmp.ALLstats.list.count1
##sed -i -z 's/^/reads\tcontig\ttotal_bases_over_10_coverage\n/1' tmp.ALLstats.list.count2
##sed -i -z 's/^/reads\tcontig\ttotal_bases_over_100_coverage\n/1' tmp.ALLstats.list.count3


#####################################
fi
#####################################



##########################################
if ([[ "${step}" =~ "B" ]] || [[ "$step" =~ "b5" ]]) && [[ ! "$reference" == "" ]]; then
###########################################

## map reads of each sample against reference genomes , to later identify the coverage of regions of interest - i.e. the RDRP primer target

#### ALL contigs vs ALL reads for each sample in the bin and count them per contig
rm -fr tmp.mapping.ref; mkdir tmp.mapping.ref
rm -f run_mapping_vs_ref; touch run_mapping_vs_ref

    #### prepare file to run mapping from
for F in $(echo ${reference} | sed 's/,/\n/g'); do
    for f in $(cat ${input}); do
        READS=$(echo $f | cut -d ";" -f 2)
        Q="sbatch --time=12:00:00 --cpus-per-task=4 --mem=48GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval \"\$(conda shell.bash hook)\"; conda activate bioinftools;"
        R="'"
        echo "$Q seqkit seq --min-len 200 ${READS} | minimap2 -t \$SLURM_CPUS_PER_TASK -ax map-ont --secondary=no $F - | samtools view -h -F 2308 | samtools sort -O BAM - > tmp.mapping.ref/$(basename $READS).$(basename $F).bam$R"

    done >> run_mapping_vs_ref
done
bash run_mapping_vs_ref

ISJOBRUNNING=$(squeue --me | grep -s "wrap" | wc -l); 
until [[ $ISJOBRUNNING -le 0 ]]; do 
    sleep 1m; ISJOBRUNNING=$(squeue --me | grep -s "wrap" | wc -l);
done;

#####################################
fi
#####################################

##########################################
if [[ "${step}" =~ "C" ]] || [[ "$step" =~ "c1" ]] && [[ ! "$reference" == "" ]]; then
###########################################

################ get the intervals by window size for sequences that are targets for mapping
for F in $(echo ${reference} | sed 's/,/\n/g'); do

    samtools faidx $F
    cut -f 1,2 ${F}.fai > tmp.$(basename $F).reffasta
    WINDOW=1000
while IFS=$'\t' read -r sequence length; do
        for ((start=1; start<=length; start+=WINDOW)); do
            end=$((start+WINDOW-1))
         if [ "$end" -gt "$length" ]; then
             end="$length"
         fi
            echo -e "$sequence\t$start\t$end"
        done
    done < tmp.$(basename $F).reffasta > $(basename ${F}).bed

done

rm -f tmp.mapping.ref/*.meanwindowdepth
rm -f tmp.mapping.ref/*.idxstats
rm -f tmp.mapping.ref/*.totalmapped

for F in $(echo ${reference} | sed 's/,/\n/g'); do

for f in tmp.mapping.ref/*.bam; do
    N=$(basename $f)
    bedtools coverage -mean -a $(basename ${F}).bed -b $f | awk -F "\t" -v N=$N '{print N"\t"$0}' > ${f}.$(basename ${F}).meanwindowdepth

    ##### get the contig length and total reads mapped for each reference sequence 

    samtools index $f
    samtools idxstats $f | sed '/^\*/d' | cut -f 1,2 | awk -F "\t" -v N=$N '{print N"\t"$0}' > ${f}.$(basename ${F}).idxstats
    NR=$(samtools view -c $f)
    sed -i "s/$/\t$NR/g" ${f}.$(basename ${F}).idxstats ## add the true num reads in library
    samtools view -c -F 2308 $f | awk -F "\t" -v N=$N '{print N"\t"$0}' > ${f}.$(basename ${F}).totalmapped ## total reads mapped per reference per sample
done

done

rm -f ALL.mapping.ref.idxstats; touch ALL.mapping.ref.idxstats

for f in tmp.mapping.ref/*.idxstats; do
    cat $f >> ALL.mapping.ref.idxstats
done

rm -f ALL.mapping.ref.meanwindowdepth; touch ALL.mapping.ref.meanwindowdepth

for f in tmp.mapping.ref/*.meanwindowdepth; do
    cat $f >> ALL.mapping.ref.meanwindowdepth
done

rm -f ALL.mapping.ref.totalmapped; touch ALL.mapping.ref.totalmapped

for f in tmp.mapping.ref/*.totalmapped; do
    cat $f >> ALL.mapping.ref.totalmapped
done

sort -Vu ALL.mapping.ref.totalmapped -o ALL.mapping.ref.totalmapped

### add headers
sed -i -z 's/^/specific_read_source\treference\tref_length\tNumReadsInSpecificReadSource\n/g' ALL.mapping.ref.idxstats
sed -i -z 's/^/specific_read_source\treference\tref_start\tref_end\tavg_coverage_by_window\n/g' ALL.mapping.ref.meanwindowdepth
sed -i -z 's/^/specific_read_source\ttotal_mapped_reads\n/g' ALL.mapping.ref.totalmapped

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

########################################################
########################################################
########################################################
########################################################

########################################## get bin info for anvio run (repeat for each binning run you want to include frmo within main folder - concatenate them into a single "bin_list" and run rest of script using that)
if [[ "$step" =~ "d1" ]]; then
###########################################
rm -f bin_list; touch bin_list

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

awk -F "\t" '{print $3"______"$6}' bin_list | sed 's/.*\///g' | paste bin_list - > tmp && mv tmp bin_list


awk -F "\t" 'BEGIN{OFS=FS} { gsub(/ .*/," ", $7) }1' bin_list > tmp.bin_list && mv tmp.bin_list bin_list
awk -F "\t" 'BEGIN{OFS=FS} { gsub(/ .*/," ", $6) }1' bin_list > tmp.bin_list && mv tmp.bin_list bin_list
sed -z -i 's/^/AnvioBin\tContigNameBin\tContigPath\tReadPath\tSampleNameBin\tContigNameCanu\tContigNameDerep\n/g' bin_list

###########################################
fi
###########################################



########################################################
########################################################
########################################################
########################################################


########################################## run checkv
if [[ "${step}" =~ "E" ]] || [[ "$step" =~ "e1" ]]; then
###########################################

checkv end_to_end all-contigs.fa checkv_out -t $THREADS -d /panfs/jay/groups/27/dcschroe/shared/checkv-db-v1.5

###########################################
fi
###########################################


########################################## download local
if [[ "${step}" == "download" ]]; then
###########################################
RProjectFolder="RProject_Spillover_FINAL"
rm -fr $RProjectFolder; mkdir $RProjectFolder
cp ALLSamples.idxstats $RProjectFolder/
cp SelfSamples.idxstats $RProjectFolder/
cp read_total_per_sample $RProjectFolder/
##cp tmp.ALLcontig $RProjectFolder/ 
##cp tmp.ALLstats.list.count1 $RProjectFolder/ 
##cp tmp.ALLstats.list.count2 $RProjectFolder/ 
##cp tmp.ALLstats.list.count3 $RProjectFolder/ 
cp bin_list $RProjectFolder/ 
cp bin_mmseq $RProjectFolder/ 
##cp ALL.mapping.ref.meanwindowdepth $RProjectFolder/ 
##cp ALL.mapping.ref.idxstats $RProjectFolder/ 
##cp ALL.mapping.ref.totalmapped $RProjectFolder/ 
cp checkv_out/*.tsv $RProjectFolder/ 
cp SeqsRemovedByClustering $RProjectFolder/
cp *_AllContigs.dmnd.blastx $RProjectFolder/
cp all-contigs.fa.kaiju.*.names $RProjectFolder/
cp AllContigLengths.tsv $RProjectFolder/

cp tmp.mapping/AllSamtoolsDepthCovMappingByWindow.tsv $RProjectFolder/AllSamtoolsDepthCovMappingByWindow.tsv
cp tmp.mapping.self/AllSamtoolsDepthCovMappingByWindow.tsv $RProjectFolder/SelfSamtoolsDepthCovMappingByWindow.tsv
#scp -r dmckeow@agate.msi.umn.edu:/panfs/jay/groups/27/dcschroe/dmckeow/data/Spillover_FINAL/$RProjectFolder/ /mnt/c/Users/Dean\ Mckeown/Downloads/Spillover_FINAL/

###########################################
fi
###########################################



########################################## bin split
if [[ "${step}" == "binsplit" ]]; then
###########################################
rm -fr ContigsBinsFinal; mkdir ContigsBinsFinal
while IFS=$'\t' read -r SuperBin RepresentativeName; do
    N=$(echo "${SuperBin}_${RepresentativeName}" | sed -E -e 's/[^a-zA-Z0-9_]/_/g' -e 's/_+/_/g' -e 's/^_|_$//g')
    touch ContigsBinsFinal/$N.fa
done < <(sed '1,1d' dfContigsManual.tsv | cut -f 8,9 | sort -Vu | sed -E '/Cellular;|ExcludedByLengthCutoff/d' )

while IFS=$'\t' read -r contig SuperBin RepresentativeName; do
    N=$(echo "${SuperBin}_${RepresentativeName}" | sed -E -e 's/[^a-zA-Z0-9_]/_/g' -e 's/_+/_/g' -e 's/^_|_$//g')
    seqkit grep -p $contig derep-contigs.fa >> ContigsBinsFinal/$N.fa
done < <(sed '1,1d' dfContigsManual.tsv | cut -f 1,8,9 | sed -E '/Cellular;|ExcludedByLengthCutoff/d' | sort -V -t $'\t' -k 2,3)



###########################################
fi
###########################################






###########################################
if [[ "${step}" == "SelfmapRun" ]]; then
###########################################

for f in $(cat ${input}); do
    READS=$(echo $f | cut -d ";" -f 2)
    cp $READS tmp.derep/
done

#### dereplicate further by specific merged barcodes, using last contig in each row as merged file name
sed -E -e 's/^(.*)\t(.*)$/\1\t\2.bctrimmedreads.fastq.gz > tmp \&\& mv tmp tmp.derep.merge\/\2/g' -e 's/^/zcat tmp.derep\//g' -e 's/$/.bctrimmedreads.fastq.gz/g' -e 's/\t/.bctrimmedreads.fastq.gz tmp.derep\//g' input_merge_bc > tmp.run_bc_merge
bash tmp.run_bc_merge

sed -E 's/^zcat (.*) >.*/rm -f \1/g' tmp.run_bc_merge > tmp.run_bc_merge.rm
bash tmp.run_bc_merge.rm
for f in tmp.derep.merge/*.bctrimmedreads.fastq.gz; do mv $f tmp.derep/ ; done

#### ALL contigs vs ALL reads for each sample in the bin and count them per contig
rm -fr tmp.mapping.self; mkdir tmp.mapping.self

    #### prepare file to run mapping from
    for f in $(cat ${input}); do
        READS=$(echo $f | cut -d ";" -f 2 | sed 's/\.bctrimmedreads\.fastq\.gz//g')
        Q="sbatch --time=24:00:00 --cpus-per-task=4 --mem=48GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval \"\$(conda shell.bash hook)\"; conda activate bioinftools;"
        R="'"
        echo "$Q seqkit seq --min-len 200 tmp.derep/$(basename $READS).bctrimmedreads.fastq.gz | minimap2 -t \$SLURM_CPUS_PER_TASK -ax map-ont --secondary=no tmp.derep/$(basename $READS).contigs.fasta - | samtools view -h -F 2308 | samtools sort -O BAM - > tmp.mapping.self/$(basename $READS).bam$R"

    done > run_mapping

bash run_mapping

###########################################
fi
###########################################

###########################################
if [[ "${step}" == "SelfmapCov" ]]; then
###########################################

#### get coverage depths across bases
CoverageDepthInput="tmp.mapping.self"
CoverageDepth

############ cat the idxstats
rm -f SelfSamples.idxstats; touch SelfSamples.idxstats
for f in tmp.mapping.self/*.bam; do
        N=$(echo $f | sed -e 's/\.bctrimmedreads\.fastq\.gz\.bam//g' -e 's/.*\///g')
        samtools idxstats $f | sed "s/^/$N\t/g" >> SelfSamples.idxstats
done

awk -F "\t" '{if($4 >= 1) print $0}' SelfSamples.idxstats | sed -z 's/^/ReadSource\tContig\tContigLength\tNumReadsMapped\tNumUnmappedReads\n/g' > tmp.SelfSamples.idxstats && mv tmp.SelfSamples.idxstats SelfSamples.idxstats

cd $CoverageDepthInput
sbatch --time=24:00:00 --cpus-per-task=12 --mem=240GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; Rscript /panfs/jay/groups/27/dcschroe/dmckeow/bioinf/bin/spillover-CovDepth.r'
cd ..
###########################################
fi
###########################################