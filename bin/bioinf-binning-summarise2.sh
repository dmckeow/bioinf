#!/bin/bash

umask 002
green='\033[32m';red='\033[31m';cyan='\033[36m';nocolor='\033[m'
echo -e "${green}RUNNING bioinf-binning-summarise with the following parameters: $@ ${nocolor}"

Help()
{
echo -e "########### This script is for generating various outputs for checking contigs within a bin. You must have ran anvio and anvi-summarise for this script to work. This script is designed to run on a single bin per job submission and within the folder for each bin. RUN this script within the BINNING_ output folder of your binning run ###########\n"
echo -e "-i -input [REQUIRED]\tFull path to the same -input file used by the bioinf-binning run"
echo -e "-s -step\tWhich script steps to run"
echo -e "-t -threads\tNumber of parallel processes"

}

### set euclidean options to default setting here

##

while getopts i:s:t:h option
do 
    case "${option}" in
    i)input=${OPTARG};;
    s)step=${OPTARG};;
    t)threads=${OPTARG};;
    h)Help; exit;;
    esac

done

if [[ -z "${input}" ]]; then echo "-i, -input REQUIRED"; Help; exit; fi


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
PROJECT_NAME=$(basename $(pwd) -SUMMARY | sed 's/^BINNING_//g')

####################### THE SCRIPT #####################################

#################################################################
for BIN in WPF_Bombus_final_no_IOC-SUMMARY/bin_by_bin/Bin_11_Eukaryote*; do
#################################################################

    cd ${BIN}
    BIN_NAME=$(basename $(pwd))
    rm -fr $tmp; mkdir $tmp
    rm -fr figures; mkdir figures

##### fix the bin contig file if anvio fucked it up by literally splitting the contigs up (it does this for binning without mapping for some reason!)
    if grep -q partial_[0-9]*_[0-9]* "${BIN_NAME}-contigs.fa"; then
        sed -E 's/(.+)_split.+/\1/' ${BIN_NAME}-original_split_names.txt | sort -Vu > ${tmp}/tmp.nomap_binfix
        seqkit grep -n -f ${tmp}/tmp.nomap_binfix ../../../${PROJECT_NAME}.fa > ${BIN_NAME}-contigs.fa
    fi

########################################## STEP A1 - prepare mapping plots for the contigs within the bin
    if [[ -z "${step}" ]] || [[ "$step" =~ "A1" ]]; then
###########################################
  
#### input prep
    grep ">" ${BIN_NAME}-contigs.fa | sed 's/>//g' | sort -Vu > ${tmp}/tmp.bin_contig_fa_deflines
    grep ">" ${BIN_NAME}-contigs.fa | sed 's/>//g' | sed -E 's/contig_(.*)_[0-9]+$/\1/g' | sort -Vu > ${tmp}/tmp.bin_contig_fa_deflines_samplesnames
    grep -w -f ${tmp}/tmp.bin_contig_fa_deflines_samplesnames $input > ${tmp}/tmp.bin_input_filtered


########### sample chromosome loci depth


    for SETS in $(cat ${tmp}/tmp.bin_input_filtered); do
        PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
        seqkit grep -rp contig_${PROFILENAME}_[0-9]+ ../../../${PROJECT_NAME}.fa > ${tmp}/tmp.bin_filter_${PROFILENAME}.fa
    done

#### ALL contigs vs ALL reads for each sample in the bin

    for SETS in $(cat ${tmp}/tmp.bin_input_filtered); do
        PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
        READS=$(echo $SETS | cut -d ";" -f 2)

        seqkit seq --min-len 200 ${READS} | minimap2 -t $THREADS -ax map-ont --secondary=no ${tmp}/tmp.bin_filter_${PROFILENAME}.fa - | samtools sort -O BAM - > ${tmp}/tmp.bin_filter_${PROFILENAME}.bam
# | samtools view -h -F 0x900 
        samtools index -@ $THREADS -b ${tmp}/tmp.bin_filter_${PROFILENAME}.bam
        samtools depth ${tmp}/tmp.bin_filter_${PROFILENAME}.bam | awk -v P="${PROFILENAME}" '{print P"\t"$0}' > ${tmp}/tmp.bin_filter_${PROFILENAME}.depth
    done


### concatenate and filter down to only the contigs present in the bin
    cat ${tmp}/tmp.bin_filter_*.depth | grep -f ${tmp}/tmp.bin_contig_fa_deflines - > ${BIN_NAME}.bioinf_binning_summarise.FINAL1.tsv

########## sample total_number_reads /path/to_bam

    for SETS in $(cat ${tmp}/tmp.bin_input_filtered); do
        PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
        READS=$(echo $SETS | cut -d ";" -f 2)
        NO_READS=$(seqkit stat -T $READS | sed '1,1d' | cut -f 4)
        realpath ${tmp}/tmp.bin_filter_${PROFILENAME}.bam | awk -v P="${PROFILENAME}" -v R="$NO_READS" '{print P"\t"R"\t"$0}' > ${tmp}/tmp.bin_read_numbers_${PROFILENAME}.txt
    done

    cat ${tmp}/tmp.bin_read_numbers_*.txt > ${BIN_NAME}.bioinf_binning_summarise.FINAL2.tsv


    sed 's/ /\t/g' ../../../${PROJECT_NAME}.deflinekey | cut -f 1,3 | sed 's/len=//g' | grep -w -f ${tmp}/tmp.bin_contig_fa_deflines - > ${BIN_NAME}.bioinf_binning_summarise.FINAL3.tsv

    Rscript /panfs/jay/groups/27/dcschroe/heske011/scripts/archive/plotRPKM_and_NumReads_4APRIL2023.r


#####################################
fi
#####################################





########################################## STEP A2 - prepare table summarising properties
if [[ -z "${step}" ]] || [[ "$step" =~ "A2" ]]; then
###########################################

cd ${tmp}
awk -F "\t" '{print $4 > "tmp."$2".a1"}' ../postbin_contigcheck_FINAL6.tsv
cd ..

for f in ${tmp}/tmp.*.a1; do
    awk -F "\t" '{ sum += $1; n++ } END { if (n > 0) printf "%0.0f\n", sum / n}' $f > ${tmp}/$(basename $f .a1).a2
done

for f in ${tmp}/tmp.*.a2; do awk -v f=$f '{print f"\t"$1}' $f | sed -E 's/tmp\.(.+)\.a2/\1/g'; done | sort -t $'\t' -k 1,1V | cut -f 2 > ${tmp}/average_coverage_depth_percontig

sort -t $'\t' -k 1,1V -o postbin_contigcheck_FINAL3.tsv postbin_contigcheck_FINAL3.tsv

cut -f 1,2 postbin_contigcheck_FINAL1.tsv | sort -Vu | sort -t $'\t' -k 2,2V | paste - postbin_contigcheck_FINAL3.tsv average_coverage_depth_percontig | cut --complement -f 2 | sort -t $'\t' -k 4,4nr > postbin_contigcheck_FINAL5.tsv
sed -i "s/$/\t${BIN_NAME}/g" postbin_contigcheck_FINAL5.tsv

sed -i "s/$/\t${minlength}/g" postbin_contigcheck_FINAL5.tsv

sed -z -i 's/^/sample\tcontig\tlength_bp\taverage_depth\tbin\tREF_genome_length\n/1' postbin_contigcheck_FINAL5.tsv
sed -i 's/\t/;/g' postbin_contigcheck_FINAL5.tsv

#####################################
fi
#####################################

########################################## STEP B1 - alignment of final curated contigs
if [[ -z "${step}" ]] || [[ "$step" =~ "B1" ]]; then
###########################################

### blastx whole splits against custom dmnd databases
### reduce to best hit per contig and VOG hit only (by evalue)
### set the max target sequences based on the longest contig, divided by 100 (at the shorter end of genes)
MAXTARGETS=$(seqkit stat -T ${BIN_NAME}-d${mindepth}-l${minlength}-finalcontigs.fa | sed '1,1d' | awk -F "\t" '{printf "%0.0f\n", $8/100}')

diamond blastx --max-target-seqs $MAXTARGETS --evalue 1e-80 --sensitive -p $THREADS -d $bioinfdb/DMND/vog.dmnd -q ${BIN_NAME}-d${mindepth}-l${minlength}-finalcontigs.fa -f 6 -o ${BIN_NAME}-d${mindepth}-l${minlength}-finalcontigs.dmnd.blastx

awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/\..*/,"",$2)}1' ${BIN_NAME}-d${mindepth}-l${minlength}-finalcontigs.dmnd.blastx | awk -F "\t" '!a[$1,$2]++' | awk -F "\t" 'BEGIN{OFS="\t"};{print $1,"blastx","CDS",$7,$8,".","-",".","blastx_subject="$2}' > ${BIN_NAME}-d${mindepth}-l${minlength}-finalcontigs.gff

#####################################
fi
#####################################


########################################## STEP test -convert annotations to the actual gff format
if [[ -z "${step}" ]] || [[ "$step" =~ "B2" ]]; then
###########################################

awk -F "\t" 'BEGIN{OFS="\t"};{print $2,"prodigal","CDS",$3,$4,".",$5,".","ID="$1}' ${BIN_NAME}-gene_calls.txt | sed -e 's/\tr\t/\t-\t/g' -e 's/\tf\t/\t+\t/g' -e '/^contig\tprodigal\tCDS\tstart/d' > ${tmp}/${BIN_NAME}-gene_calls.tmp1

cut --complement -f 1-5 ${BIN_NAME}-gene_calls.txt > ${tmp}/${BIN_NAME}-gene_calls.tmp2

head -1 ${tmp}/${BIN_NAME}-gene_calls.tmp2 | tr "\t" "\n" > ${tmp}/${BIN_NAME}-gene_calls.tmp3
sed -i -E -e 's/[^a-zA-Z0-9_]/_/g' -e 's/_+/_/g' -e 's/^_|_$//g' ${tmp}/${BIN_NAME}-gene_calls.tmp3

sed -i '1,1d' ${tmp}/${BIN_NAME}-gene_calls.tmp2

rm -f ${tmp}/${BIN_NAME}-*.tmp4

id=0
while read filename; do
  id=$((id+1))
 tail -n +1 ${tmp}/${BIN_NAME}-gene_calls.tmp2 | cut -f $id >> ${tmp}/${BIN_NAME}-"$filename".tmp4
done < ${tmp}/${BIN_NAME}-gene_calls.tmp3

rm -f ${tmp}/${BIN_NAME}-dna_sequence.tmp4

sed -i 's/^$/-/g' ${tmp}/${BIN_NAME}-*.tmp4

for f in ${tmp}/${BIN_NAME}-*.tmp4; do
    awk -v f=$(basename $f .tmp4) '{print f"="$0}' $f > ${tmp}/$(basename $f .tmp4).tmp5
done

paste -d ";" ${tmp}/${BIN_NAME}-gene_calls.tmp1 ${tmp}/${BIN_NAME}-*.tmp5 > ${BIN_NAME}-gene_calls.gff

rm -f ${tmp}/${BIN_NAME}-*.tmp4

#####################################
fi
#####################################



#################################################################
done
#################################################################