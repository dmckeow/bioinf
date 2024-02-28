#!/bin/bash

umask 002
green='\033[32m';red='\033[31m';cyan='\033[36m';nocolor='\033[m'
echo -e "${green}RUNNING bioinf-binning-summarise with the following parameters: $@ ${nocolor}"

Help()
{
echo -e "########### This script is for generating various outputs for checking contigs within a bin. You must have ran anvio and anvi-summarise for this script to work. This script is designed to run on a single bin per job submission and within the folder for each bin. It some of the same inputs needed to run the anvio script (-i, -p, -o) ###########\n"
echo -e "########### REQUIRED arguments ###########\n-i --input\tfull path to the file provided via --input to the bioinf-binning.sh script"
echo -e "-s --step\tWhich script steps to run"
echo -e "\n-R [optional] --reference\tpath(s) to fasta files of reference sequences to include in alignment. Provide multiple paths separated by commas\n"
echo -e "\n\nTo proceed to steps B, you should manually inspect postbin_contigcheck_FINAL5.tsv and decide on a cutoff read coverage depth and contig length to filter out contigs may be artefacts or of poor assembly quality. To use only one of these parameters, simply use 0 for either parameter. You must then provide these values through the following arguments:"
echo -e "\n-d --mindepth\t the minimum average coverage depth a contig in your bin must have to pass contig curation\n"
echo -e "\n-l --minlength\t the minimum length for contigs and reference sequence cutoff filters\n"

}

### set euclidean options to default setting here

##

while getopts i:s:R:d:l:L:t:h option
do 
    case "${option}" in
    i)input=${OPTARG};;
    s)step=${OPTARG};;
    R)reference=${OPTARG};;
    d)mindepth=${OPTARG};;
    l)minlength=${OPTARG};;
    L)maxlength=${OPTARG};;
    t)threads=${OPTARG};;
    h)Help; exit;;
    esac

done

if [[ ! "$step" == "fetch" ]]; then
    if [[ -z "${input}" ]]; then echo "-i, --input REQUIRED"; Help; exit; fi
    if [[ -z "${mindepth}" ]]; then echo "-d, --mindepth REQUIRED"; Help; exit; fi
fi

if [[ "$step" == "fetch" ]]; then
    if [[ -z "${input}" ]]; then echo "-i, --input REQUIRED"; Help; exit; fi
    if [[ -z "${mindepth}" ]]; then echo "-d, --mindepth REQUIRED"; Help; exit; fi
fi

if [[ -z "${minlength}" ]]; then echo "-l, --minlength REQUIRED"; Help; exit; fi
if [[ -z "${maxlength}" ]]; then echo "-L, --maxlength REQUIRED"; Help; exit; fi

### THREADS

if [[ ! -z "${SLURM_CPUS_PER_TASK}" ]]; then THREADS="${SLURM_CPUS_PER_TASK}"; echo -e "${green}Using threads set by SLURM_CPUS_PER_TASK:${nocolor}"; fi

if [[ ! -z "${threads}" ]]; then THREADS="${threads}"; echo -e "${green}Using threads set by --threads flag:${nocolor}"; fi

if [[ -z "${threads}" ]] && [[ -z "${SLURM_CPUS_PER_TASK}" ]]; then THREADS="1"; echo -e "${green}No SLURM_CPUS_PER_TASK or --threads set, using default threads:${nocolor}"; fi
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

bin=$(pwd | sed 's/.*\///g')
tmp=$(echo -e "tmp.bioinf_binning_summarise_$bin") 
prefix=$(realpath ../../../ | sed 's/.*\/ANVIO_//g')

####################### THE SCRIPT #####################################


#################### STEP fetch - fetch references genomes from NCBI #################################
### NEEDS fixed
if [[ "$step" == "fetch" ]]; then

taxonomy="Iflavirus"
minsize="8000"
maxsize="11000"
cutoff="0.90"

#### first use a taxonomic name to fetch fastas and then filter it by min and max sizes range
esearch -db nuccore -query "${taxonomy}" | efetch -format fasta | seqkit seq -j $THREADS -m $minsize -M $maxsize - | seqkit rmdup -j $THREADS -s - | sed -E -e 's/[^>a-zA-Z0-9.]/_/g' -e 's/_+/_/g' -e 's/>/>REF_/g' > tmp.${taxonomy}.fetched_reference_fastas

### remove redundant sequences by a similarity threshold
cd-hit-est -c ${cutoff} -i tmp.${taxonomy}.fetched_reference_fastas -o tmp.${taxonomy}.fetched_reference_fastas.cd-hit

grep ">" tmp.${taxonomy}.fetched_reference_fastas.cd-hit | sed 's/>//g' > curated_sequence_list

#### remove unwanted sequences from curated_sequence_list then use seqkit to pull out a subset you want for analyses from this file
seqkit grep -n -f curated_sequence_list tmp.${taxonomy}.fetched_reference_fastas.cd-hit > tmp.${taxonomy}.fetched_reference_fastas.final

#### now fix the sequence names as you want to look good on your tree, etc

grep ">" tmp.${taxonomy}.fetched_reference_fastas.final > og_ref_names
cp og_ref_names new_ref_names
paste -d "/" og_ref_names new_ref_names | awk -F "/" '{print "s/"$1"/"$2"/g"}' | sed -f - tmp.${taxonomy}.fetched_reference_fastas.final > ${taxonomy}.fetched_curated_reference.fa

exit

##############################################################
fi
##############################################################




##### fix the bin contig file if anvio fucked it up by literally splitting the contigs up (it does this for binning without mapping for some reason!)
if grep -q partial_[0-9]*_[0-9]* "${bin}-contigs.fa"; then
    sed -E 's/(.+)_split.+/\1/' ${bin}-original_split_names.txt | sort -Vu > ${tmp}/tmp.nomap_binfix
    seqkit grep -n -f ${tmp}/tmp.nomap_binfix ../../../${prefix}.fa > ${bin}-contigs.fa
fi

########################################## STEP A1 - prepare mapping plots for the contigs within the bin
if [[ -z "${step}" ]] || [[ "$step" =~ "A1" ]]; then
###########################################
  

rm -fr $tmp; mkdir $tmp
rm -fr figures; mkdir figures


#### input prep
grep ">" ${bin}-contigs.fa | sed 's/>//g' | sort -Vu > ${tmp}/tmp.bin_contig_fa_deflines
grep ">" ${bin}-contigs.fa | sed 's/>//g' | sed -E 's/contig_(.*)_[0-9]+$/\1/g' | sort -Vu > ${tmp}/tmp.bin_contig_fa_deflines_samplesnames
grep -w -f ${tmp}/tmp.bin_contig_fa_deflines_samplesnames $input > ${tmp}/tmp.bin_input_filtered


########### sample chromosome loci depth


for SETS in $(cat ${tmp}/tmp.bin_input_filtered); do
  PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
    seqkit grep -rp contig_${PROFILENAME}_[0-9]+ ../../../${prefix}.fa > ${tmp}/tmp.bin_filter_${PROFILENAME}.fa
done

#### ALL contigs vs ALL reads for each sample in the bin

for SETS in $(cat ${tmp}/tmp.bin_input_filtered); do
    PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
    READS=$(echo $SETS | cut -d ";" -f 2)
    minimap2 -t $THREADS -ax map-ont ${tmp}/tmp.bin_filter_${PROFILENAME}.fa $READS | samtools sort -O BAM - > ${tmp}/tmp.bin_filter_${PROFILENAME}.bam
    samtools index -@ $THREADS -b ${tmp}/tmp.bin_filter_${PROFILENAME}.bam
    samtools depth ${tmp}/tmp.bin_filter_${PROFILENAME}.bam | awk -v P="${PROFILENAME}" '{print P"\t"$0}' > ${tmp}/tmp.bin_filter_${PROFILENAME}.depth
done


### concatenate and filter down to only the contigs present in the bin
cat ${tmp}/tmp.bin_filter_*.depth | grep -f ${tmp}/tmp.bin_contig_fa_deflines - > ${bin}.bioinf_binning_summarise.FINAL1.tsv

########## sample total_number_reads /path/to_bam

for SETS in $(cat ${tmp}/tmp.bin_input_filtered); do
    PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
    READS=$(echo $SETS | cut -d ";" -f 2)
    NO_READS=$(seqkit stat -T $READS | sed '1,1d' | cut -f 4)
    realpath ${tmp}/tmp.bin_filter_${PROFILENAME}.bam | awk -v P="${PROFILENAME}" -v R="$NO_READS" '{print P"\t"R"\t"$0}' > ${tmp}/tmp.bin_read_numbers_${PROFILENAME}.txt
done

cat ${tmp}/tmp.bin_read_numbers_*.txt > ${bin}.bioinf_binning_summarise.FINAL2.tsv


sed 's/ /\t/g' ../../../*.deflinekey | cut -f 1,3 | sed 's/len=//g' | grep -w -f ${tmp}/tmp.bin_contig_fa_deflines - > ${bin}.bioinf_binning_summarise.FINAL3.tsv

Rscript /panfs/jay/groups/27/dcschroe/heske011/scripts/archive/plotRPKM_and_NumReads_4APRIL2023.r


#####################################
fi
#####################################





########################################## STEP A2 - prepare table summarising properties
if [[ -z "${step}" ]] || [[ "$step" =~ "A2" ]]; then
###########################################

#### if user provides no references fastas, then just set the input for the blast to the bin_contigs
if [ -z "$reference" ]; then
    cp ${bin}-contigs.fa ${tmp}/tmp.fastas_to_blast

 #### otherwise, combine the reference fastas and bin contigs
elif [ ! -z "$reference" ]; then
    cp ${bin}-contigs.fa ${tmp}/tmp.fastas_to_blast
    echo -e "${reference}" | sed 's/,/\n/g' > ${tmp}/tmp.reference
    for f in $(cat ${tmp}/tmp.reference); do
        cat $f >> ${tmp}/tmp.fastas_to_blast
    done
fi

cd ${tmp}
awk -F "\t" '{print $4 > "tmp."$2".a1"}' ../postbin_contigcheck_FINAL6.tsv
cd ..

for f in ${tmp}/tmp.*.a1; do
    awk -F "\t" '{ sum += $1; n++ } END { if (n > 0) printf "%0.0f\n", sum / n}' $f > ${tmp}/$(basename $f .a1).a2
done

for f in ${tmp}/tmp.*.a2; do awk -v f=$f '{print f"\t"$1}' $f | sed -E 's/tmp\.(.+)\.a2/\1/g'; done | sort -t $'\t' -k 1,1V | cut -f 2 > ${tmp}/average_coverage_depth_percontig

sort -t $'\t' -k 1,1V -o postbin_contigcheck_FINAL3.tsv postbin_contigcheck_FINAL3.tsv

cut -f 1,2 postbin_contigcheck_FINAL1.tsv | sort -Vu | sort -t $'\t' -k 2,2V | paste - postbin_contigcheck_FINAL3.tsv average_coverage_depth_percontig | cut --complement -f 2 | sort -t $'\t' -k 4,4nr > postbin_contigcheck_FINAL5.tsv
sed -i "s/$/\t${bin}/g" postbin_contigcheck_FINAL5.tsv

sed -i "s/$/\t${minlength}/g" postbin_contigcheck_FINAL5.tsv

sed -z -i 's/^/sample\tcontig\tlength_bp\taverage_depth\tbin\tREF_genome_length\n/1' postbin_contigcheck_FINAL5.tsv
sed -i 's/\t/;/g' postbin_contigcheck_FINAL5.tsv

#####################################
fi
#####################################

########################################## STEP B1 - alignment of final curated contigs
if [[ -z "${step}" ]] || [[ "$step" =~ "B1" ]]; then
###########################################

awk -F ";" -v minlength=${minlength} -v mindepth=${mindepth} '{if($3 > minlength && $4 > mindepth && $0 !~ "sample;contig;length") print $0";yes;"minlength";"mindepth; else if($0 ~ "sample;contig;length") print $0";minlengthdepthpass;minlength;mindepth";  else print $0";no;"minlength";"mindepth}' postbin_contigcheck_FINAL5.tsv > tmp && mv tmp postbin_contigcheck_FINAL5.tsv

sed '1,1d' postbin_contigcheck_FINAL5.tsv | awk -F ";" -v minlength=${minlength} -v mindepth=${mindepth} '{if($3 > minlength && $4 > mindepth) print $2}' > ${tmp}/tmp.B001

grep ">REF_" ${tmp}/tmp.fastas_to_blast | sed 's/>//g' >> ${tmp}/tmp.B001

seqkit grep -n -f ${tmp}/tmp.B001 ${tmp}/tmp.fastas_to_blast > ${bin}-d${mindepth}-l${minlength}-finalcontigs.fa

mafft --adjustdirectionaccurately --thread $THREADS --auto ${bin}-d${mindepth}-l${minlength}-finalcontigs.fa > ${bin}-d${mindepth}-l${minlength}-finalcontigs.aln

FastTree -gtr -nt < ${bin}-d${mindepth}-l${minlength}-finalcontigs.aln > ${bin}-d${mindepth}-l${minlength}-finalcontigs.fasttree


### blastx whole splits against custom dmnd databases
### reduce to best hit per contig and VOG hit only (by evalue)
### set the max target sequences based on the longest contig, divided by 100 (at the shorter end of genes)
MAXTARGETS=$(seqkit stat -T ${bin}-d${mindepth}-l${minlength}-finalcontigs.fa | sed '1,1d' | awk -F "\t" '{printf "%0.0f\n", $8/100}')

diamond blastx --max-target-seqs $MAXTARGETS --evalue 1e-80 --sensitive -p $THREADS -d $bioinfdb/DMND/vog.dmnd -q ${bin}-d${mindepth}-l${minlength}-finalcontigs.fa -f 6 -o ${bin}-d${mindepth}-l${minlength}-finalcontigs.dmnd.blastx

awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/\..*/,"",$2)}1' ${bin}-d${mindepth}-l${minlength}-finalcontigs.dmnd.blastx | awk -F "\t" '!a[$1,$2]++' | awk -F "\t" 'BEGIN{OFS="\t"};{print $1,"blastx","CDS",$7,$8,".","-",".","blastx_subject="$2}' > ${bin}-d${mindepth}-l${minlength}-finalcontigs.gff

#####################################
fi
#####################################


########################################## STEP test -convert annotations to the actual gff format
if [[ -z "${step}" ]] || [[ "$step" =~ "B2" ]]; then
###########################################

awk -F "\t" 'BEGIN{OFS="\t"};{print $2,"prodigal","CDS",$3,$4,".",$5,".","ID="$1}' ${bin}-gene_calls.txt | sed -e 's/\tr\t/\t-\t/g' -e 's/\tf\t/\t+\t/g' -e '/^contig\tprodigal\tCDS\tstart/d' > ${tmp}/${bin}-gene_calls.tmp1

cut --complement -f 1-5 ${bin}-gene_calls.txt > ${tmp}/${bin}-gene_calls.tmp2

head -1 ${tmp}/${bin}-gene_calls.tmp2 | tr "\t" "\n" > ${tmp}/${bin}-gene_calls.tmp3
sed -i -E -e 's/[^a-zA-Z0-9_]/_/g' -e 's/_+/_/g' -e 's/^_|_$//g' ${tmp}/${bin}-gene_calls.tmp3

sed -i '1,1d' ${tmp}/${bin}-gene_calls.tmp2

rm -f ${tmp}/${bin}-*.tmp4

id=0
while read filename; do
  id=$((id+1))
 tail -n +1 ${tmp}/${bin}-gene_calls.tmp2 | cut -f $id >> ${tmp}/${bin}-"$filename".tmp4
done < ${tmp}/${bin}-gene_calls.tmp3

rm -f ${tmp}/${bin}-dna_sequence.tmp4

sed -i 's/^$/-/g' ${tmp}/${bin}-*.tmp4

for f in ${tmp}/${bin}-*.tmp4; do
    awk -v f=$(basename $f .tmp4) '{print f"="$0}' $f > ${tmp}/$(basename $f .tmp4).tmp5
done

paste -d ";" ${tmp}/${bin}-gene_calls.tmp1 ${tmp}/${bin}-*.tmp5 > ${bin}-gene_calls.gff

rm -f ${tmp}/${bin}-*.tmp4

#####################################
fi
#####################################



