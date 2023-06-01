#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=120GB
#SBATCH --partition ag2tb
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

echo -e "COMMAND LINE JOB SUBMISSSION:\n\tsbatch /home/dcschroe/dmckeow/projects/DWV/script/postbin_contigcheck.sh $@"

Help()
{
echo -e "########### This script is for generating various outputs for checking contigs within a bin. You must have ran anvio and anvi-summarise for this script to work. This script is designed to run on a single bin per job submission. It some of the same inputs needed to run the anvio script (-i, -p, -o) ###########\n"
echo -e "########### REQUIRED arguments ###########\n-i --input\tfull path to a file that lists input file locations, semi-colon delimited, 3 fields:\n"
echo -e "\tFIELD1=contig fasta, separate line per assembly (full path). Can be individual assemblies or co-assemblies."
echo -e "\tFIELD2=reads fq.gz, corresponding to each assembly contig fasta (full path)"
echo -e "\tFIELD3=profile name; a short name for each profile which will appear in all outputs and figures."
echo -e "!!! profile names MUST NOT begin with a digit, and MUST ONLY include letters, digits, and underscores and MUST be UNIQUE within your samples being run!!!\n"
echo -e "! example file for --input:\n\n/home/data/19AUG22_1_DM_barcode01.contigs.fasta;/home/data/19AUG22_1_DM_barcode01-trimmedbcs.fastq.gz;location1\n/home/data/19AUG22_1_DM_barcode03.contigs.fasta;/home/data/19AUG22_1_DM_barcode03-trimmedbcs.fastq.gz;location2\n/home/data/19AUG22_1_DM_barcode06.contigs.fasta;/home/data/19AUG22_1_DM_barcode06-trimmedbcs.fastq.gz;location3\n\n"
echo -e "-s --step\tWhich script steps to run"
echo -e "\n-o --output\tabsolute path for final output directory that ALREADY EXISTS. An output folder named after ANVIO_(prefix) will be created here\n"
echo -e "\n-p --prefix\ta meaningful name for output folder that must NOT ALREADY EXIST. Will also be the name of the merged analyses files. It should be a UNIQUE name that cannot possibly be identical to any other folder within your output destination e.g. -p 11OCT22project will make the folder ANVIO_11OCT22project\n"
echo -e "\n-B --bin\tthe exact name of your bin - must be identical to bin folder name\n"
echo -e "\n-m --minsize\tthe min size (in bp) of your full genome size expected for this bin. Will be used to generate statistics\n"
echo -e "\n-M --maxsize\tthe max size (in bp) of your full genome size expected for this bin. Will be used to generate statistics\n"
echo -e "\n-T --taxonomy\ta taxonomic keyword(s) that will match exactly something in the full taxonomy info for reference genomes for this bin. Helps ensure appropriate taxa are used. If using multiple keywords, separate them with commas, and separate phrases or multi-part names with underscores\n"
echo -e "\n-R [optional] --reference\ta genome name keyword(s) that will match exactly reference genomes names for this bin. You can use taxids or accession numbers too! Can ensure that specific reference genomes are retrieved. If using multiple keywords, separate them with commas, and separate phrases or multi-part names with underscores\n"
echo -e "\n\nSteps in A are automatic. To proceed to steps B, you should manually inspect postbin_contigcheck_FINAL5.tsv and decide on a cutoff read coverage depth and contig length to filter out contigs may be artefacts or of poor assembly quality. To use only one of these parameters, simply use 0 for either parameter. You must then provide these values through the following arguments:"
echo -e "\n-d --mindepth\t the minimum average coverage depth a contig in your bin must have to pass contig curation\n"
echo -e "\n-l --minlength\t the minimum length of a contig in your bin must have to pass contig curation\n"




}


while getopts i:s:o:p:B:m:M:T:R:d:l:h option
do 
    case "${option}" in
    i)input=${OPTARG};;
    s)step=${OPTARG};;
    o)output=${OPTARG};;
    p)prefix=${OPTARG};;
    B)bin=${OPTARG};;
    m)minsize=${OPTARG};;
    M)maxsize=${OPTARG};;
    T)taxonomy=${OPTARG};;
    R)reference=${OPTARG};;
    d)mindepth=${OPTARG};;
    l)minlength=${OPTARG};;
    h)Help; exit;;
    esac

done

if [[ -z "${input}" ]]; then echo "-i, --input REQUIRED"; Help; exit; fi
if [[ -z "${step}" ]]; then echo "-s, --step REQUIRED"; Help; exit; fi
if [[ -z "${output}" ]]; then echo "-o, --output REQUIRED"; Help; exit; fi
if [[ -z "${prefix}" ]]; then echo "-p, --prefix REQUIRED"; Help; exit; fi
if [[ -z "${minsize}" ]]; then echo "-m, --minsize REQUIRED"; Help; exit; fi
if [[ -z "${maxsize}" ]]; then echo "-M, --maxsize REQUIRED"; Help; exit; fi
if [[ -z "${bin}" ]]; then echo "-B, --bin REQUIRED"; Help; exit; fi
if [[ -z "${taxonomy}" ]]; then echo "-T, --taxonomy REQUIRED"; Help; exit; fi

if [[ "$step" =~ "B1" ]] && [[ -z "${mindepth}" ]]; then echo "-d, --mindepth REQUIRED"; Help; exit; fi
if [[ "$step" =~ "B1" ]] && [[ -z "${minlength}" ]]; then echo "-l, --minlength REQUIRED"; Help; exit; fi


####################### PREREQUISITES #####################################

####################### SCRIPT PURPOSE #####################################

####################### HOW TO RUN SCRIPT #####################################

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################

O_N="postbin_contigcheck" #### output suffix for all outputs from this script
tmp="tmp.${O_N}"

####################### OUTPUT FILES #####################################

####################### SOFTWARE, DATABSES, ETC #####################################
### LOAD software available via shared environment on server:
module purge
eval "$(conda shell.bash hook)"

CUSTOM_DMND="/panfs/jay/groups/27/dcschroe/dmckeow/custom_DMND"

########################################################################
########################################################################

######### SOFTWARE/OTHER SCRIPTS THAT NEED SETUP BEFORE RUNNING:
seqkit="/home/dcschroe/dmckeow/seqkit"

##### Diamond - installed locally
DIAMOND="/panfs/jay/groups/27/dcschroe/shared/tools/diamond"

##### NCBI DIRECT ENTREZ TOOLS
EDIRECT="/home/dcschroe/dmckeow/edirect"

######################################################################
#################VOGDB#####################################################


####################### THE SCRIPT #####################################
cd ${output}/ANVIO_${prefix}/${prefix}-SUMMARY/bin_by_bin/${bin}

##### fix the bin contig file if anvio fucked it up by literally splitting the contigs up (it does this for binning without mapping for some reason!)
if grep -q partial_[0-9]*_[0-9]* "${bin}-contigs.fa"; then
    sed -E 's/(.+)_split.+/\1/' ${bin}-original_split_names.txt | sort -Vu > ${tmp}.000
    $seqkit grep -n -f ${tmp}.000 ${output}/ANVIO_${prefix}/${prefix}.fa > ${bin}-contigs.fa
fi

########################################## STEP A1 - prepare mapping pltos for the contigs within the bin
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A1" ]]; then
###########################################
module purge; conda deactivate
  

rm -f tmp.*${O_N}*


#### input prep
grep ">" ${bin}-contigs.fa | sed 's/>//g' | sort -Vu > tmp.${O_N}.00
grep ">" ${bin}-contigs.fa | sed 's/>//g' | sed -E 's/contig_(.*)_[0-9]+$/\1/g' | sort -Vu > tmp.${O_N}.0
grep -w -f tmp.${O_N}.0 $input > tmp.${O_N}.1


########### sample chromosome loci depth

module purge; conda deactivate
  module load minimap2/2.17
  module load samtools


for SETS in $(cat tmp.${O_N}.0); do
  PROFILENAME=$(echo $SETS)
    $seqkit grep -rp contig_${PROFILENAME}_[0-9]+ ${output}/ANVIO_${prefix}/${prefix}.fa > tmp.${PROFILENAME}.${O_N}.fa
done

#### ALL contigs vs ALL reads for each sample in the bin

for SETS in $(cat tmp.${O_N}.1); do
    PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
    READS=$(echo $SETS | cut -d ";" -f 2)
    minimap2 -t $SLURM_CPUS_PER_TASK -ax map-ont tmp.${PROFILENAME}.${O_N}.fa $READS | samtools sort -O BAM - > tmp.${PROFILENAME}.${O_N}.bam
    samtools index -@ $SLURM_CPUS_PER_TASK -b tmp.${PROFILENAME}.${O_N}.bam
    samtools depth tmp.${PROFILENAME}.${O_N}.bam | awk -v P="${PROFILENAME}" '{print P"\t"$0}' > tmp.${PROFILENAME}.${O_N}.2
done


### cocnatenate and filter down to only the contigs present in the bin
cat tmp.*.${O_N}.2 | grep -f tmp.${O_N}.00 - > ${O_N}_FINAL1.tsv

########## sample total_number_reads /path/to_bam

for SETS in $(cat tmp.${O_N}.1); do
    PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
    READS=$(echo $SETS | cut -d ";" -f 2)
    NO_READS=$($seqkit stat -T $READS | sed '1,1d' | cut -f 4)
    realpath ${PROFILENAME}.${O_N}.bam | awk -v P="${PROFILENAME}" -v R="$NO_READS" '{print P"\t"R"\t"$0}' > tmp.${PROFILENAME}.${O_N}.3
done

cat tmp.*.${O_N}.3 > ${O_N}_FINAL2.tsv

######### contig length
######cat tmp.*.${O_N}.fa | grep ">" | sed 's/>//g' | sort -Vu > tmp.${O_N}.4

sed 's/ /\t/g' ${output}/ANVIO_${prefix}/${prefix}.deflinekey | cut -f 1,3 | sed 's/len=//g' | grep -w -f tmp.${O_N}.00 - > ${O_N}_FINAL3.tsv

rm -fr figures; mkdir figures
##conda activate /panfs/jay/groups/27/dcschroe/heske011/.conda/envs/R-env
##Rscript /panfs/jay/groups/27/dcschroe/heske011/scripts/archive/plotRPKM_and_NumReads_4APRIL2023.r


#####################################
fi
#####################################

########################################## STEP A2 - prepare downsampling filters for plots to test mapping
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A2" ]]; then
###########################################
module purge; conda deactivate
module load samtools
module load minimap2/2.17

##### this step includes downsampling by contig size
#### ALL contigs vs ALL reads for each sample in the bin


#### input prep
grep ">" ${bin}-contigs.fa | sed 's/>//g' | sort -Vu > tmp.${O_N}.00
grep ">" ${bin}-contigs.fa | sed 's/>//g' | sed -E 's/contig_(.*)_[0-9]+$/\1/g' | sort -Vu > tmp.${O_N}.0
grep -w -f tmp.${O_N}.0 $input > tmp.${O_N}.1


for i in {1..10}; do
  LENGTH=$(($minsize/10*$i))
for SETS in $(cat tmp.${O_N}.1); do
    PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
    READS=$(echo $SETS | cut -d ";" -f 2)
    $seqkit seq -m $LENGTH tmp.${PROFILENAME}.${O_N}.fa | minimap2 -t $SLURM_CPUS_PER_TASK -ax map-ont - $READS | samtools sort -O BAM - > tmp.${LENGTH}.${PROFILENAME}.${O_N}.bam
    samtools index -@ $SLURM_CPUS_PER_TASK -b tmp.${LENGTH}.${PROFILENAME}.${O_N}.bam
    samtools depth tmp.${LENGTH}.${PROFILENAME}.${O_N}.bam | awk -v P="${PROFILENAME}" -v L="${LENGTH}" '{print P"\t"$0"\t"L}' > tmp.${LENGTH}.${PROFILENAME}.${O_N}.2
done
done

### with no filter for a later step
for SETS in $(cat tmp.${O_N}.1); do
    PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
    READS=$(echo $SETS | cut -d ";" -f 2)
    $seqkit seq -m 0 tmp.${PROFILENAME}.${O_N}.fa | minimap2 -t $SLURM_CPUS_PER_TASK -ax map-ont - $READS | samtools sort -O BAM - > tmp.0.${PROFILENAME}.${O_N}.bam
    samtools index -@ $SLURM_CPUS_PER_TASK -b tmp.0.${PROFILENAME}.${O_N}.bam
    samtools depth tmp.0.${PROFILENAME}.${O_N}.bam | awk -v P="${PROFILENAME}" '{print P"\t"$0}' | grep -f tmp.${O_N}.00 - > tmp.0.${PROFILENAME}.${O_N}.0
    mv tmp.0.${PROFILENAME}.${O_N}.bam ${PROFILENAME}.${O_N}.bam
done

### cocnatenate and filter down to only the contigs present in the bin
cat tmp.*.${O_N}.2 | grep -f tmp.${O_N}.00 - > ${O_N}_FINAL4.tsv


### cocnatenate and filter down to only the contigs present in the bin
cat tmp.0.*.${O_N}.0 | grep -f tmp.${O_N}.00 - > ${O_N}_FINAL6.tsv

rm -f tmp.*${O_N}*.bam tmp.*${O_N}*.bai tmp.*.${O_N}.2


conda activate /panfs/jay/groups/27/dcschroe/heske011/.conda/envs/R-env
Rscript /panfs/jay/groups/27/dcschroe/heske011/scripts/plot_postbin_mappingcheck.r

#####################################
fi
#####################################

#################### STEP A4 - Diamond blastx whole contigs against VOGDB (optional) #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A4" ]]; then
module purge; conda deactivate
conda activate $ANVIO

tmp="tmp.${O_N}"
rm -f ${tmp}.*

grep ">" ${bin}-contigs.fa | sed 's/>//g' | sort -Vu > ${tmp}.000

#### first get the taxonomy ids for this bin, and fetch potential viral reference genomes for this bin
cut -f 1 ${bin}-gene_calls.txt | sed -E 's/(.*)/^C\t\1\t/g' | grep -f - ${output}/ANVIO_${prefix}/${prefix}_gene_calls.rvdb_nr_euk.merged.kaiju | cut -f 2,3,8 > taxonomy-per-gene.txt

awk -F "\t" '{print "s/___"$1"___/"$2"/g"}' ${bin}-gene_calls.txt > ${tmp}.001

awk -F "\t" '{print "___"$1"___"}' taxonomy-per-gene.txt | sed -f ${tmp}.001 - | paste - taxonomy-per-gene.txt | sort -t $'\t' -V -k 1,2 | awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/[^a-zA-Z0-9;]/,"_",$4)}1' | sed -E 's/;_|_;/;/g' > ${tmp}.002 && mv ${tmp}.002 taxonomy-per-gene.txt

##### get 10 random genomes within taxonomic group identified
    echo -e "${taxonomy}" | sed 's/,/\n/g' > tmp.taxonomy
    for f in $(sed '/;NA;$/d' taxonomy-per-gene.txt | grep -f tmp.taxonomy - | cut -f 3 | sort -Vu | sort -VR | head -10); do
        ${EDIRECT}/esearch -db nuccore -query "txid${f}" | ${EDIRECT}/efetch -format fasta | $seqkit seq -j $SLURM_CPUS_PER_TASK -m $minsize -M $maxsize - | $seqkit rmdup -j $SLURM_CPUS_PER_TASK -s - | sed -E -e 's/[^>a-zA-Z0-9.]/_/g' -e 's/_+/_/g' -e 's/>/>reference_/g' > ${tmp}.${f}.003
    done


find ${tmp}.*.003 -type f -empty -delete
### here the taxa specific reference genomes are randomly shuffled
cat ${tmp}.*.003 | $seqkit shuffle > ${tmp}.006

#### without any specific keywords for genome provided, 3 random genomes are chosen, with priority for user selected taxa/genomes
if [ ! -z "$reference" ]; then
    echo -e "${reference}" | sed 's/,/\n/g' > tmp.reference

#### get fastas for user provided genomes on NCBI
    for f in $(cat tmp.reference); do
        ${EDIRECT}/esearch -db nuccore -query "${f}" | ${EDIRECT}/efetch -format fasta | sed -E -e 's/[^>a-zA-Z0-9.]/_/g' -e 's/_+/_/g' -e 's/>/>reference_/g'
    done > ${tmp}.003.ref.fa

    cat ${tmp}.003.ref.fa ${tmp}.006 | $seqkit rmdup -j $SLURM_CPUS_PER_TASK -n - > ${tmp}.007

#### select 3 reference genomes, with priority for user-provided ones, then followed by randomly ordered taxa specific ones

    grep ">" ${tmp}.007 | sed 's/>//g' | awk '!a[$0]++' > ${tmp}.003.ref

head -3 ${tmp}.003.ref | $seqkit grep -n -f - ${tmp}.007 | cat - ${bin}-contigs.fa > ${tmp}.005

fi

#### if user provides no references to fetch, then 3 random genomes are chosen based on the kaiju taxid and size range
if [ -z "$reference" ]; then

    rm -f ${tmp}.003.ref; touch ${tmp}.003.ref
    grep ">" ${tmp}.006 | sort -VR | sed 's/>//g' | awk '!a[$0]++' >> ${tmp}.003.ref
    head -3 ${tmp}.003.ref | $seqkit grep -n -f - ${tmp}.006 | cat - ${bin}-contigs.fa > ${tmp}.005
fi

### blastx whole splits against custom dmnd databases
### reduce to best hit per contig and VOG hit only (by evalue)
### set the max target sequences based on the longest contig, divided by 100 (at the shorter end of genes)
MAXTARGETS=$($seqkit stat -T ${tmp}.005 | sed '1,1d' | awk -F "\t" '{printf "%0.0f\n", $8/100}')

$DIAMOND blastx --max-target-seqs $MAXTARGETS --evalue 1e-120 --sensitive -p $SLURM_CPUS_PER_TASK -d ${CUSTOM_DMND}/vog.dmnd -q ${tmp}.005 -f 6 -o contigs-refs.dmnd.blastx

awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/\..*/,"",$2)}1' contigs-refs.dmnd.blastx | awk -F "\t" '!a[$1,$2]++' | cut -f 1,2,4,12 > ${tmp}.004

cut -f 2 ${tmp}.004 | sort -Vu > ${tmp}.008
grep ">" ${tmp}.005 | sort -Vu | sed 's/>//g' > ${tmp}.009

for f in $(cat ${tmp}.008); do awk -F "\t" -v F=${f} '{print $0"\t"F"\t0\t0"}' ${tmp}.009 ; done > ${tmp}.010

#### add everything back together, keeping the highest bitscore, which removes no hit results for those with hits 
cat ${tmp}.004 ${tmp}.010 | sort -t $'\t' -k 1,2V -k 3,3nr | awk -F "\t" '!a[$1,$2]++' | sed -z 's/^/contig\tVOG\tbitscore\talignlength\n/g' > contigs-refs-VOGs.tsv
##sed -i 's/\t0\t0$/\tNA\tNA/g' contigs-refs-VOGs.tsv

#### add VOG totals per query
cut -f 1,3 contigs-refs-VOGs.tsv | sed -e '1,1d' -e '/\t0$/d' | cut -f 1 | uniq -c | sed -E 's/ +([0-9]+) (.+)/\2\t\1/g' | awk '{print "s/"$1"/"$2"/g"}' > ${tmp}.011
cut -f 1,3 contigs-refs-VOGs.tsv | sed -e '1,1d' | awk '/\t0$/' | awk '{print "s/"$1"/0/g"}' >> ${tmp}.011
cut -f 1 contigs-refs-VOGs.tsv | sed -f ${tmp}.011 - | sed 's/^contig$/totalVOGs/g' | paste contigs-refs-VOGs.tsv - > ${tmp}.012 && mv ${tmp}.012 contigs-refs-VOGs.tsv

### reduce query and VOG names to max 50 characters
sed -i -E 's/^(..................................................).*\t(.*)\t(.*)\t(.*)\t(.*)/\1\t\2\t\3\t\4\t\5/g' contigs-refs-VOGs.tsv
sed -i -E 's/^(.*)\t(..................................................).*\t(.*)\t(.*)\t(.*)/\1\t\2\t\3\t\4\t\5/g' contigs-refs-VOGs.tsv

cut -f 1 contigs-refs-VOGs.tsv | sed '1,1d' | sort -Vu | wc -l > tmp
cut -f 2 contigs-refs-VOGs.tsv | sed '1,1d' | sort -Vu | wc -l | cat tmp - | sort -nr | head -1 > contigs-refs-VOGs.sc

conda activate /panfs/jay/groups/27/dcschroe/heske011/.conda/envs/R-env
Rscript /panfs/jay/groups/27/dcschroe/dmckeow/projects/DWV/script/postbin_heatmap.r

##############################################################
fi
##############################################################



########################################## STEP A5 - prepare table summarising properties
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A5" ]]; then
###########################################
module purge; conda deactivate


awk -F "\t" '{print $4 > "tmp."$2".a1"}' postbin_contigcheck_FINAL6.tsv

for f in tmp.*.a1; do
    awk -F "\t" '{ sum += $1; n++ } END { if (n > 0) printf "%0.0f\n", sum / n}' $f > $(basename $f .a1).a2
done

for f in tmp.*.a2; do awk -v f=$f '{print f"\t"$1}' $f | sed -E 's/tmp\.(.+)\.a2/\1/g'; done | sort -t $'\t' -k 1,1V | cut -f 2 > average_coverage_depth_percontig

sed -e '1,1d' -e '/^reference_/d' contigs-refs-VOGs.tsv | cut -f 1,5 | sort -Vu | sort -t $'\t' -k 1,1V | cut -f 2 > total_no_VOGs_percontig
sort -t $'\t' -k 1,1V -o postbin_contigcheck_FINAL3.tsv postbin_contigcheck_FINAL3.tsv

cut -f 1,2 postbin_contigcheck_FINAL1.tsv | sort -Vu | sort -t $'\t' -k 2,2V | paste - postbin_contigcheck_FINAL3.tsv average_coverage_depth_percontig total_no_VOGs_percontig | cut --complement -f 2 | sort -t $'\t' -k 4,4nr > postbin_contigcheck_FINAL5.tsv
sed -i "s/$/\t${bin}/g" postbin_contigcheck_FINAL5.tsv

  LENGTH=$((($minsize+$maxsize)/2))
sed -i "s/$/\t${LENGTH}/g" postbin_contigcheck_FINAL5.tsv

sed -z -i 's/^/sample\tcontig\tlength_bp\taverage_depth\tno_VOG_hits\tbin\treference_genome_length\n/1' postbin_contigcheck_FINAL5.tsv
sed -i 's/\t/;/g' postbin_contigcheck_FINAL5.tsv

#####################################
fi
#####################################

########################################## STEP A6 - alignment
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A6" ]]; then
###########################################
module purge; conda deactivate

module load mafft

mafft --adjustdirectionaccurately --thread $SLURM_CPUS_PER_TASK --auto ${tmp}.005 > postbin_contigcheck_FINAL5.aln
sed -i -E 's/^(>..................................................).*/\1/g' postbin_contigcheck_FINAL5.aln

#module load lastz/1.03.02
#lastz ${tmp}.005[multiple] ${tmp}.005[multiple] --notrivial --format=maf > postbin_contigcheck_FINAL5.maf

##conda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/R-env
##Rscript /panfs/jay/groups/27/dcschroe/heske011/scripts/plot_multiple_alignment.r

#####################################
fi
#####################################

########################################## STEP B1 - alignment of final curated contigs
if [[ "$step" =~ "B0" ]] || [[ "$step" =~ "B1" ]]; then
###########################################
module purge; conda deactivate
module load mafft
module load fasttree/2.1.8


awk -F ";" -v minlength=${minlength} -v mindepth=${mindepth} '{if($3 > minlength && $4 > mindepth && $0 !~ "sample;contig;length") print $0";yes;"minlength";"mindepth; else if($0 ~ "sample;contig;length") print $0";minlengthdepthpass;minlength;mindepth";  else print $0";no;"minlength";"mindepth}' postbin_contigcheck_FINAL5.tsv > tmp && mv tmp postbin_contigcheck_FINAL5.tsv

sed '1,1d' postbin_contigcheck_FINAL5.tsv | awk -F ";" -v minlength=${minlength} -v mindepth=${mindepth} '{if($3 > minlength && $4 > mindepth) print $2}' > ${tmp}.B001

grep ">reference_" ${tmp}.005 | sed 's/>//g' >> ${tmp}.B001

$seqkit grep -n -f ${tmp}.B001 ${tmp}.005 > ${bin}-d${mindepth}-l${minlength}-finalcontigs.fa
sed -i -E 's/^(>..................................................).*/\1/g' ${bin}-d${mindepth}-l${minlength}-finalcontigs.fa


mafft --adjustdirectionaccurately --thread $SLURM_CPUS_PER_TASK --auto ${bin}-d${mindepth}-l${minlength}-finalcontigs.fa > ${bin}-d${mindepth}-l${minlength}-finalcontigs.aln

FastTree -gtr -nt < ${bin}-d${mindepth}-l${minlength}-finalcontigs.aln > ${bin}-d${mindepth}-l${minlength}-finalcontigs.fasttree


#####################################
fi
#####################################


########################################## STEP B2 - HYPHY GARD
if [[ "$step" =~ "B0" ]] || [[ "$step" =~ "B2" ]]; then
###########################################
module purge; conda deactivate
module load hyphy/2.5.33


hyphy CPU=$SLURM_CPUS_PER_TASK gard --alignment ${bin}-d${mindepth}-l${minlength}-finalcontigs.aln

#####################################
fi
#####################################

########################################## STEP test -convert annotations to the actual gff format
if [[ "$step" =~ "B0" ]] || [[ "$step" =~ "B3" ]]; then
###########################################

awk -F "\t" 'BEGIN{OFS="\t"};{print $2,"prodigal","CDS",$3,$4,".",$5,".","ID="$1}' ${bin}-gene_calls.txt | sed -e 's/\tr\t/\t-\t/g' -e 's/\tf\t/\t+\t/g' -e '/^contig\tprodigal\tCDS\tstart/d' > ${bin}-gene_calls.tmp1

cut --complement -f 1-5 ${bin}-gene_calls.txt > ${bin}-gene_calls.tmp2

head -1 ${bin}-gene_calls.tmp2 | tr "\t" "\n" > ${bin}-gene_calls.tmp3
sed -i -E -e 's/[^a-zA-Z0-9_]/_/g' -e 's/_+/_/g' -e 's/^_|_$//g' ${bin}-gene_calls.tmp3

sed -i '1,1d' ${bin}-gene_calls.tmp2

rm -f ${bin}-*.tmp4

id=0
while read filename; do
  id=$((id+1))
 tail -n +1 ${bin}-gene_calls.tmp2 | cut -f $id >> ${bin}-"$filename".tmp4
done < ${bin}-gene_calls.tmp3

rm -f ${bin}-dna_sequence.tmp4

sed -i 's/^$/-/g' ${bin}-*.tmp4

for f in ${bin}-*.tmp4; do
    awk -v f=$(basename $f .tmp4) '{print f"="$0}' $f > $(basename $f .tmp4).tmp5
done

paste -d ";" ${bin}-gene_calls.tmp1 ${bin}-*.tmp5 > ${bin}-gene_calls.gff

rm -f ${bin}-*.tmp4

#####################################
fi
#####################################