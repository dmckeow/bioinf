#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --mem=48GB
#SBATCH --cpus-per-task=12
#SBATCH --partition ag2tb
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

echo -e "COMMAND LINE JOB SUBMISSSION:\n\tsbatch /home/dcschroe/dmckeow/projects/DWV/script/amplicont.sh $@"

Help()
{
echo -e "########### REQUIRED arguments ###########\n-r, --reads\t file listing the full path(s) to fastq.gz files OR the absolute path to a single folder or file. Can be multiple or single, and files or folders. MUST be .fastq.gz"
echo -e "-f, --fasta\t reference sequence to map reads against - only used to remove non-target reads - none of these sequence processing steps are reference-based"
echo -e "-s, --step\tWhich script steps to run: A1 prepares and polishes reads, A2 extracts and clips are final reads"
echo -e "-o, --output\tfull path for final output directory\n"
echo -e "-p, --prefix\tprefix name for output files \n"
echo -e "-l, --length\tthe length for your amplicon. This will be used as a cutoff value to remove shorter or partial amplicons\n"

echo -e "-h, --help\tshow this help message and exit\n"
}

while getopts r:f:s:o:p:l:h option
do 
    case "${option}" in
        r)reads=${OPTARG};;
        f)fasta=${OPTARG};;
        s)step=${OPTARG};;
    o)output=${OPTARG};;
    p)prefix=${OPTARG};;
    l)length=${OPTARG};;

    h)Help; exit;;
    esac
done

if [[ -z "${reads}" ]]; then echo "-r, --reads REQUIRED"; Help; exit; fi
if [[ -z "${fasta}" ]]; then echo "-f, --fasta REQUIRED"; Help; exit; fi
if [[ -z "${step}" ]]; then echo "-s, --step REQUIRED"; Help; exit; fi
if [[ -z "${output}" ]]; then echo "-o, --output REQUIRED"; Help; exit; fi
if [[ -z "${prefix}" ]]; then echo "-p, --prefix REQUIRED"; Help; exit; fi
if [[ -z "${length}" ]]; then echo "-l, --length REQUIRED"; Help; exit; fi


####################### PREREQUISITES #####################################

####################### SCRIPT PURPOSE #####################################

####################### HOW TO RUN SCRIPT #####################################

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################

####################### OUTPUT FILES #####################################

####################### LOAD SOFTWARE #####################################
module purge
eval "$(conda shell.bash hook)"
conda deactivate

SEQKIT="/home/dcschroe/dmckeow/seqkit"
####MINIASM="/panfs/jay/groups/27/dcschroe/dmckeow/miniasm/miniasm"

####################### SET VARIABLES/ARGUMENTS #####################################


################################################################################################
####################### RUNNING THE SCRIPT #####################################################
################################################################################################

##############################################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A1" ]]; then
##############################################################
module purge; conda deactivate
module load minimap2/2.17
module load samtools
module load porechop
module load racon/1.4.20

### single/multiple fastq.gz file(s) OR single/multiple directories containing fastq.gz file(s)
mkdir ${output}
rm -fr ${output}/amplicont-${prefix}
mkdir ${output}/amplicont-${prefix}
touch ${output}/amplicont-${prefix}/${prefix}.bctrimmedreads.fastq.tmp
touch ${output}/amplicont-${prefix}/${prefix}.bctrimmedreads.list

######### concatenate reads into a single file
for f in $(cat $reads)
do
    find $f -name "*.fastq.gz" >> ${output}/amplicont-${prefix}/${prefix}.bctrimmedreads.list
done

#### check if a folder or file path was given instead of readlist file
i1=$(grep -s "*.fastq.gz" ${output}/amplicont-${prefix}/${prefix}.bctrimmedreads.list | wc -l)
if [[ $i1 -le 0 ]]; then
    find $reads -name "*.fastq.gz" > ${output}/amplicont-${prefix}/${prefix}.bctrimmedreads.list
fi

for f in $(cat ${output}/amplicont-${prefix}/${prefix}.bctrimmedreads.list)
do
    zcat $f >>  ${output}/amplicont-${prefix}/${prefix}.bctrimmedreads.fastq.tmp
done

#### remove barcodes
porechop -t $SLURM_CPUS_PER_TASK -i ${output}/amplicont-${prefix}/${prefix}.bctrimmedreads.fastq.tmp > ${output}/amplicont-${prefix}/${prefix}-trimmedbcs.fastq


#### polish ALL reads
minimap2 -t $SLURM_CPUS_PER_TASK -ax ava-ont ${output}/amplicont-${prefix}/${prefix}-trimmedbcs.fastq ${output}/amplicont-${prefix}/${prefix}-trimmedbcs.fastq > ${output}/amplicont-${prefix}/${prefix}-trimmedbcs.sam

racon -f ${output}/amplicont-${prefix}/${prefix}-trimmedbcs.fastq ${output}/amplicont-${prefix}/${prefix}-trimmedbcs.sam ${output}/amplicont-${prefix}/${prefix}-trimmedbcs.fastq -t $SLURM_CPUS_PER_TASK -q 9 > ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads.fa

rm -f ${output}/amplicont-${prefix}/${prefix}-trimmedbcs.sam

##############################################################
fi
##############################################################

##############################################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A2" ]]; then
##############################################################

module purge; conda deactivate
module load minimap2/2.17
module load samtools
conda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/fgbio
module load cd-hit

##### pull out only reads with that align to reference
minimap2 -t $SLURM_CPUS_PER_TASK --sam-hit-only --secondary=no -ax map-ont $fasta ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads.fa > ${output}/amplicont-${prefix}/${prefix}.sam
samtools sort -@ $SLURM_CPUS_PER_TASK -O BAM ${output}/amplicont-${prefix}/${prefix}.sam > ${output}/amplicont-${prefix}/${prefix}.bam
rm -f ${output}/amplicont-${prefix}/${prefix}.sam
samtools index -@ $SLURM_CPUS_PER_TASK -b ${output}/amplicont-${prefix}/${prefix}.bam

##### convert soft clips of reads into hard clips (bad sequences at the ends are removed)
samtools sort -nu -@ $SLURM_CPUS_PER_TASK -O BAM ${output}/amplicont-${prefix}/${prefix}.bam > ${output}/amplicont-${prefix}/${prefix}-tmp.bam
fgbio ClipBam -H -i ${output}/amplicont-${prefix}/${prefix}-tmp.bam -r $fasta -o ${output}/amplicont-${prefix}/${prefix}-hardclip.bam

samtools sort -@ $SLURM_CPUS_PER_TASK -O BAM ${output}/amplicont-${prefix}/${prefix}-hardclip.bam > ${output}/amplicont-${prefix}/${prefix}-hardclip-tmp.bam
mv ${output}/amplicont-${prefix}/${prefix}-hardclip-tmp.bam ${output}/amplicont-${prefix}/${prefix}-hardclip.bam
samtools index -@ $SLURM_CPUS_PER_TASK -b ${output}/amplicont-${prefix}/${prefix}-hardclip.bam

##### extract hardclipped fastas AND remove 100% identical sequences
samtools fasta -@ $SLURM_CPUS_PER_TASK ${output}/amplicont-${prefix}/${prefix}-hardclip.bam | $SEQKIT rmdup -s - > ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip.fa

##### remove shorter sequences that are 100 % overlapped by a longer one
cd-hit-est -c 1.0 -i ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip.fa -o ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup.fa

#### remove sequences are shorter than your amplicon length threshold
$SEQKIT seq -m ${length} ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup.fa > ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup-tmp.fa
mv ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup-tmp.fa ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup.fa

rm -f ${output}/amplicont-${prefix}/${prefix}-tmp.bam

##############################################################
fi
##############################################################


##############################################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A3" ]]; then
##############################################################
module purge; conda deactivate
module load mafft
module load fasttree/2.1.8

cat $fasta ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup.fa > ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup-tmp.fa

mafft --adjustdirectionaccurately --thread $SLURM_CPUS_PER_TASK --reorder --auto ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup-tmp.fa > ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup.aln

FastTree -gtr -nt < ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup.aln > ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup.fasttree

##############################################################
fi
##############################################################

##############################################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A4" ]]; then
##############################################################
module purge; conda deactivate
module load cd-hit

##### remove shorter sequences that are overlapped by a longer one by various %s
cd-hit-est -c 0.95 -i ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup.fa -o ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup-95.fa
cd-hit-est -c 0.9 -i ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup.fa -o ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup-90.fa
cd-hit-est -c 0.85 -i ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup.fa -o ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup-85.fa
cd-hit-est -c 0.80 -i ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup.fa -o ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup-80.fa

sed -i "s/>/>${prefix}__/g" ${output}/amplicont-${prefix}/${prefix}-trimmedbcs-polishedreads-mapped-hardclip-dedup*.fa



##############################################################
fi
##############################################################

##############################################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A5" ]]; then
##############################################################
module purge; conda deactivate
module load mafft
module load fasttree/2.1.8

mkdir ${output}/SUMMARY

a=$(find . -name "*-dedup.fa")
b=$(find . -name "*-dedup-95.fa")
c=$(find . -name "*-dedup-90.fa")
d=$(find . -name "*-dedup-85.fa")
e=$(find . -name "*-dedup-80.fa")

cat $fasta $a > ${output}/SUMMARY/ALLSAMPLES_trimmedbcs-polishedreads-mapped-hardclip-dedup.fa
cat $fasta $b > ${output}/SUMMARY/ALLSAMPLES_trimmedbcs-polishedreads-mapped-hardclip-dedup-95.fa
cat $fasta $c > ${output}/SUMMARY/ALLSAMPLES_trimmedbcs-polishedreads-mapped-hardclip-dedup-90.fa
cat $fasta $d > ${output}/SUMMARY/ALLSAMPLES_trimmedbcs-polishedreads-mapped-hardclip-dedup-85.fa
cat $fasta $e > ${output}/SUMMARY/ALLSAMPLES_trimmedbcs-polishedreads-mapped-hardclip-dedup-80.fa

rm -f ${output}/SUMMARY/stats_all_fastas.tsv; touch ${output}/SUMMARY/stats_all_fastas.tsv

for f in ${output}/SUMMARY/ALLSAMPLES_trimmedbcs-polishedreads-mapped-hardclip-dedup*.fa; do
$SEQKIT stat -T $f >> ${output}/SUMMARY/stats_all_fastas.tsv
mafft --adjustdirectionaccurately --thread $SLURM_CPUS_PER_TASK --reorder --auto $f > ${f}.aln

FastTree -gtr -nt < ${f}.aln > ${f}.fasttree
done 
##############################################################
fi
##############################################################

