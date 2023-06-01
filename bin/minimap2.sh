#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8
#SBATCH --partition agsmall
#SBATCH -o /panfs/jay/groups/27/dcschroe/dmckeow/logs/slurm.%N.%j.out
#SBATCH -e /panfs/jay/groups/27/dcschroe/dmckeow/logs/slurm.%N.%j.err

echo -e "COMMAND LINE JOB SUBMISSSION:\n\tsbatch /home/dcschroe/dmckeow/projects/DWV/script/minimap2.sh $@"

Help()
{
echo -e "########### REQUIRED arguments ###########\n-r, --reads\t file listing the full path(s) to fastq.gz files OR the absolute path to a single folder or file. Can be multiple or single, and files or folders. MUST be .fastq.gz"
echo -e "-f, --fasta\t reference sequence to map reads against "
echo -e "-s, --step\tWhich script steps to run: A1 with porechop, A2 no porechop"
echo -e "-o, --output\tfull path for final output directory\n"
echo -e "-p, --prefix\tprefix name for output files \n"

echo -e "-h, --help\tshow this help message and exit\n"
}

while getopts r:f:s:o:p:h option
do 
    case "${option}" in
        r)reads=${OPTARG};;
        f)fasta=${OPTARG};;
        s)step=${OPTARG};;
    o)output=${OPTARG};;
    p)prefix=${OPTARG};;
    h)Help; exit;;
    esac
done

if [[ -z "${reads}" ]]; then echo "-r, --reads REQUIRED"; Help; exit; fi
if [[ -z "${fasta}" ]]; then echo "-f, --fasta REQUIRED"; Help; exit; fi
if [[ -z "${step}" ]]; then echo "-s, --step REQUIRED"; Help; exit; fi
if [[ -z "${output}" ]]; then echo "-o, --output REQUIRED"; Help; exit; fi
if [[ -z "${prefix}" ]]; then echo "-p, --prefix REQUIRED"; Help; exit; fi

####################### PREREQUISITES #####################################
####################### SCRIPT PURPOSE #####################################

####################### HOW TO RUN SCRIPT #####################################

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################

####################### OUTPUT FILES #####################################

####################### LOAD SOFTWARE #####################################
module purge
eval "$(conda shell.bash hook)"
module load minimap2/2.17
module load samtools
module load pigz
module load porechop
SEQKIT="/home/dcschroe/dmckeow/seqkit"

####################### SET VARIABLES/ARGUMENTS #####################################
################################################################################################
####################### RUNNING THE SCRIPT #####################################################
################################################################################################
if [[ "$step" =~ "A1" ]]; then

### single/multiple fastq.gz file(s) OR single/multiple directories containing fastq.gz file(s)

rm -fr ${output}/minimap${prefix}
mkdir ${output}/minimap${prefix}
touch ${output}/minimap${prefix}/${prefix}.bctrimmedreads.fastq.tmp
touch ${output}/minimap${prefix}/${prefix}.bctrimmedreads.list

for f in $(cat $reads)
do
    find $f -name "*.fastq.gz" >> ${output}/minimap${prefix}/${prefix}.bctrimmedreads.list
done

#### check if a folder or file path was given instead of readlist file
i1=$(grep -s "*.fastq.gz" ${output}/minimap${prefix}/${prefix}.bctrimmedreads.list | wc -l)
if [[ $i1 -le 0 ]]; then
    find $reads -name "*.fastq.gz" > ${output}/minimap${prefix}/${prefix}.bctrimmedreads.list
fi

for f in $(cat ${output}/minimap${prefix}/${prefix}.bctrimmedreads.list)
do
    zcat $f >>  ${output}/minimap${prefix}/${prefix}.bctrimmedreads.fastq.tmp
done

porechop -t $SLURM_CPUS_PER_TASK -i  ${output}/minimap${prefix}/${prefix}.bctrimmedreads.fastq.tmp -o ${output}/minimap${prefix}/${prefix}.bctrimmedreads.fastq
rm -f ${output}/minimap${prefix}/${prefix}.bctrimmedreads.fastq.tmp
pigz -p $SLURM_CPUS_PER_TASK ${output}/minimap${prefix}/${prefix}.bctrimmedreads.fastq


#### mapping, including sorting and indexing of bam alignment
minimap2 -t $SLURM_CPUS_PER_TASK --sam-hit-only -ax map-ont $fasta ${output}/minimap${prefix}/${prefix}.bctrimmedreads.fastq.gz > ${output}/minimap${prefix}/${prefix}.sam
samtools sort -@ $SLURM_CPUS_PER_TASK -O BAM ${output}/minimap${prefix}/${prefix}.sam > ${output}/minimap${prefix}/${prefix}.bam
rm -f ${output}/minimap${prefix}/${prefix}.sam
samtools index -@ $SLURM_CPUS_PER_TASK -b ${output}/minimap${prefix}/${prefix}.bam

fi

##############################################################
##############################################################

################################################################################################
if [[ "$step" =~ "A2" ]]; then

### single/multiple fastq.gz file(s) OR single/multiple directories containing fastq.gz file(s)

rm -fr ${output}/minimap${prefix}
mkdir ${output}/minimap${prefix}
touch ${output}/minimap${prefix}/${prefix}.bctrimmedreads.fastq.tmp
touch ${output}/minimap${prefix}/${prefix}.bctrimmedreads.list

for f in $(cat $reads)
do
    find $f -name "*.fastq.gz" >> ${output}/minimap${prefix}/${prefix}.bctrimmedreads.list
done

#### check if a folder or file path was given instead of readlist file
i1=$(grep -s "*.fastq.gz" ${output}/minimap${prefix}/${prefix}.bctrimmedreads.list | wc -l)
if [[ $i1 -le 0 ]]; then
    find $reads -name "*.fastq.gz" > ${output}/minimap${prefix}/${prefix}.bctrimmedreads.list
fi

for f in $(cat ${output}/minimap${prefix}/${prefix}.bctrimmedreads.list)
do
    zcat $f >>  ${output}/minimap${prefix}/${prefix}.bctrimmedreads.fastq.tmp
done

mv ${output}/minimap${prefix}/${prefix}.bctrimmedreads.fastq.tmp ${output}/minimap${prefix}/${prefix}.bctrimmedreads.fastq
rm -f ${output}/minimap${prefix}/${prefix}.bctrimmedreads.fastq.tmp
pigz -p $SLURM_CPUS_PER_TASK ${output}/minimap${prefix}/${prefix}.bctrimmedreads.fastq


#### mapping, including sorting and indexing of bam alignment
minimap2 -t $SLURM_CPUS_PER_TASK --sam-hit-only -ax map-ont $fasta ${output}/minimap${prefix}/${prefix}.bctrimmedreads.fastq.gz > ${output}/minimap${prefix}/${prefix}.sam
samtools sort -@ $SLURM_CPUS_PER_TASK -O BAM ${output}/minimap${prefix}/${prefix}.sam > ${output}/minimap${prefix}/${prefix}.bam
rm -f ${output}/minimap${prefix}/${prefix}.sam
samtools index -@ $SLURM_CPUS_PER_TASK -b ${output}/minimap${prefix}/${prefix}.bam

fi

##############################################################
##############################################################