#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8
#SBATCH --partition agsmall
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

echo -e "COMMAND LINE JOB SUBMISSSION:\n\tsbatch /home/dcschroe/dmckeow/projects/DWV/script/minimap2_simple.sh $@"

Help()
{
echo -e "########### REQUIRED arguments ###########\n-r, --reads\t path to a single file"
echo -e "-f, --fasta\t reference sequence to map reads against "
echo -e "-o, --output\tname for output files\n"

echo -e "-h, --help\tshow this help message and exit\n"
}

while getopts r:f:o:h option
do 
    case "${option}" in
        r)reads=${OPTARG};;
        f)fasta=${OPTARG};;
    o)output=${OPTARG};;
    h)Help; exit;;
    esac
done

if [[ -z "${reads}" ]]; then echo "-r, --reads REQUIRED"; Help; exit; fi
if [[ -z "${fasta}" ]]; then echo "-f, --fasta REQUIRED"; Help; exit; fi
if [[ -z "${output}" ]]; then echo "-o, --output REQUIRED"; Help; exit; fi

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

####################### SET VARIABLES/ARGUMENTS #####################################
################################################################################################
####################### RUNNING THE SCRIPT #####################################################
################################################################################################

### single/multiple fastq.gz file(s) OR single/multiple directories containing fastq.gz file(s)

#### mapping, including sorting and indexing of bam alignment
minimap2 -t $SLURM_CPUS_PER_TASK --sam-hit-only -ax map-ont $fasta $reads > ${output}.sam
samtools sort -@ $SLURM_CPUS_PER_TASK -O BAM ${output}.sam > ${output}.bam
rm -f ${output}.sam
samtools index -@ $SLURM_CPUS_PER_TASK -b ${output}.bam


##############################################################
##############################################################

