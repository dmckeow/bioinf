#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=150GB
#SBATCH --partition agsmall
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

Help()
{
echo -e "########### REQUIRED arguments ###########\n-i, --input\tfull path to directory containing reads fastq.gz files (e.g. the barcode folder)"
echo -e "-s, --step\tWhich script steps to run. A0 = run the whole script. "
echo -e "-o, --output\tfull path for final output directory\n"
echo -e "-t, --tmp\tfull path for directory in which a temporary output directory will be created. I recommend you do it in scratch storage, e.g.: /scratch.global/dcschroe\n"
echo -e "-p, --prefix\tprefix name for output files - e.g. sequencingrun11OCT22_barcode5\n"

echo -e "-h, --help\tshow this help message and exit\n"
}

while getopts i:s:o:t:p:h option
do 
    case "${option}" in
        i)input=${OPTARG};;
        s)step=${OPTARG};;
    o)output=${OPTARG};;
    t)tmp=${OPTARG};;
    p)prefix=${OPTARG};;
    h)Help; exit;;
    esac
done

if [[ -z "${input}" ]]; then echo "-i, --input REQUIRED"; Help; exit; fi
if [[ -z "${step}" ]]; then echo "-s, --step REQUIRED"; Help; exit; fi
if [[ -z "${output}" ]]; then echo "-o, --output REQUIRED"; Help; exit; fi
if [[ -z "${tmp}" ]]; then echo "-t, --tmp REQUIRED"; Help; exit; fi
if [[ -z "${prefix}" ]]; then echo "-p, --prefix REQUIRED"; Help; exit; fi



####################### PREREQUISITES #####################################
####################### SCRIPT PURPOSE #####################################

####################### HOW TO RUN SCRIPT #####################################

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################

####################### OUTPUT FILES #####################################

####################### LOAD SOFTWARE #####################################
module purge
eval "$(conda shell.bash hook)"
#conda activate flye
FLYE="python /panfs/jay/groups/27/dcschroe/dmckeow/data/Flye/bin/flye"
module load minimap2/2.17
module load samtools
module load porechop/0.2.4
module load pigz

####################### SET VARIABLES/ARGUMENTS #####################################

################################################################################################
####################### RUNNING THE SCRIPT #####################################################
################################################################################################

################################ STEP 001 ###########################################

######### Prepare the reads for assembly

if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A1" ]]; then

rm -fr ${tmp}/tmp_flye_assembly_${prefix}
mkdir ${tmp}/tmp_flye_assembly_${prefix}
touch ${tmp}/tmp_flye_assembly_${prefix}/tmp_allreads_${prefix}.fastq
find $input -name "*.fastq.gz" >  ${tmp}/tmp_flye_assembly_${prefix}/tmp_allreads_${prefix}.list
sed -i '/\/unclassified\//d'  ${tmp}/tmp_flye_assembly_${prefix}/tmp_allreads_${prefix}.list

for f in $(cat ${tmp}/tmp_flye_assembly_${prefix}/tmp_allreads_${prefix}.list)
do
    zcat $f >>  ${tmp}/tmp_flye_assembly_${prefix}/tmp_allreads_${prefix}.fastq
done

porechop -t $SLURM_CPUS_PER_TASK -i  ${tmp}/tmp_flye_assembly_${prefix}/tmp_allreads_${prefix}.fastq -o  ${tmp}/tmp_flye_assembly_${prefix}/tmp
$SEQKIT seq -j $SLURM_CPUS_PER_TASK -m 50  ${tmp}/tmp_flye_assembly_${prefix}/tmp >  ${tmp}/tmp_flye_assembly_${prefix}/tmp_allreads_${prefix}.fastq
rm -f  ${tmp}/tmp_flye_assembly_${prefix}/tmp
pigz -p $SLURM_CPUS_PER_TASK  ${tmp}/tmp_flye_assembly_${prefix}/tmp_allreads_${prefix}.fastq

fi

################################ STEP 002 ###########################################

####### DE NOVO ASSEMBLY
##### with flag A2, you will rnu only this step, which means CANU will restart previous run attempt, without losing any progress
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A2" ]]; then

#READMIN=$($SEQKIT stat -T -j $SLURM_CPUS_PER_TASK  ${tmp}/tmp_flye_assembly_${prefix}/tmp_allreads_${prefix}.fastq.gz | cut -f 7 | sed "1,1d" | sed -E "s/(.+)/\1-49/g" | bc | sed -E "s/([0-9]+)\.[0-9]+/\1/g")

#### run assembly
rm -fr ${output}/$prefix
mkdir ${output}/$prefix

$FLYE --nano-raw  ${tmp}/tmp_flye_assembly_${prefix}/tmp_allreads_${prefix}.fastq.gz --out-dir  ${tmp}/tmp_flye_assembly_${prefix} --meta -t $SLURM_CPUS_PER_TASK -m 50

#cp  ${tmp}/tmp_flye_assembly_${prefix}/tmp_allreads_${prefix}.fastq.gz ${output}/$prefix/tmp_allreads_${prefix}.fastq.gz

fi


