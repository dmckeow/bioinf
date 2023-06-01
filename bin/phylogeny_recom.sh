#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=120GB
#SBATCH --partition ag2tb
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err


Help()
{
echo -e "########### REQUIRED arguments ###########\n-i, --input\tfull path to fasta file"
echo -e "########### OPTIONAL arguments ###########\n"
echo -e "-s, --step\tWhich script steps to run. leave blank = run the whole script. A1 is mafft alignment, A2 is hyphy-gard, A3 is phylogeny"
echo -e "-h, --help\tshow this help message and exit\n"
}

while getopts i:s:h option
do 
    case "${option}" in
        i)input=${OPTARG};;
        s)step=${OPTARG};;
    h)Help; exit;;
    esac
done

if [[ -z "${input}" ]]; then echo "-i, --input REQUIRED"; Help; exit; fi

####################### PREREQUISITES #####################################
## concatenated fasta of sequences for analyses
####################### SCRIPT PURPOSE #####################################

####################### HOW TO RUN SCRIPT #####################################

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################

####################### OUTPUT FILES #####################################

####################### SOFTWARE, DATABSES, ETC #####################################
### LOAD software available via shared environment on server:
module purge
eval "$(conda shell.bash hook)"

########################################################################
########################################################################

######### SOFTWARE/OTHER SCRIPTS THAT NEED SETUP BEFORE RUNNING:

SEQKIT="/home/dcschroe/dmckeow/seqkit" ## seqkit - installed locally


####################### THE SCRIPT #####################################
S_N="phylogeny_recom_sh" ## script name
O_N=$(echo ${input} | sed 's/\./_/g') ## output name
C_S="A1" ## current step name
T_N="${O_N}-tmp-${S_N}.${C_S}" ## temporary file name
F_N="${O_N}-final-${S_N}" ## final output file name

################################################################################################
################################################################################################
########################################################################
########################################################################

#################### STEP A1 - alignment #################################
C_S="A1"
if [[ -z "$step" ]] || [[ "$step" =~ "$C_S" ]]; then
  module purge
  conda deactivate
  module load mafft

mafft --adjustdirectionaccurately --thread $SLURM_CPUS_PER_TASK --reorder --auto $input > ${F_N}.aln

fi ########################################################################


#################### STEP A2 - fast phylogeny #################################
C_S="A2"
if [[ -z "$step" ]] || [[ "$step" =~ "$C_S" ]]; then
  module purge
  conda deactivate
  module load fasttree/2.1.8

  if [[ -f "${F_N}.aln" ]]; then
    FastTree ${F_N}.aln > ${F_N}.fasttree
  else
    FastTree ${input}.aln > ${F_N}.fasttree
  fi

fi ########################################################################



#################### STEP A3 - slow phylogeny #################################
C_S="A3"
if [[ -z "$step" ]] || [[ "$step" =~ "$C_S" ]]; then
  module purge
  conda deactivate
  module load raxml/8.2.9_pthread

raxmlHPC-PTHREADS -T $SLURM_CPUS_PER_TASK -d -f ae -m PROTGAMMAAUTO -p $RANDOM -x $RANDOM -N 200 -s ${input}.aln -n ${input}.200bs.raxml

fi ########################################################################


#################### STEP A4 - hyphy GARD analyses #################################
C_S="A4"
if [[ -z "$step" ]] || [[ "$step" =~ "$C_S" ]]; then
  
module load hyphy/2.5.33

  if [[ -f "${F_N}.aln" ]]; then
    hyphy CPU=$SLURM_CPUS_PER_TASK gard --alignment ${F_N}.aln
  else
    hyphy CPU=$SLURM_CPUS_PER_TASK gard --alignment $input
  fi

fi ########################################################################