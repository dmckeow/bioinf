#!/usr/bin/env python3
# -*- coding: utf-8

####################### PYTHON IMPORTS ###########################

import argparse
import os
import subprocess
import shlex
import sys
import pysam

####################### PARSE THE ARGUMENTS AND INPUT FLAGS ###########################

# create the argparse parser
parser = argparse.ArgumentParser()

# add the arguments
parser.add_argument("-i", "--input", required=True, 
                    help="absolute or relative path to a fasta or alignment file")
parser.add_argument("-s", "--step", 
                    help="Which script steps to run. ", 
                    default="A0")
parser.add_argument("-s", "--suffix", 
                    help="A meaningful name for output folder. Any pre-existing driectory at this location and identical name will be deleted. Will also be the name of the merged analyses files. It should be a UNIQUE name that cannot possibly be identical to any other folder within your output destination e.g. -p 11OCT22project will make the folder ANVIO_11OCT22project", default=".canu-read-fate"suffix)
parser.add_argument("-t", "--threads", 
                    help="\tNumber of threads (default: 12)", default=12)

# parse the arguments
args = parser.parse_args()

input = args.input
step = args.step
suffix = args.suffix
threads = args.threads

####################### VARIABLES ###########################
path = os.getcwd() ## set output directory

####################### SOFTWARE ###########################

# KAIJU databases
kaiju_rvdb = "/panfs/jay/groups/27/dcschroe/shared/tools/kaijudb/rvdb"
kaiju_nr_euk = "/panfs/jay/groups/27/dcschroe/shared/tools/kaijudb/nr_euk"

##### Seqkit - must be installed locally, e.g.:
SEQKIT="/panfs/jay/groups/27/dcschroe/shared/tools/seqkit"

##### Diamond - installed locally
DIAMOND="/panfs/jay/groups/27/dcschroe/shared/tools/diamond"

##### Kaiju - INSTALLED via CONDA
KAIJU="/panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/kaiju"

####################### DATABASES ###########################



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