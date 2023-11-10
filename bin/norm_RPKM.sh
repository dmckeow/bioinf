#!/bin/bash

umask 002
green='\033[32m';red='\033[31m';cyan='\033[36m';nocolor='\033[m'
echo -e "${green}RUNNING bioinf-binning with the following parameters: $@ ${nocolor}"

Help()
{
echo -e "Bioinf-binning is a pipeline that uses anvi'o to perform binning on sequence datasets using custom databases setup by the user. It is intended to run on contigs from an assembly, and to maps the reads against this assembly. It can also be run without reads. Please note that currently the default mapping uses minimap and expects Oxford Nanopore reads. This script is designed to run on a single project per job submission, i.e. every sample in your input will be analysed together and end up in a single anvi'o figure/analyses\n"
echo -e "${green}########### REQUIRED parameters ###########${nocolor}\n-i --input\tAbsolute path to a plain text file with 3 columns separated by semicolons or tabs:\n"
echo -e "\tCOLUMN 1: /absolute_path/to/each/contig_fasta - separate line per assembly"
echo -e "\tCOLUMN 2: /absolute_path/to/each/reads_file - separate line per assembly, on same line with corresponding contigs"
echo -e "\tCOLUMN 3: profile name; a short name for each profile which will appear in all outputs and figures ${red} MUST NOT begin with a digit, MUST ONLY include letters, digits, and underscores, and MUST be unique${nocolor}"
echo -e "\t--input file example:\n\t\t${green}/home/data/19AUG22_1_DM_barcode01.contigs.fasta;/home/data/19AUG22_1_DM_barcode01-trimmedbcs.fastq.gz;location_1\n\t\t/home/data/19AUG22_1_DM_barcode03.contigs.fasta;/home/data/19AUG22_1_DM_barcode03-trimmedbcs.fastq.gz;location_2\n\t\t/home/data/19AUG22_1_DM_barcode06.contigs.fasta;/home/data/19AUG22_1_DM_barcode06-trimmedbcs.fastq.gz;location_3${nocolor}"
echo -e "\n-p --project\tA meaningful and unique name for your project. Output files will contain this name\n"
echo -e "${cyan}########### OPTIONAL parameters ###########${nocolor}"
echo -e "\n-m --mincontigsize\tinteger to set minimum size of contigs in bp to include in analyses (don't go below 1000)"
echo -e "-s --step\tWhich script steps to run. If not provided, the whole script is run. \nMultiple steps can be specified, e.g.: A1_A2_A3 will only run steps A1-A3.\n\t${green}STEPS AVAILABLE:${nocolor}"
awk '/^###### STEP-/ {print "\t\t"$0}' $(which bioinf-binning.sh)
echo -e "\n-n, --nomap [optional]\trun without mapping reads to contigs - allows you to run pipeline if you only have an assembly but no reads\n"
echo -e "-S --splitlength\tinteger to set minimum size of split in bp (default is 20000)\n"
echo -e "\n${green}example job submissions:${nocolor}\n"
echo -e "\tSUBMIT TO SLURM (see README for more info): ${cyan}sbatch --time=96:00:00 --cpus-per-task=24 --mem=240GB --partition long -o slurm.%N.%j.out -e slurm.%N.%j.err bioinf-binning.sh -i /path/to/samplelistfile -p project_name${nocolor}"
echo -e "\tSUBMIT TO SLURM, running only steps A1 to A3, and without mapping: ${cyan}sbatch --time=96:00:00 --cpus-per-task=24 --mem=240GB --partition long -o slurm.%N.%j.out -e slurm.%N.%j.err bioinf-binning.sh -i /path/to/samplelistfile -p project_name -s A1_A2_A3 -n${nocolor}"
echo -e "\tRUN LOCALLY (not recommended): ${cyan}bioinf-binning.sh -i /path/to/samplelistfile -p project_name --threads 24${nocolor}"
echo -e "\n-C, --concoct\trun auto binning with CONCOCT. Requires this script be run with mapping. Requires concoct setup through anvi'o as detailed in anvi'o installation\n"
echo -e "-t, --threads\tnumber of threads for job (SLURM --cpus-per-task does the same thing); default 1\n"
echo -e "-h, --help\tshow this help message and exit"
}

### set euclidean options to default setting here
##

while getopts i:s:p:H:t:h option
do 
    case "${option}" in
        i)input=${OPTARG};;
        p)project=${OPTARG};;
        ### optional:
        s)step=${OPTARG};;
        H)HMM=${OPTARG};;
        t)threads=${OPTARG};;
        h)Help; exit;;
    esac

done

####################### SET AND CHECK VARIABLES/ARGUMENTS #####################################
OUTDIR="${PWD}/BINNING_${project}"
mkdir -p ${bioinftmp} ## make tmp data storage in scratch

######### from commandline arugments ######
if [[ -z "${input}" ]]; then echo -e "${red}-i, --input missing"; Help; exit; fi
if [[ -z "${project}" ]]; then echo -e "${red}-p, --project missing"; Help; exit; fi

### set defaults for optional arguments

### THREADS

if [[ ! -z "${SLURM_CPUS_PER_TASK}" ]]; then THREADS="${SLURM_CPUS_PER_TASK}"; echo -e "${green}Using threads set by SLURM_CPUS_PER_TASK:${nocolor}"; fi

if [[ ! -z "${threads}" ]]; then THREADS="${threads}"; echo -e "${green}Using threads set by --threads flag:${nocolor}"; fi

if [[ -z "${threads}" ]] && [[ -z "${SLURM_CPUS_PER_TASK}" ]]; then THREADS="1"; echo -e "${green}No SLURM_CPUS_PER_TASK or --threads set, using default threads:${nocolor}"; fi
echo -e "\t${THREADS} threads"

#######
if [[ -z $(conda env list | grep "bioinftools") ]]; then echo "NO conda environment for bioinftools found - see README"; else echo "conda environment for bioinftools FOUND"; fi

####################### SOFTWARE #####################################
### LOAD software available via shared environment on server:
module purge
eval "$(conda shell.bash hook)"
conda activate bioinftools

####################### DATABASES #####################################

########################################################################
########################################################################

####################### THE SCRIPT #####################################

################################################################################################
################################################################################################
########################################################################
########################################################################

###### STEP - build the hmm profiles
if [[ "$step" =~ "build" ]]; then
echo -e "${cyan}\t\tRUNNING build${nocolor}"

rm -f ${project}.hmm; touch ${project}.hmm

for f in ${input}/*; do
    N=$(echo $f | cut -d "." -f 1)
    clustalo -i $f -o ${N}.sto -v --outfmt=st
    hmmbuild ${N}.hmm ${N}.sto
    cat ${N}.hmm >> ${project}.hmm
    rm -f ${N}.hmm
done

##############################################################
fi
##############################################################

###### STEP - pre-process your contig fasta(s)
if [[ -z "${step}" ]] || [[ "$step" =~ "translate" ]]; then
echo -e "${cyan}\t\tRUNNING translate${nocolor}"

seqkit fq2fa -j $THREADS ${input} | esl-translate - > ${project}-trans.faa

echo -e "${cyan}\t\tRUNNING search${nocolor}"

hmmsearch --cpu $THREADS --tblout ${project}.pfamout ${HMM} ${project}-trans.faa

awk '/\tViruses;/' 091523_BIP23plus___barcode_61.bctrimmedreads.fastq.allcompare.kaiju | cut -f 12 | sort -V | uniq -c

##############################################################
fi
##############################################################