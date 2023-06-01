#!/usr/bin/env python3
# -*- coding: utf-8

####################### PYTHON IMPORTS ###########################

import argparse
import os
import subprocess
import shlex
import sys
import pysam

####################### VARIABLES TO SET BEFORE PARSING THE ARGUMENTS ###########################

###os.path.basename('guyana_primate_barcode14.contigs.fasta').replace(".", "_")

# define a function to compute the new file name
##def get_file_name_with_underscores(filename):
    ##basename = os.path.basename(filename)
    ##name_without_path, ext = os.path.splitext(basename)
    ##name_with_underscores = name_without_path.replace(".", "_")
    ##return name_with_underscores + ext

####################### PARSE THE ARGUMENTS AND INPUT FLAGS ###########################

# create the argparse parser
parser = argparse.ArgumentParser()

# add the arguments
parser.add_argument("-i", "--input", required=True, 
                    help="Full path to a file that lists input file locations, semi-colon delimited, 3 fields:\n \tFIELD1=contig fasta, separate line per assembly (full path). Can be individual assemblies or co-assemblies\n \tFIELD2=reads fq.gz, corresponding to each assembly contig fasta (full path)\n \tFIELD3=profile name; a short name for each profile which will appear in all outputs and figures\n\n\nprofile names MUST NOT begin with a digit, and MUST ONLY include letters, digits, and underscores and MUST be UNIQUE within your samples being run!!!\n!Example:\n\n/home/data/19AUG22_1_DM_barcode01.contigs.fasta;/home/data/19AUG22_1_DM_barcode01-trimmedbcs.fastq.gz;location1\n/home/data/19AUG22_1_DM_barcode03.contigs.fasta;/home/data/19AUG22_1_DM_barcode03-trimmedbcs.fastq.gz;location2\n/home/data/19AUG22_1_DM_barcode06.contigs.fasta;/home/data/19AUG22_1_DM_barcode06-trimmedbcs.fastq.gz;location3\n\n")
parser.add_argument("-s", "--step", 
                    help="Which script steps to run. A0 = run the whole script.\n\n\tMultiple steps can be specified, e.g.: A1_A2_A3 will only run steps A1-A3.\n\nIf running WITHOUT mapping - i.e. you only have contigs, but no reads, then add nomap to -s , e.g. A0_nomap\n\nALSO, some steps can be resumed without deleting previous progress for that step (see below), by adding the flag resume to -s but if you resume I recommend that you delete the last output file of that step, in case it is an incomplete output file\n\n IF using nomap then submit one job per sample as they cannot be merged and set -prefix to be the same as the profile name\n\n!!! STEPS AVAILABLE: !!! \n\n (Default A0)", 
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


####################### THE SCRIPT #####################################
############## STEP A0 - runs STEPS A1-A9 & B1 ##############
#############################################################

############## STEP A1 - pre-process your contig fasta(s) ##############
os.system(f"{DIAMOND} blastx --max-target-seqs {MAXTARGETS} --evalue 1e-180 --sensitive -p {threads} -d {CUSTOM_DMND}/{f} -q {prefix}-SPLITS.fa -f 6 -o {f}.evalue_min.dmnd.blastx")

if "A0" in step or "A1" in step:
    print("\n\n##########################################", file=sys.stderr)
    print("Step A1 - sequences must be trimmed and then filtered down to a uniform size (100 bp e.g. from metavir)", file=sys.stderr)
    print("##########################################\n\n", file=sys.stderr)

### to normalize the number of sequences across all samples, subsample to a certain number of sequences (50,000)

### Each sample is compared to every other sample using tBLASTx (contigs vs contigs or reads vs reads)

### the best hit between each unique pair of sequences is retained and their total bitscores are summed per sample pair comparison

### instead of subsampling before tblastx, we could instead divide the total bitscore by the number of bitscores per sample pair comparison

### Finally, the resulting score matrix (i.e. similarity scores for all virome pairs) is used to cluster viromes using R software and the pvclust package, working with default parameters and 100 bootstraps






print("\n\n#############################################################################", file=sys.stderr)
print("Workflow finished, check the log file that each step was complete as expected", file=sys.stderr)
print("#############################################################################\n\n", file=sys.stderr)
