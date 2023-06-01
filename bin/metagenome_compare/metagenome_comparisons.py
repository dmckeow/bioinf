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
                    help="Full path to a file that lists input file locations, semi-colon delimited, 3 fields:\n \tFIELD1=contig fasta, separate line per assembly (full path). Can be individual assemblies or co-assemblies\n \tFIELD2=reads fq.gz, corresponding to each assembly contig fasta (full path)\n \tFIELD3=profile name; a short name for each profile which will appear in all outputs and figures\n\n\nprofile names MUST NOT begin with a digit, and MUST ONLY include letters, digits, and underscores and MUST be UNIQUE within your samples being run!!!\n!Example:\n\n/home/data/19AUG22_1_DM_barcode01.contigs.fasta;/home/data/19AUG22_1_DM_barcode01-trimmedbcs.fastq.gz;location1\n/home/data/19AUG22_1_DM_barcode03.contigs.fasta;/home/data/19AUG22_1_DM_barcode03-trimmedbcs.fastq.gz;location2\n/home/data/19AUG22_1_DM_barcode06.contigs.fasta;/home/data/19AUG22_1_DM_barcode06-trimmedbcs.fastq.gz;location3\n\n")
parser.add_argument("-s", "--step", 
                    help="Which script steps to run. A0 = run the whole script.\n\n\tMultiple steps can be specified, e.g.: A1_A2_A3 will only run steps A1-A3.\n\nIf running WITHOUT mapping - i.e. you only have contigs, but no reads, then add nomap to -s , e.g. A0_nomap\n\nALSO, some steps can be resumed without deleting previous progress for that step (see below), by adding the flag resume to -s but if you resume I recommend that you delete the last output file of that step, in case it is an incomplete output file\n\n IF using nomap then submit one job per sample as they cannot be merged and set -prefix to be the same as the profile name\n\n!!! STEPS AVAILABLE: !!! \n\n (Default A0)", 
                    default="A0")
parser.add_argument("-m", "--mincontigsize", 
                    help="Integer to set minimum size of contigs in bp to include in analyses, don't go below 1000! (Default: 2000)", 
                    default=2000)
parser.add_argument("-o", "--output", required=True, 
                    help="Absolute path for final output directory. Will be created if it does not already exist. An output folder named after ANVIO_(prefix) will be created here\n")
parser.add_argument("-p", "--prefix", required=True, 
                    help="A meaningful name for output folder. Any pre-existing driectory at this location and identical name will be deleted. Will also be the name of the merged analyses files. It should be a UNIQUE name that cannot possibly be identical to any other folder within your output destination e.g. -p 11OCT22project will make the folder ANVIO_11OCT22project\n\n")
parser.add_argument("-S", "--splitlength", 
                    help="\tInteger to set minimum size of split in bp. For RNA honeybee viruses we have used:-S 5000 (default: 20000)", 
                    default=20000)
parser.add_argument("-t", "--threads", 
                    help="\tNumber of threads (default: 24)", default=24)
parser.add_argument("--check-databases", action="store_true",
                    help="Check if all required databases are present and exit")


# parse the arguments
args = parser.parse_args()

input = args.input
output = args.output
step = args.step
mincontigsize = args.mincontigsize
prefix = args.prefix
splitlength = args.splitlength
threads = args.threads

####################### SOFTWARE, DATABSES, ETC ###########################

##### Seqkit - must be installed locally, e.g.:
SEQKIT="/panfs/jay/groups/27/dcschroe/shared/tools/seqkit"

##### Diamond - installed locally
DIAMOND="/panfs/jay/groups/27/dcschroe/shared/tools/diamond"

####################### THE SCRIPT #####################################
############## STEP A0 - runs STEPS A1-A9 & B1 ##############
#############################################################

############## STEP A1 - pre-process your contig fasta(s) ##############
os.system(f"{DIAMOND} blastx --max-target-seqs {MAXTARGETS} --evalue 1e-180 --sensitive -p {threads} -d {CUSTOM_DMND}/{f} -q {prefix}-SPLITS.fa -f 6 -o {f}.evalue_min.dmnd.blastx")

if "A0" in step or "A1" in step:
    print("\n\n##########################################", file=sys.stderr)
    print("Step A1 - sequences must be trimmed and then filtered down to a uniform size (100 bp e.g. from metavir)", file=sys.stderr)
    print("##########################################\n\n", file=sys.stderr)

### split ALL sequences to a minimum length, removing any below that length

### Each sample is compared to every other sample using tBLASTx (contigs vs contigs or reads vs reads)

### the best hit between each unique pair of sequences is retained and their total bitscores are summed per sample pair comparison

### instead of subsampling before tblastx, we could instead divide the total bitscore by the number of bitscores per sample pair comparison

### Finally, the resulting score matrix (i.e. similarity scores for all virome pairs) is used to cluster viromes using R software and the pvclust package, working with default parameters and 100 bootstraps

os.system(f"{DIAMOND} blastx --max-target-seqs {MAXTARGETS} --evalue 1e-180 --sensitive -p {threads} -d {CUSTOM_DMND}/{f} -q {prefix}-SPLITS.fa -f 6 -o {f}.evalue_min.dmnd.blastx")

    # make output parent directories if they do not already exist
    os.makedirs(output, exist_ok=True)
    
    os.system(f"dos2unix {input}")
    os.system(f"sed -i 's/\t/;/g' {input}")
    
    # make the final output directory
    final_output_dir = os.path.join(output, f"ANVIO_{prefix}")
    os.makedirs(final_output_dir, exist_ok=True)
    os.chdir(final_output_dir)
    
    os.makedirs(f"{prefix}-reformat", exist_ok=True)
    os.chdir(f"{prefix}-reformat")
    
    ##### fix fasta deflines and save "deflinekey" files to show corresponding original format and new format deflines
    with open(input, "r") as f:
        for sets in f:
            contigs = sets.split(";")[0].strip() # remove any leading/trailing whitespace
            profilename = sets.split(";")[2].strip()
        # run the anvio reformat script now
            os.system(f"anvi-script-reformat-fasta -o {profilename}.fa --simplify-names --report-file {profilename}.deflinekey -l {mincontigsize} --seq-type NT --prefix 'contig_{profilename}' {contigs}")
    
    ##### concatenate fasta files for the project (but ignores the file if it doesnt exist)
    os.chdir(final_output_dir)
    if os.path.exists(f"{prefix}.fa"):
        os.remove(f"{prefix}.fa")
        open(f"{prefix}.fa", "w").close()
    if os.path.exists(f"{prefix}.deflinekey"):
        os.remove(f"{prefix}.deflinekey")
        open(f"{prefix}.deflinekey", "w").close()

    with open(f"{prefix}.fa", "a") as outfile:
        for file in sorted(os.listdir(f"{prefix}-reformat")):
            if file.endswith(".fa"):
                with open(os.path.join(f"{prefix}-reformat", file), "r") as infile:
                    outfile.write(infile.read())

    with open(f"{prefix}.deflinekey", "a") as outfile:
        for file in sorted(os.listdir(f"{prefix}-reformat")):
            if file.endswith(".deflinekey"):
                with open(os.path.join(f"{prefix}-reformat", file), "r") as infile:
                    outfile.write(infile.read())

os.system(f"rm -fr {prefix}-reformat")




print("\n\n#############################################################################", file=sys.stderr)
print("Workflow finished, check the log file that each step was complete as expected", file=sys.stderr)
print("#############################################################################\n\n", file=sys.stderr)
