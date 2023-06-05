#!/usr/bin/env python3

import os
import argparse

parser = argparse.ArgumentParser(description="Perform multiple sequence alignment and phylogenetic tree inference from FASTA sequences using mafft and IQtree. Autodetects input sequence type and phylogenetic model to use", epilog="""USAGE EXAMPLE: sbatch --time=8:00:00 --cpus-per-task=8 --mem=24GB --partition partition_name -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; bioinf-phylogeny.py -i file.fa -p bee_virus_project'""")
parser.add_argument("-i", "--input", help="Path to input fasta file(s). You can provide multiple paths to files, separated by spaces", required=True, nargs='+')
parser.add_argument("-p", "--project", help="A name for your project - all output files will include this name", required=True)
parser.add_argument("-b", "--bootstraps", help="Specifies the number of bootstrap replicates (default 1000).", type=int, default=1000)
parser.add_argument("-N", "--nophylogeny", action='store_true', help="skip running iqtree phylogeny (default false)")
parser.add_argument("-H", "--hyphy", action='store_true', help="Run hyphy GARD recombination analyses (default false)")
parser.add_argument("-t", "--threads", help="number of threads to use for mafft (default=4)", type=int, default=4)

args = parser.parse_args()

## specific settings
input_files = args.input
project = args.project
bootstraps = args.bootstraps

## step settings
nophylogeny = args.nophylogeny
hyphy = args.hyphy

## generic settings
THREADS = args.threads ## use --threads whether provided or not (uses default set in argparse)
SLURM_THREADS = os.getenv("SLURM_CPUS_PER_TASK") ## if slurm --cpus-per-task provided, overwrite --threads with it instead
if SLURM_THREADS:
	THREADS = SLURM_THREADS
	print(f'Threads set by SLURM --cpus-per-task: {THREADS}')

## set in-script variables
tmp_file = "tmp." + project + ".fa"
aln_file = project + ".aln"


# Concatenate reference and input files
def concatenate_files(input, output):
    with open(output, 'w') as outfile:
        for file in input:
            with open(file, 'r') as infile:
                outfile.write(infile.read())

concatenate_files(input_files, tmp_file)

# Perform multiple sequence alignment
os.system(f"mafft --thread {THREADS} --adjustdirectionaccurately --auto --reorder --maxiterate 1000 {tmp_file} > {aln_file}")

# Remove temporary files
os.remove(tmp_file)

# Do phylogeny
if not nophylogeny:
	os.system(f"iqtree -s {aln_file} -redo --prefix {project} -m MFP -B {bootstraps} -T {THREADS}")

## RUN hyphy gard if flag provided
if hyphy:
	os.system(f"hyphy CPU={THREADS} gard --alignment {aln_file}")