#!/usr/bin/env python3

import os
import argparse

parser = argparse.ArgumentParser(description="Perform multiple sequence alignment and phylogenetic tree inference from FASTA sequences using mafft and IQtree. Autodetects input sequence type and phylogenetic model to use", epilog="""USAGE EXAMPLE: sbatch --time=8:00:00 --cpus-per-task=8 --mem=24GB --partition partition_name -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; bioinf-phylogeny.py -i file.fa -p bee_virus_project'""")
parser.add_argument("-i", "--input", help="Path to a single input fasta file(s) OR a file containing a list of paths to multiple fasta files", required=True)
parser.add_argument("-p", "--project", help="A name for your project - all output files will include this name", required=True)
parser.add_argument("-b", "--bootstraps", help="Specifies the number of bootstrap replicates (default 1000).", type=int, default=1000)
parser.add_argument("-N", "--nophylogeny", action='store_true', help="just run alignment (default false)")
parser.add_argument("-P", "--phylogenyonly", action='store_true', help="only run the phylogeny iqtree (default false)")
parser.add_argument("-H", "--hyphy", action='store_true', help="Run hyphy GARD recombination analyses (default false)")
parser.add_argument("-t", "--threads", help="number of threads to use for mafft (default=4)", type=int, default=4)

args = parser.parse_args()

## specific settings
input_arg = args.input
project = args.project
bootstraps = args.bootstraps

## step settings
nophylogeny = args.nophylogeny
phylogenyonly = args.phylogenyonly
hyphy = args.hyphy

## generic settings
THREADS = args.threads ## use --threads whether provided or not (uses default set in argparse)
SLURM_THREADS = os.getenv("SLURM_CPUS_PER_TASK") ## if slurm --cpus-per-task provided, overwrite --threads with it instead
if SLURM_THREADS:
	THREADS = SLURM_THREADS
	print(f'Threads set by SLURM --cpus-per-task: {THREADS}')

## set in-script variables
concat_fasta = project + ".concat" + ".fa"
aln_file = project + ".aln"

# Concatenate reference and input files
def concatenate_files(input_files, output):
	with open(output, 'w') as outfile:
		for file in input_files:
			with open(file, 'r') as infile:
				outfile.write(infile.read())

def parse_input(input_arg):
	# Check if input_arg is a file
	if os.path.isfile(input_arg):
		with open(input_arg, 'r') as file:
			first_line = file.readline()
			if not first_line.startswith('>'): # check if it is a fasta file
				file.seek(0) # reset back to start of file
				fasta_files = [line.strip() for line in file if line.strip()] # Read paths to multiple fasta files
			else:
				# If input_arg is a single fasta file
				fasta_files = [input_arg]
	return fasta_files

# Parse input argument
input_files = parse_input(input_arg)

concatenate_files(input_files, concat_fasta)

# Perform multiple sequence alignment
if not phylogenyonly:
	os.system(f"seqkit rmdup {concat_fasta} > {concat_fasta}.tmp && mv {concat_fasta}.tmp {concat_fasta}")
	os.system(f"mafft --thread {THREADS} --adjustdirectionaccurately --auto --reorder --maxiterate 1000 {concat_fasta} > {aln_file}")

# Do phylogeny
if not nophylogeny:
	os.system(f"iqtree -s {aln_file} --prefix {project} -m MFP -B {bootstraps} -T {THREADS}")

## RUN hyphy gard if flag provided
if hyphy:
	os.system(f"hyphy CPU={THREADS} gard --alignment {aln_file}")