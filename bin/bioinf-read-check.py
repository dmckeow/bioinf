#!/usr/bin/env python3

import os
import argparse
import gzip
import re
import glob
import gzip
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Perform multiple sequence alignment and phylogenetic tree inference from FASTA sequences using mafft and IQtree. Autodetects input sequence type and phylogenetic model to use", epilog="""USAGE EXAMPLE: sbatch --time=8:00:00 --cpus-per-task=12 --mem=36GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; bioinf-read-check.py -c /panfs/jay/groups/27/dcschroe/dmckeow/data/CANU_molecular_survey/CANU_molecular_survey_barcode11 -p testreadcheck'""")
parser.add_argument("-i", "--input", help="Path to input fastq and/or fasta file(s). You can provide multiple paths to files, separated by spaces", nargs='+')
parser.add_argument("-d", "--differ", help="If any of your files provided in --input have headers that do not match the first file - for example if your reads were assigned different names in a processed read file, or if you want to summarise reads and contigs, then you can provide a tab-separated file with 2 columns: column 1 is the name in your first sequence file, and column 2 is an equivalent name found in any other sequence file")
parser.add_argument("-c", "--canu", help="Path to directory containing canu output file(s)for a single assembly")
parser.add_argument("-p", "--project", help="A name for your project - all output files will include this name", required=True)
parser.add_argument("-t", "--threads", help="number of threads to use for mafft (default=4)", type=int, default=4)

args = parser.parse_args()

## check for variables
if (args.input is None and args.canu is None) or (args.input is not None and args.canu is not None):
	parser.error("Provide either --input or --canu, and not both")

## specific settings
## if input (non-canu) provided, then it just uses whatever input and differ is given
if args.input:
	input_files = args.input
	differ_file = args.differ

## if canu provided, then it globs specific expected files wthin the directory given
if args.canu:
	patterns = ['*.bctrimmedreads.fastq*', '*.correctedReads.fasta*', '*.trimmedReads.fasta*', '*.contigs.fasta*'] # list specific canu files in order
	files = []
	for pattern in patterns:
		file_pattern = args.canu + '/' + pattern
		files.extend(glob.glob(file_pattern))
	# Set sorted canu files as input
	input_files = sorted(files, key=lambda x: [patterns.index(p) for p in patterns if p in x]) ## sort but keep the order specified earlier

	## the differ file is set to the read to tig file output of canu
	differ_pattern = args.canu + '/*.contigs.layout.readToTig.assembledonly.origNames'
	differ_file = glob.glob(differ_pattern)

project = args.project

## step settings

## generic settings
THREADS = args.threads ## use --threads whether provided or not (uses default set in argparse)
SLURM_THREADS = os.getenv("SLURM_CPUS_PER_TASK") ## if slurm --cpus-per-task provided, overwrite --threads with it instead
if SLURM_THREADS:
	THREADS = SLURM_THREADS
	print(f'Threads set by SLURM --cpus-per-task: {THREADS}')


	
def process_input_files(input_files):
	counter = 1 # set counter variable for file output names
	for input_file in input_files:
		output_file = get_output_filename(input_file, counter)
		# Check if the file is gzipped
		is_gzipped = False
		if input_file.endswith('.gz'):
			is_gzipped = True
		# Process the file based on its format
		if input_file.endswith('.fasta') or input_file.endswith('.fa'):
			extract_fasta_headers(input_file, output_file + '.fastaheaders')
		elif input_file.endswith('.fastq') or input_file.endswith('.fq'):
			convert_fastq_to_fasta(input_file, output_file + '.fasta.tmp')
			extract_fasta_headers(output_file  + '.fasta.tmp', output_file + '.fastaheaders')
		elif is_gzipped and (input_file.endswith('.fastq.gz') or input_file.endswith('.fq.gz')):
			uncompress_gzip(input_file, output_file + '.fastq.tmp')
			convert_fastq_to_fasta(output_file + '.fastq.tmp', output_file + '.fasta.tmp')
			extract_fasta_headers(output_file  + '.fasta.tmp', output_file + '.fastaheaders')
		elif is_gzipped and (input_file.endswith('.fasta.gz') or input_file.endswith('.fa.gz')):
			uncompress_gzip(input_file, output_file + '.fasta.tmp')
			extract_fasta_headers(output_file + '.fasta.tmp', output_file + '.fastaheaders')
		else:
			print(f"Unsupported file format: {input_file}")
		counter += 1 # the counter increases by one after each file

def extract_fasta_headers(input_file, output_file):
	with open(input_file, 'r') as fasta_file, open(output_file, 'w') as headers_file:
		for record in SeqIO.parse(fasta_file, 'fasta'):
			headers_file.write(f'{record.id}\n')

def convert_fastq_to_fasta(input_file, output_file):
	with open(output_file, 'w') as fasta_file:
		for record in SeqIO.parse(input_file, 'fastq'):
			SeqIO.write(record, fasta_file, 'fasta')

def uncompress_gzip(input_file, output_file):
	with gzip.open(input_file, 'rt') as gzipped_file, open(output_file, 'w') as output:
		output.writelines(gzipped_file)

def get_output_filename(input_file, counter):
	file_name = os.path.basename(input_file)
	file_name_without_ext = os.path.splitext(file_name)[0]
	return f"{counter}___{file_name_without_ext}"

##################
def compare_headers(search_file_pattern):
	current_directory = os.getcwd()  # Get the current working directory

	# Find the search file 1___ for fasta headers
	search_file = None
	for file_name in os.listdir(current_directory):
		if re.match(search_file_pattern, file_name):
			search_file = os.path.join(current_directory, file_name)
			break

	if not search_file:
		print(f"Search file matching pattern '{search_file_pattern}' not found.")
		return

	# Read the list of fastaheaders from the search file
	with open(search_file, 'r') as f:
		search_words = [line.strip() for line in f]

	# Iterate over each fastaheader list file in the directory
	for file_name in os.listdir(current_directory):
		if not re.match(r'^[0-9]___.*\.fastaheaders$', file_name):
			continue  # Skip files that do not match the pattern
		file_path = os.path.join(current_directory, file_name)
		output_file = f"{file_name}.compare.1.tmp"  # Output file named after the searched file

		# Open the file and read its contents
		with open(file_path, 'r') as f:
			found_words = []  # List to store found words in the file
			for line in f:
				# Check if any word from the search file is found in the line
				if any(search_word in line for search_word in search_words):
					# Add the found word to the list
					found_words.extend([search_word for search_word in search_words if search_word in line])

			# Create a set of found words for faster comparison
			found_words_set = set(found_words)

		# Write "1" if a header is found, otherwise write "0" 
		with open(output_file, 'w') as f:
			for search_word in search_words:
				if search_word in found_words_set:
					f.write("1\n")
				else:
					f.write("0\n")

def search_first_vs_differ(differ_file, first_sequence_file):
	# Read lines from File 2 and store them in a set
	with open(first_sequence_file, 'r') as file2:
		lines_file2 = {line.strip() for line in file2}

	file2_name = os.path.splitext(os.path.basename(first_sequence_file))[0]
	output_path = f"{file2_name}_search_first_vs_differ.tmp"

	# Search for matching lines in File 1 and write them to the output file
	with open(differ_file, 'r') as file1, open(output_path, 'w') as output:
		for line in file1:
			read, contig = line.strip().split('\t')
			if read in lines_file2:
				output.write(line)

def search_and_append_matches(file1_path, file_regex_pattern):
	# Read column 2 of File 1 and store it in a set
	with open(file1_path, 'r') as file1:
		column2_file1 = {line.strip().split('\t')[1] for line in file1}

	# Compile the regex pattern
	file_regex = re.compile(file_regex_pattern)

	# Search for matching lines in the set of files and append column 2 whenever a match occurs
	files_to_search = [filename for filename in os.listdir('.') if re.match(file_regex, filename)]
	for filename in files_to_search:
		with open(filename, 'a') as file:
			for line in column2_file1:
				if line == filename:
					continue  # Skip appending if the line matches the filename
				file.write('\t' + line + '\n')

## prepare the sequence headers
process_input_files(input_files)

## delete tmp files
#for f in glob.glob("*.tmp"):
#	os.remove(f)

if differ_file:
	search_file_pattern = r'^1___.*\.fastaheaders$'
	current_directory = os.getcwd()  # Get the current working directory
	# Find the search file 1___ for fasta headers
	first_sequence_file = None
	for file_name in os.listdir(current_directory):
		if re.match(search_file_pattern, file_name):
			first_sequence_file = os.path.join(current_directory, file_name)
			break
	for differ in differ_file:
		differ = differ
	search_first_vs_differ(differ, first_sequence_file)

# append alternate read names into appropriate fastaheader files before counting read presence
	file2_name = os.path.splitext(os.path.basename(first_sequence_file))[0]
	search_first_vs_differ_file = f"{file2_name}_search_first_vs_differ.tmp"
	file_regex_pattern = r'^[0-9]___.*\.fastaheaders$'
	search_and_append_matches(search_first_vs_differ_file, file_regex_pattern)

## compare sequence headers
search_file_pattern = r'^1___.*\.fastaheaders$' # Pattern to match the search file names
compare_headers(search_file_pattern)

# Perform multiple sequence alignment
#os.system(f"mafft --thread {THREADS} --adjustdirectionaccurately --auto --reorder --maxiterate 1000 {tmp_file} > {aln_file}")

# Remove temporary files
#os.remove(tmp_file)

