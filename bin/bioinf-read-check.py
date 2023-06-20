#!/usr/bin/env python3

import os
import argparse
import gzip
import re
import glob
import gzip
import subprocess
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Compare one sequence file to any number of other sequence files, providing summaries of the presence of every sequence of that file in every other file provided. Whilst it can be used for any type of sequence, an example of its use is to compare raw reads by sequence name to show the proportions of reads that made it through to various steps in an assembly - i.e. barcode trimming, correction, and finally contigs", epilog="""USAGE EXAMPLE: sbatch --time=24:00:00 --cpus-per-task=12 --mem=240GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; bioinf-read-check.py -c /panfs/jay/groups/27/dcschroe/dmckeow/data/CANU_molecular_survey/CANU_molecular_survey_barcode11 -p testreadcheck'""")
parser.add_argument("-i", "--input", help="Paths to input fastq and/or fasta file(s). The first (leftmost) file will be comapred to all the other files, to assess the presence or absence of each sequence in the first file in each of the other files. You must provide at least 2 files and can provide as many as you want, separated by spaces", nargs='+')
parser.add_argument("-d", "--differ", help="If any of your files provided in --input have headers that do not match the first file - for example if your reads were assigned different names in a processed read file, or if you want to summarise reads and contigs, then you can provide a tab-separated file with 2 columns: column 1 is the name in your first sequence file, and column 2 is an equivalent name found in any other sequence file")
parser.add_argument("-c", "--canu", help="Path to directory containing canu output file(s)for a single assembly. With this option, the script will automatically set --input and --differ using the output files from canu")
parser.add_argument("-t", "--threads", help="number of threads (default=4)", type=int, default=4)

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
	differ_file = glob.glob(differ_pattern)[0]

## step settings

## generic settings
THREADS = args.threads ## use --threads whether provided or not (uses default set in argparse)
SLURM_THREADS = os.getenv("SLURM_CPUS_PER_TASK") ## if slurm --cpus-per-task provided, overwrite --threads with it instead
if SLURM_THREADS:
	THREADS = SLURM_THREADS
	print(f'Threads set by SLURM --cpus-per-task: {THREADS}')
bioinfdb = os.environ["bioinfdb"]
bioinftmp = os.environ["bioinftmp"]

	
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
		if not re.match(r'^[0-9]+___.*\.fastaheaders$', file_name):
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

def compare_and_append(file1_path, file_set_pattern):
	# Read File1
	with open(file1_path, 'r') as file1:
		file1_lines = file1.readlines()

	# Find files in the set using the provided pattern
	file_set = []
	for file_name in os.listdir('.'):
		if re.match(file_set_pattern, file_name) and os.path.isfile(file_name):
			file_set.append(file_name)

	# Compare File1 to each file in the set
	for file_name in file_set:
		appended_lines = []
		with open(file_name, 'r') as file:
			file_lines = file.readlines()

		# Iterate over lines in File1
		for line in file1_lines:
			column1, column2 = line.strip().split('\t')

			# Check if column2 of File1 matches column1 in the current file
			if any(column2 == line.strip().split('\t')[0] for line in file_lines):
				# Append column1 from the matching line in File1 to the current file
				appended_lines.append(column1 + '\n')

		appended_content = ''.join(appended_lines)

		with open(file_name, 'a') as file:
			file.write(appended_content)

## prepare the sequence headers
process_input_files(input_files)

## delete tmp files
#for f in glob.glob("*.tmp"):
	#os.remove(f)


## differ fix
if differ_file:
	compare_and_append(differ_file, r'^[0-9]+___.*\.fastaheaders$')


## compare sequence headers to count read naem occurances in the headers of each sequence file
search_file_pattern = r'^1___.*\.fastaheaders$' # Pattern to match the search file names
compare_headers(search_file_pattern)

## delete the comparison between sequence file 1 and itself
for f in glob.glob("1___*.fastaheaders.compare.1.tmp"):
	os.remove(f)

#### run kaiju on the first sequence file
for f in os.listdir(os.path.join(bioinfdb, "KAIJU")):
	f_path = os.path.join(bioinfdb, "KAIJU", f)
	if os.path.isdir(f_path):
		kaiju_input = input_files[0]
		file_name = os.path.basename(kaiju_input)
		kaiju_output = os.path.join(os.path.splitext(file_name)[0] + "_" + os.path.basename(f_path))

		nodes_dmp = os.path.join(f_path, "nodes.dmp")
		fmi_files = os.path.join(f_path, "*.fmi")
		names_dmp = os.path.join(f_path, "names.dmp")

		# Execute kaiju command
		kaiju_command = f"kaiju -t {nodes_dmp} -f {fmi_files} -i {kaiju_input} -o {kaiju_output}.kaiju.1.tmp -z {THREADS} -v"
		#subprocess.run(kaiju_command, shell=True)

		# Execute sort command
		sort_command = f"sort -t $'\t' -V -k 2,2 {kaiju_output}.kaiju.1.tmp -o {kaiju_output}.kaiju.2.tmp"
		#subprocess.run(sort_command, shell=True)

		# Execute kaiju-addTaxonNames command
		add_taxon_names_command = f"kaiju-addTaxonNames -t {nodes_dmp} -n {names_dmp} -i {kaiju_output}.kaiju.2.tmp -o {kaiju_output}.kaiju.names.tmp -r superkingdom,phylum,order,class,family,genus,species"
		#subprocess.run(add_taxon_names_command, shell=True)


# Reduce to the best kaiju hit per query across all databases queried
kaiju_names_pattern = "*.kaiju.names.tmp"

kaiju_input = input_files[0]
file_name = os.path.basename(kaiju_input)
kaiju_output = os.path.splitext(file_name)[0]

reduce_to_best_kaiju = f"sed 's/$/\t0\tNA\tNA\tNA\tNA/g' {kaiju_names_pattern} | cut -f 1-8 | sort -t $'\t' -k 2,2V -k 4,4nr | awk -F '\t' '!a[$2]++' | sed 's/0\tNA\tNA\tNA\tNA//g' > {kaiju_output}.kaiju.merge"
subprocess.run(reduce_to_best_kaiju, shell=True)


def sort_by_column(file_path, column_index_for_sort):
	# Read the file into a DataFrame
	df = pd.read_csv(file_path, sep='\t')

	# Perform natural sort by the specified column
	df_sorted = df.sort_values(by=df.columns[column_index_for_sort])

	# Overwrite the input file with the sorted version
	df_sorted.to_csv(file_path, sep='\t', index=False)

	print("File sorted and overwritten:", file_path)


def concatenate_files(patterns, output_file, column_index_for_sort):
	file_list = []  # List to store file paths

	for pattern in patterns:
		file_list.extend(glob.glob(pattern))  # Get a list of files matching each pattern

	if not file_list:
		print("No files found matching the patterns.")
		return

	dfs = []  # List to store DataFrames

	# Read each file and append it to the list
	for file_name in file_list:
		df = pd.read_csv(file_name, sep='\t')
		dfs.append(df)
		print("Concatenating the following files:", file_name)

	# Concatenate the DataFrames horizontally
	merged_df = pd.concat(dfs, axis=1)

	# Add file names as column headers
	merged_df.columns = file_list

	# sort by column specified before adding headers
	merged_df_sorted = merged_df.sort_values(by=merged_df.columns[column_index_for_sort])

	# Save the merged DataFrame as a tab-separated file
	merged_df_sorted.to_csv(output_file, sep='\t', index=False)

	print("Merged DataFrame saved to", output_file)


# sort kaiju merge by column 2 (read names)
kaiju_input = input_files[0]
file_name = os.path.basename(kaiju_input)
kaiju_merge = os.path.splitext(file_name)[0] + ".kaiju.merge"
sort_by_column(kaiju_merge, 1)

### horizontal concatenation of the comparison files and the fasta headers of the first file, and kaiju merge
file_name = os.path.basename(input_files[0])
output_file = os.path.splitext(file_name)[0] + ".concatcompare.tmp"

patterns = ['1___*fastaheaders', '*___*fastaheaders.compare.1.tmp', '*.kaiju.merge']
concatenate_files(patterns, output_file, 0)



