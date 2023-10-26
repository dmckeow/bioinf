#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import gzip
import re
import glob
import gzip
import subprocess
import pandas as pd
from Bio import SeqIO
import shlex

parser = argparse.ArgumentParser(description="Bioinf-read-check uses a successful CANU assembly (i.e. with contigs) and creates a series of plots to show read classifications (by Kaiju) and a breakdown of how many reads for taxa failed or passed each step of the CANU assembly process", epilog="""USAGE EXAMPLE (DO NOT run for multiple assemblies from the same directory - it will overwrite any previous outputs from this script): sbatch --time=24:00:00 --cpus-per-task=12 --mem=240GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; bioinf-read-check.py -c /panfs/jay/groups/27/dcschroe/dmckeow/data/CANU_molecular_survey/CANU_molecular_survey_barcode11'""")
parser.add_argument("-i", "--input", help="NO LONGER IN USE", nargs='+')
parser.add_argument("-d", "--differ", help="NO LONGER IN USE")
parser.add_argument("-c", "--canu", help="Path to directory containing canu output file(s)for a single assembly. With this option, the script will automatically set --input and --differ using the output files from canu", required=True)
parser.add_argument("-t", "--threads", help="number of threads (default=4)", type=int, default=4)

args = parser.parse_args()

## delete any files from previous runs if needed
for f in glob.glob("*___*.fastaheaders*"):
	os.remove(f)
for f in glob.glob("*.kaiju.*.tmp"):
	os.remove(f)
for f in glob.glob("*.kaiju.merge"):
	os.remove(f)
for f in glob.glob("*allcompare.*"):
	os.remove(f)
for f in glob.glob("read_header"):
	os.remove(f)

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

################## original
def compare_headers_OLD(search_file_pattern):
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
		if not re.match(r'^[0-9]+___.*\.fastaheaders$', file_name) and not re.match(r'^1___.*\.fastaheaders$', file_name):
			continue  # Skip files that do not match the pattern
		file_path = os.path.join(current_directory, file_name)
		output_file = f"{file_name}.compare.1.tmp"  # Output file named after the searched file
		print(f"BEGIN Compare header on '{file_name}'")
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
		print(f"Compare header on '{file_name}' DONE")



####### TESTING!
def compare_headers(search_file_pattern):
	current_directory = os.getcwd()

	# Find the search file for fasta headers
	search_file = None
	for file_name in os.listdir(current_directory):
		if re.match(search_file_pattern, file_name):
			search_file = os.path.join(current_directory, file_name)
			break

	if not search_file:
		print(f"Search file matching pattern '{search_file_pattern}' not found.")
		return

	# Read the list of fasta headers from the search file
	with open(search_file, 'r') as f:
		search_words = {line.strip() for line in f}  # Use a set for faster membership tests

	# Compile the regex pattern outside the loop
	header_file_pattern = re.compile(r'^[0-9]+___.*\.fastaheaders$')

	# Iterate over each fasta header list file in the directory
	for file_name in os.listdir(current_directory):
		if not header_file_pattern.match(file_name):
			continue
		if file_name == os.path.basename(search_file):
			continue ## skip if it is the search file 1__

		file_path = os.path.join(current_directory, file_name)
		output_file = f"{file_name}.compare.1.tmp"
		print(f"BEGIN Compare header on '{file_name}'")
		# Read the entire file into memory
		with open(file_path, 'r') as f:
			file_content = f.read()

		found_words = {search_word for search_word in search_words if search_word in file_content}

		with open(output_file, 'w') as f:
			for search_word in search_words:
				f.write("1\n" if search_word in found_words else "0\n")
		print(f"Compare header on '{file_name}' DONE")

#######
def concatenate_files(pattern, output_file):
	with open(output_file, 'w') as outfile:
		file_list = glob.glob(pattern)
		for file_name in file_list:
			with open(file_name, 'r') as infile:
				outfile.write(infile.read())
## original
def compare_and_append_OLD(file1_path, file_set_pattern):
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
		print(f"BEGIN Compare and append on '{file_name}'")
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
		print(f"Compare and append on '{file_name}' DONE")

##### NEW - faster!
def compare_and_append(file1_path, file_set_pattern):
	# Read File1
	with open(file1_path, 'r') as file1:
		file1_lines = file1.readlines()

	# Find files in the set using the provided pattern
	file_set = []
	for file_name in os.listdir('.'):
		if re.match(file_set_pattern, file_name) and os.path.isfile(file_name):
			file_set.append(file_name)
	# Compile the regex pattern outside the loop
	file_set_pattern_compiled = re.compile(file_set_pattern)

	# Create a dictionary from file_lines for faster matching
	file_lines_dict = {}
	for file_name in file_set:
		print(f"BEGIN Compare and append on '{file_name}'")
		with open(file_name, 'r') as file:
			file_lines = file.readlines()
			file_lines_dict[file_name] = {line.strip().split('\t')[0] for line in file_lines}

	# Iterate over lines in File1
	appended_lines = []
	for line in file1_lines:
		column1, column2 = line.strip().split('\t')
		for file_name in file_set:
			if column2 in file_lines_dict[file_name]:
				# Append column1 from the matching line in File1 to the current file
				appended_lines.append(column1 + '\n')
				break  # Break once a match is found in a file

	appended_content = ''.join(appended_lines)

	# Append the content to each file in the file_set
	for file_name in file_set:
		with open(file_name, 'a') as file:
			file.write(appended_content)
		print(f"Compare and append on '{file_name}' DONE")

## prepare the sequence headers
process_input_files(input_files)

## delete tmp files
for f in glob.glob("*.tmp"):
	os.remove(f)


## differ fix
if differ_file:
	compare_and_append(differ_file, r'^[0-9]+___.*\.fastaheaders$')

## compare sequence headers to count read name occurances in the headers of each sequence file
search_file_pattern = r'^1___.*\.fastaheaders$' # Pattern to match the search file names
compare_headers(search_file_pattern)

#### run kaiju on the first sequence file
for f in os.listdir(os.path.join(bioinfdb, "KAIJU")):
	f_path = os.path.join(bioinfdb, "KAIJU", f)
	if os.path.isdir(f_path):
		kaiju_input = input_files[0]
		file_name = os.path.basename(kaiju_input)
		kaiju_output = os.path.join(os.path.splitext(file_name)[0] + "_" + os.path.basename(f_path))

		nodes_dmp = os.path.join(f_path, "nodes.dmp")
		kaiju_fmi = os.path.join(f_path, os.path.basename(f_path) + ".fmi")
		names_dmp = os.path.join(f_path, "names.dmp")

		# Execute kaiju command
		cmd = "kaiju -v -t "+ nodes_dmp +" -f "+ kaiju_fmi +" -i "+ kaiju_input +" -o "+ kaiju_output +".kaiju.1.tmp -z "+ THREADS
		kaiju_command = shlex.split(cmd)
		subprocess.call(kaiju_command)

		# Execute kaiju-addTaxonNames command
		add_taxon_names_command = f"kaiju-addTaxonNames -t {nodes_dmp} -n {names_dmp} -i {kaiju_output}.kaiju.1.tmp -o {kaiju_output}.kaiju.names.tmp -r superkingdom,phylum,order,class,family,genus,species"
		subprocess.run(add_taxon_names_command, shell=True)


# Reduce to the best kaiju hit per query across all databases queried
file_pattern = "*.kaiju.names.tmp"
for file_path in glob.glob(file_pattern):
	# Read the input file and store the modified rows in a list
	modified_rows = []
	with open(file_path, "r") as file:
		for line in file:
			row = line.rstrip("\n").split("\t")
			row.extend(['0', '-', '-', '-', '-'])
			modified_row = "\t".join(row[:8]) + "\n"
			modified_rows.append(modified_row)

	# Overwrite the file with the modified data
	with open(file_path, "w") as file:
		file.writelines(modified_rows)

	print("File overwritten:", file_path)

file_name = os.path.basename(input_files[0])
kaiju_merge = os.path.splitext(file_name)[0] + ".kaiju.merge"
concatenate_files("*.kaiju.names.tmp", kaiju_merge)

#### be CAREFUL - this way of sorting adds an extra column AND row
df = pd.read_csv(kaiju_merge, sep='\t', header=None)
col_to_convert = [3]
df[col_to_convert] = df[col_to_convert].apply(pd.to_numeric)
sorted_df = df.sort_values(by=[1, 3], ascending=[True, False])
sorted_df.drop_duplicates(subset=[1], keep='first', inplace=True)
sorted_df.to_csv(kaiju_merge, sep='\t', index=False, header=False)


def paste_files(patterns, output_file_path):
	with open(output_file_path, 'w') as output_file:
		# Find files matching the pattern

		file_paths = []
		for pattern in patterns:
			file_paths.extend(glob.glob(pattern))
		
		# Zip the file iterators together
		for lines in zip(*[open(file) for file in file_paths]):
			# Concatenate the lines and write to the output file
			concatenated_line = '\t'.join(line.strip() for line in lines)
			output_file.write(concatenated_line + '\n')

### horizontal concatenation of the comparison files and the fasta headers of the first file

file_name = os.path.basename(input_files[0])
output_file = os.path.splitext(file_name)[0] + ".allcompare.tmp"
patterns = ['1___*fastaheaders', '2___*fastaheaders.compare.1.tmp', '3___*fastaheaders.compare.1.tmp', '4___*fastaheaders.compare.1.tmp']
paste_files(patterns, output_file)

### sort the allcompare
df = pd.read_csv(output_file, sep='\t', header=None)
sorted_df = df.sort_values(by=[0], ascending=[True])
sorted_df.to_csv(output_file, sep='\t', index=False, header=False)


### sort the kaiju merge
file_name = os.path.basename(input_files[0])
kaiju_merge = os.path.splitext(file_name)[0] + ".kaiju.merge"
df = pd.read_csv(kaiju_merge, sep='\t', header=None)
sorted_df = df.sort_values(by=[1], ascending=[True])
sorted_df.to_csv(kaiju_merge, sep='\t', index=False, header=False)

### horizontal concatenation of the comparison files and the fasta headers of the first file, and kaiju merge
file_name = os.path.basename(input_files[0])
output_file = os.path.splitext(file_name)[0] + ".allcompare.kaiju"

patterns = ['*.allcompare.tmp', '*.kaiju.merge']
paste_files(patterns, output_file)


## get names for the files used as input to be headers
file_patterns = ['1___*fastaheaders', '2___*fastaheaders.compare.1.tmp', '3___*fastaheaders.compare.1.tmp', '4___*fastaheaders.compare.1.tmp']
output_file = "read_header"

with open(output_file, 'w') as f:
	for pattern in file_patterns:
		file_names = glob.glob(pattern)
		for file_name in file_names:
			f.write(file_name + '\t')
	f.write('\n')

## fix header names

# Read the content of the file
with open('read_header', 'r') as file:
	file_content = file.read()
# Remove ".fastaheaders" and ".fastaheaders.compare.1.tmp" from the content using replace()
modified_content = file_content.replace(".fastaheaders", "").replace(".compare.1.tmp", "")
# Write the modified content back to the file
with open('read_header', 'w') as file:
	file.write(modified_content)


###### add kaiju header info
# Read the content of the file
with open('read_header', 'r') as file:
	lines = file.readlines()

modified_lines = [line.strip() + "\tclassified_unclassified\tread_name\tNCBI_taxid_best\tkaiju_match_score\tNCBI_taxid_all\tNCBI_accession_all\tmatching_seq_frag\tNCBI_full_taxon_lineage\n" for line in lines]

# Write the modified lines back to the file
with open('read_header', 'w') as file:
	file.writelines(modified_lines)

##### add header to final file
# Define the names of the input files
file1_name = 'read_header'
file2_pattern  = '*.allcompare.kaiju'

# Define the name of the output (concatenated) file
output_file_name = 'bioinf_readcheck_allcompare.final'

# Read the content of the first file
with open(file1_name, 'r') as file1:
	content_file1 = file1.read()

# Initialize an empty string to store the concatenated content
concatenated_content = content_file1

# Find files matching the glob pattern for file2
file2_matches = glob.glob(file2_pattern)

# Iterate through matching files and concatenate their content
for file2_name in file2_matches:
	with open(file2_name, 'r') as file2:
		content_file2 = file2.read()
	concatenated_content += content_file2

# Write the concatenated content to the output file
with open(output_file_name, 'w') as output_file:
	output_file.write(concatenated_content)

## make csv file and replace possible false delimiters

# Read the content of the file
with open('bioinf_readcheck_allcompare.final', 'r') as file:
	file_content = file.read()
# Remove ".fastaheaders" and ".fastaheaders.compare.1.tmp" from the content using replace()
modified_content = file_content.replace("; ", ";").replace(",", ";").replace("\t", ",")
# Write the modified content back to the file
with open('bioinf_readcheck_allcompare.final', 'w') as file:
	file.write(modified_content)


##### run the R script
subprocess.call("bioinf-read-check.r")