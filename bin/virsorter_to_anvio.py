#!/usr/bin/env python
import collections
import argparse
import sys
import os
import csv
import glob

### A useful function
def load_hallmark_functions(db_dir,hash):
	for hallmark_gene_list in glob.glob(os.path.join(db_dir,'group/*/hallmark-gene.list')):
		print("... reading {}".format(hallmark_gene_list))
		with open(hallmark_gene_list, newline='') as csvfile:
			hallmark_reader = csv.reader(csvfile,delimiter='\t',quotechar='"')
			for row in hallmark_reader:
				if row[0].startswith("#"): continue
				if row[0] == "": continue
				### Note - MAY HAVE TO CLEAN UP row[0] (which has faa)
				hash[row[0]]=row[1]


"""
Parse viral contig prediction results from VirSorter and generates three files usable by Anvi'o:
1. "virsorter_additional_info.txt" which you can import as an additional data file for splits. This can be helpful while binning using anvi-interactive or anvi-refine to identify contigs that are phages.
2. "virsorter_collection.txt" which is a collection file that will automatically group splits belonging to each VirSorter phage prediction into its own bin. This can be imported using anvi-import-collection.
3. "virsorter_annotations.txt" which is an annotations file containing predicted functions for hallmark genes. This can be imported using anvi-import-functions.
"""
"""
Parse viral contig prediction results from VirSorter and generates three files usable by Anvi'o:
1. "virsorter_additional_info.txt" which you can import as an additional data file for splits. This can be helpful while binning using anvi-interactive or anvi-refine to identify contigs that are phages.
2. "virsorter_collection.txt" which is a collection file that will automatically group splits belonging to each VirSorter phage prediction into its own bin. This can be imported using anvi-import-collection.
3. "virsorter_annotations.txt" which is an annotations file containing predicted functions for hallmark genes. This can be imported using anvi-import-functions.
"""
parser = argparse.ArgumentParser(description="Parses VirSorter predictions for Anvi\'o.")
parser.add_argument("-i","--vs2_dir", help = "REQUIRED INPUT. Output directory from VirSorter 2\n\n")
parser.add_argument("-s","--splits-info", help = "REQUIRED INPUT. splits_basic_info.txt file.\n\n")
parser.add_argument("-n","--anvio-gene-calls", help = "OPTIONAL INPUT, but REQUIRED if you want an output functions file. Anvi'o gene calls exported from the contigs database")
parser.add_argument("-d","--db_dir", help = "OPTIONAL INPUT, but REQUIRED if you want an output functions file. Path to the VirSorter2 database directory, to extract the functions associated with hallmark genes.")
parser.add_argument("-A","--addl-info", default = "virsorter_additional_info.txt", help = "OUTPUT. Additional info file. Default = \"virsorter_additional_info.txt\". You can import this as an additional data file for splits. This can be helpful while binning using anvi-interactive or anvi-refine to identify contigs that are phages.")
parser.add_argument("-C","--phage-collection", default = "virsorter_collection.txt", help = "OUTPUT. An Anvi\'o collections file with splits for each phage gathered into a separate bin. The default name is \"virsorter_collection.txt\".")
parser.add_argument("-F","--output-functions", default = "virsorter_annotations.txt", help = "OPTIONAL OUTPUT. Annotations for hallmark genes in VirSorter-predicted phages / prophages. This file is importable using anvi-import-functions.\n\n")
parser.add_argument("-L","--min-phage-length", type = int, default = 1000, help = "PARAMETER, AFFECTS ALL OUTPUTS. Specify the minimum phage length to report. Default is 1000 bp.")
parser.add_argument("--min_score", action = "store_true", help = "PARAMETER, AFFECTS ALL OUTPUTS. Excludes all predictions with score lower than this parameter [0-1].")
parser.add_argument("--exclude-prophages", action = "store_true", help = "PARAMETER, AFFECTS ALL OUTPUTS. Exclude all prophage predictions.")
args = parser.parse_args()
arg_dict = vars(args)

if arg_dict['vs2_dir'] == None or arg_dict['splits_info'] == None:
	print("\n***A required option is missing. Try again.***\n\n")
	parser.print_help()
	sys.exit()
if (arg_dict['anvio_gene_calls'] == None and arg_dict['db_dir'] != None) or (arg_dict['anvio_gene_calls'] != None and arg_dict['db_dir'] == None):
	print("\n***You only specified one of the files required to generate an annotations file for predicted phages.\n")
	parser.print_help()
	sys.exit()

## Listing files from VirSorter2, and checking they are all here
arg_dict['global_file'] = os.path.join(arg_dict['vs2_dir'],"final-viral-score.tsv")
arg_dict['boundary_file'] = os.path.join(arg_dict['vs2_dir'],"final-viral-boundary.tsv")
arg_dict['affi_file'] = os.path.join(arg_dict['vs2_dir'],"for-dramv/viral-affi-contigs-for-dramv.tab")
if not(os.path.exists(arg_dict['global_file'])):
	print("\n*** problenm -> {} does not exit, and we need it, sorry.".format(arg_dict['global_file']))
	sys.exit()
if not(os.path.exists(arg_dict['boundary_file'])):
	print("\n*** problenm -> {} does not exit, and we need it, sorry.".format(arg_dict['boundary_file']))
	sys.exit()
if not(os.path.exists(arg_dict['global_file'])):
	print("\n*** problenm -> {} does not exit, and we need it, sorry. Did you run VirSorter2 with the option '--prep-for-dramv' ? If not, you should, because we need this file.".format(arg_dict['affi_file']))
	sys.exit()


#PART ZERO
#Here I need to read in functional annotations for each hallmark phage cluster
f_dict = collections.defaultdict(dict)
if arg_dict['db_dir'] != None:
	load_hallmark_functions(arg_dict['db_dir'],f_dict)


#PART ONE

## Pre-load boundaries if we take prophage
print("######################## Reading all relevant VirSorter2 data")
info_contigs = {} ## Store prediction info for each contig, which will also get us a dictionary of contigs we are interested in
with open(arg_dict['global_file'], newline='') as csvfile:
	global_reader = csv.reader(csvfile,delimiter='\t',quotechar='"')
	col_score = -1
	for row in global_reader:
		if row[0].startswith("#"): continue
		if row[0] == "": continue
		if row[0] == "seqname":
			for i in range(2, 9):
				if row[i] == "max_score":
					col_score = i
			continue
		if col_score == -1: ## We need this because, depending on the models selected, the score column is variable, and we don't want to add a pandas dependency when it's not 100% needed
			print("We had some problem reading the header of the VirSorter2 score file, sorry ... ")
			sys.exit()
		vir_tab = row[0].split("||")
		contig = vir_tab[0]
		if vir_tab[1]=="full" or vir_tab[1]=="lt2gene" or not(arg_dict['exclude_prophages']):
			## We are interested by this prediction
			if contig not in info_contigs:
				info_contigs[contig] = {} ## Initiate the table of prediction if this is the first time we see this contig
			info_contigs[contig][row[0]] = {} ## Initiate the table for this specific prediction
			info_contigs[contig][row[0]]["score"] = row[col_score]
			if row[col_score] == "nan":
				info_contigs[contig][row[0]]["score"] = 0 ## We set the "nan" to 0, this is not super pretty, but otherwise Anvi'o does not interpret the column as a float, and that's annoying
			info_contigs[contig][row[0]]["length"] = row[col_score+2]
			info_contigs[contig][row[0]]["nb_hallmark"] = row[col_score+3]
			if vir_tab[1]=="full" or vir_tab[1]=="lt2gene":
				## Full predictionsif contig not in info_contigs:
				info_contigs[contig][row[0]]["type"] = "full"
				if vir_tab[1]=="lt2gene":
					info_contigs[contig][row[0]]["nb_genes"] = row[col_score+3] ### the "lt2gene" are not listed in the "boundary" files which is where the total # of genes is indicated. Instead, we use the total # of hallmark genes as total number of genes for these short contigs (doesn't really matter, this nb_genes is 1 or 2 anyway)
			elif not(arg_dict['exclude_prophages']):
				info_contigs[contig][row[0]]["type"] = "prophage"
#####################
with open(arg_dict['boundary_file'], newline='') as csvfile:
	boundary_reader = csv.reader(csvfile,delimiter='\t',quotechar='"')
	for row in boundary_reader:
		if row[0].startswith("#"): continue
		if row[0] == "": continue
		if row[0] == "seqname": continue
		virus = row[28]
		virtab = virus.split("||")
		contig = virtab[0]
		if contig not in info_contigs: continue ## We did not select this contig
		if virus not in info_contigs[contig]: continue ## We did not select this prediction
		info_contigs[contig][virus]["nb_genes"] = row[15]
		if virtab[1]=="full" or virtab[1]=="lt2gene": ## We don't need to know the boundaries for predictions that include the whole contig
			info_contigs[contig][virus]["topology"] = row[27] ## But we can take the topology (i.e. if the contig is circular or linear according to VS2 )
		else:
			info_contigs[contig][virus]["start"] = row[3] ## This is a prophage, so we take the coordinates, and we don't take the contig's topology because it's not relevant for prophages
			info_contigs[contig][virus]["end"] = row[4]


#Next, I need to read in virsorter_affi, to get the gene prediction and annotation performed by VirSorter, which we will try to match to the ones done by Anvi'o
## Read affi
# If contig is of interest
## If full -> Take the gene
## If prophage -> Check if interesting gene, if so, take the gene

affi_dict = {}
affi_dict = collections.defaultdict(dict)
hallmark_dict = collections.defaultdict(dict)
c_contig = ""
c_virus = ""
with open(arg_dict['affi_file'], newline='') as csvfile:
	affi_reader = csv.reader(csvfile,delimiter='|',quotechar='"')
	for row in affi_reader:
		if row[0].startswith(">"):
			row[0] = row[0].replace(">","")
			virtab = row[0].split("__")
			if len(virtab)>2:
				last = len(virtab) - 2
				virtab[0] = "__".join(virtab[0:last])
				virtab[1]=virtab[last+1]
			c_contig = virtab[0]
			if c_contig in info_contigs:
				c_virus = c_contig + "||" + virtab[1]
				affi_dict[c_virus] = {}
				hallmark_dict[c_virus] = {}
		elif c_contig in info_contigs:
			virtab = row[0].split("__")
			num_gene = virtab[len(virtab) - 1]
			gene_id = c_virus+"||"+str(num_gene)
			start = row[1]
			stop = row[2]
			strand = row[4]
			cat = row[8]
			affi_dict[c_virus][gene_id] = {}
			affi_dict[c_virus][gene_id]['start'] = start
			affi_dict[c_virus][gene_id]['stop'] = stop
			affi_dict[c_virus][gene_id]['strand'] = strand
			affi_dict[c_virus][gene_id]['cat'] = cat
			if cat == str(0) or cat == str(3):
				#This next if statement deals with the fact that VirSorter doesn't count hallmark genes
				#if a hallmark gene has a PFAM hit with a better e-value than the phage cluster hit e-value.
				hit = row[5]
				score = row[6]
				pfam_score = row[10]
				affi_dict[c_virus][gene_id]['hallmark_function'] = f_dict[hit]
				hallmark_dict[c_virus][gene_id] = {}
				hallmark_dict[c_virus][gene_id]['hallmark_function'] = f_dict[hit]
				hallmark_dict[c_virus][gene_id]['phage_cluster'] = hit
				hallmark_dict[c_virus][gene_id]['start'] = start
				hallmark_dict[c_virus][gene_id]['stop'] = stop
				hallmark_dict[c_virus][gene_id]['cat'] = cat
				hallmark_dict[c_virus][gene_id]['score'] = score
				hallmark_dict[c_virus][gene_id]['pfam_score'] = pfam_score


#PART THREE
#This writes all phages and prophages, to additional-info and collections files.
splits_input = open(arg_dict['splits_info'], "r")
splits_output = open(arg_dict['addl_info'],'w')
collection_output = open(arg_dict['phage_collection'],'w')

#splits_basic_info column format:
#  0            1             2     3      4          5               6              7
#split   order_in_parent   start   end   length   gc_content   gc_content_parent   parent

#is_prev_phage = False
#n = 1

splits_output.write("split\tphage_name\tphage_score\tphage_length\tphage_pred_type\tnum_genes_in_phage\tnum_phage_hallmark_genes_in_phage\n")
line = splits_input.readline()
line = splits_input.readline()

split_length_dict = collections.defaultdict(dict)
while line != "":
	line = line.strip().split('\t')
	if 'length' not in split_length_dict[line[7]]:
		split_length_dict[line[7]]['length'] = 0
	split_length_dict[line[7]]['length'] += int(line[4])
	line = splits_input.readline()
splits_input.close()

print("######################## Printing collection file and additional data file")
splits_input = open(arg_dict['splits_info'], "r")
line = splits_input.readline()
line = splits_input.readline()
while line != "":
	line = line.strip().split('\t')
	split_name = line[0]
	split_start = int(line[2])
	split_stop = int(line[3])
	split_parent = line[7]
	split_length = line[4]
	parent_length = split_length_dict[split_parent]['length']
	if split_parent in info_contigs:
		virtab = list(info_contigs[split_parent].keys())
		if (len(virtab) == 1) & (info_contigs[split_parent][virtab[0]]['type'] == "full"):
			## This is the easy case - the whole contig was predicted as viral
			virus = virtab[0]
			virus_name = virus
			virus_name = virus.replace("||","__")
			nb_hallmark = info_contigs[split_parent][virus]['nb_hallmark']
			if ('nb_genes' not in  info_contigs[split_parent][virus]):
				print("{} - {} doe not have a n_genes ?".format(split_parent, virus))
			nb_genes = info_contigs[split_parent][virus]['nb_genes']
			score = info_contigs[split_parent][virus]['score']
			type = "Full"
			virus_length = info_contigs[split_parent][virus]['length']
			if split_length < virus_length:
				# print("Split_length {} is shorter than virus length {} for split {} and virus {} - we adjust".format(split_length,virus_length,split_name,virus))
				virus_length = split_length
			splits_output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (split_name, virus_name, str(score), str(virus_length), type, str(nb_genes), str(nb_hallmark)))
			collection_output.write("%s\t%s\n" %(split_name, virus_name))
		else:
			## Now we have prophage prediction(s), so we need to go through these
			for prophage in info_contigs[split_parent]:
				nb_hallmark = info_contigs[split_parent][prophage]['nb_hallmark']
				nb_genes = info_contigs[split_parent][prophage]['nb_genes']
				score = info_contigs[split_parent][prophage]['score']
				virus_length = info_contigs[split_parent][prophage]['length']
				type = "Prophage"
				start_pos = int(info_contigs[split_parent][prophage]['start'])
				end_pos = int(info_contigs[split_parent][prophage]['end'])
				if start_pos >= split_start and end_pos <= split_stop:
					## Prophage all in this split
					prophage_name = prophage.replace("||","__") + "_" + str(start_pos) + "_" + str(end_pos) ## Include prophage boundaries in prophage name
					splits_output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (split_name, prophage_name, str(score), str(virus_length), type, str(nb_genes), str(nb_hallmark)))
					collection_output.write("%s\t%s\n" %(split_name, prophage_name))
				elif (start_pos >= split_start and start_pos <= split_stop) or (end_pos >= split_start and end_pos <= split_stop):
					## Prophage partially in this split
					prophage_name = prophage.replace("||","__") + "_" + str(start_pos) + "_" + str(end_pos) ## Add prophage boundaries to name
					# if start_pos < split_start: start_pos = split_start ## Calculate reduced boundaries in the split
					# if end_pos > split_stop: end_pos = split_stop
					# prophage_name = prophage_name + "_" + str(start_pos) + "_" + str(end_pos) ## Add prophage boundaries in the split to the same name
					virus_length = end_pos - start_pos + 1
					splits_output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (split_name, prophage_name, str(score), str(virus_length), type, str(nb_genes), str(nb_hallmark)))
					collection_output.write("%s\t%s\n" %(split_name, prophage_name))
	line = splits_input.readline()
splits_input.close()
splits_output.close()
collection_output.close()


#Part 4a:
#Here I need to read in gene calls from Anvi'o
#File Format:
#gene_callers_id contig  start   stop    direction       partial source  version
#We just read in contigs to gene_dict that were annotated as being all or partially viral.
gene_dict = collections.defaultdict(dict)
if arg_dict['anvio_gene_calls'] != None:
	print("######################## Now working on the predicted gene functions")
	print("Loading gene calls from Anvi'o")
	with open(arg_dict['anvio_gene_calls'], newline='') as csvfile:
		anvio_reader = csv.reader(csvfile,delimiter='\t',quotechar='"')
		for line in anvio_reader:
			if line[1] in info_contigs:
				gcid = str(line[0])
				contig = line[1]
				start = line[2]
				stop = line[3]
				direction = line[4]
				gene_dict[contig][gcid] = {}
				gene_dict[contig][gcid]['start'] = start
				gene_dict[contig][gcid]['stop'] = stop
				gene_dict[contig][gcid]['direction'] = direction
	#Part 4b
	#I need to read in the splits file I just generated to get a set of just the contigs we're reporting on here
	#rather than integrating the below code into the many cases in part 3 above.
	#I'll only export functions for contigs/regions declared interesting in the command line.
	#open additional data file, ead in first column splitting on 'split', add to set
	contigs_set = set()
	print("Loading phage list")
	vir_contigs = open(arg_dict['phage_collection'],"r")
	line = vir_contigs.readline()
	while line != "":
		line = line.strip().split('\t')
		contig = line[0].split('_split_')[0]
		contigs_set.add(contig)
		line = vir_contigs.readline()
	vir_contigs.close()

	#write output functions file
	print("Writing the output function file")
	func_out = open(arg_dict['output_functions'],"w")
	func_out.write("gene_callers_id\tsource\taccession\tfunction\te_value\n")
	tmp_matched = 0
	tmp_total = 0
	matched = 0
	total = 0
	for virus in hallmark_dict:
		#if contig in additional data file:... tab the stuff below after turhing this on
		tmp = virus.split("||")
		contig = tmp[0]
		tmp_matched = 0
		tmp_total = 0
		tmp_test = 0
		for gene in hallmark_dict[virus]:
			tmp_total += 1
			h_gene_start = int(hallmark_dict[virus][gene]['start'])
			h_gene_stop = int(hallmark_dict[virus][gene]['stop'])
			for gcid in gene_dict[contig]:
				tmp_test += 1
				if abs(h_gene_start-int(gene_dict[contig][gcid]['start']))<2 or abs(h_gene_stop-int(gene_dict[contig][gcid]['stop']))<2:
					tmp_matched += 1
					pc = hallmark_dict[virus][gene]['phage_cluster']
					hf = hallmark_dict[virus][gene]['hallmark_function']
					score = hallmark_dict[virus][gene]['score']
					if hallmark_dict[virus][gene]['cat'] == str(3):
						hf = hf+(" (non-Caudovirales)")
					func_out.write("%s\t%s\t%s\t%s\t%s\n" % (str(gcid), "VirSorter 2 Viral hallmark genes", pc, hf, str(score)))
					break
		matched += tmp_matched
		total += tmp_total
		if tmp_total > 0:
			print("{} of {} hallmark gene(s) matched for virus {} - {} tested for contig {}".format(tmp_matched, tmp_total, virus, tmp_test, contig))
	func_out.close()
	print("{} hallmark genes matched total out of {}".format(matched, total))
