#!/usr/bin/env Rscript

#### setwd("C:/Users/Dean Mckeown/Downloads/Spillover_FINAL/phylogeny")

library(ggrepel)
library(googlesheets4)
library(tidyverse)
library(viridis)
library(TDbook)
library(Biostrings)
library(igraph)
library(RColorBrewer)
library(ggraph)

# Import data
ref.metadata <- read_sheet("https://docs.google.com/spreadsheets/d/14BDmFgMYqHksvIRgdRG_jLISXi7MqqAoBVyAn96cm8g/edit#gid=1598046758")
sam.metadata <- read_sheet("https://docs.google.com/spreadsheets/d/1kzg8noVgySe9xLxQBW2BXR1ospcgSqcZAIQNwI8Ygn4/edit#gid=1459467912")

# reduce to relevant metadata columns only
FilterMD <- function(x) {
    x %>% select(contig, sample_info1, sample_info2, sample_info3, SuperBin, RepresentativeName, contig_length, miuvig_quality, Project, Sample_metadata_code, BioRep, genus, species, collection_month, collection_year, apiary, distance, flower_genus, flower_species, num_of_bees)
}

ref.metadata <- FilterMD(ref.metadata)
sam.metadata <- FilterMD(sam.metadata)
ref.metadata$tip.name <- paste(ref.metadata$RepresentativeName,ref.metadata$contig)
sam.metadata$tip.name <- NA

#####################################################################################################
#####################################################################################################

PrepMetadata <- function() {

## reshape metadata
ref.sam.metadata <<- rbind(ref.metadata, sam.metadata)
ref.sam.metadata <<- lapply(ref.sam.metadata, as.character)
ref.sam.metadata <<- data.frame(ref.sam.metadata)
ref.sam.metadata$contig_length <<- as.numeric(ref.sam.metadata$contig_length)
ref.sam.metadata$genus <<- ifelse(ref.sam.metadata$genus != "Apis" & ref.sam.metadata$genus != "Bombus", "Other", ref.sam.metadata$genus)
ref.sam.metadata$genus <<- ifelse(is.na(ref.sam.metadata$genus), "Other", ref.sam.metadata$genus)
ref.sam.metadata$Project <<- gsub("GenBank","NCBI",ref.sam.metadata$Project)
ref.sam.metadata$Project <<- ifelse(ref.sam.metadata$Project != "NCBI", "Current study", ref.sam.metadata$Project)
ref.sam.metadata$Project <<- ifelse(is.na(ref.sam.metadata$Project), "NCBI", ref.sam.metadata$Project)
}

ReadMatrix <- function(input_matrix) {
# Read the lines from the file
	lines <- readLines(input_matrix)
	lines <- lines[-1]

# Split each line by whitespace and convert "NA" strings to actual NA values
	split_lines <- lapply(strsplit(lines, "\\s+"), function(x) { x[x == "NA"] <- NA; x })

# Find the maximum number of elements in any line
	max_length <- max(sapply(split_lines, length)) + 1

# Fill the shorter lines with NA values to make them equal length
	filled_lines <- lapply(split_lines, function(x) { c(x, rep(NA, max_length - length(x))) })

# Convert the filled lines into a matrix	
	fastANI_matrix <<- as.data.frame(do.call(rbind, filled_lines))
	fastANI_matrix[, 1] <<- gsub(".*/", "", fastANI_matrix[, 1])
	fastANI_matrix[, 1] <<- gsub(".fa$", "", fastANI_matrix[, 1])
	fastANI_matrix <<- fastANI_matrix %>% 
		column_to_rownames("V1")
	colnames(fastANI_matrix) <<- rownames(fastANI_matrix)
	fastANI_matrix[is.na(fastANI_matrix)] <<- 0 ## we keep the zeroes in the data, but we will remove unconnected vertices when we plot
	fastANI_matrix <<- as.matrix(fastANI_matrix)
}

PlotNetwork <- function(color_var) {
	set.seed(4)
	network <<- graph_from_adjacency_matrix(fastANI_matrix, weighted=T, mode="lower", diag=F)
	isolated = which(degree(network)==0) ## get isolated vertices
	network = delete_vertices(network, isolated) ## remoev the isolated verteices
	layout = layout_with_fr(network) ## get network layout without the isolated vertices
	layout = layout[-isolated,] ## get network layout without the isolated vertices

	is_simple(simplify(network, remove.loops = TRUE, remove.multiple = TRUE))
	cluster_label_prop(network)

	network_df <<- as_data_frame(network)
	ref.sam.metadata.sample <<- network_df %>%
		left_join(ref.sam.metadata, by = c('from' = 'contig')) %>%
		select(-to, -weight)

	# Map the color to color_var
	coul <<- brewer.pal(nlevels(as.factor(color_var)), "Set1")
	my_color <<- coul[as.numeric(as.factor(color_var))]
	par(bg="black", mar=c(0,0,0,0))
	sn_color_range <- colorRampPalette(c("grey50","grey40","grey30","grey20"))
	sn_color <- sn_color_range(length(network_df$weight))

	plot.igraph(network, 
    	vertex.size=4,
    	vertex.color=my_color,
    	vertex.label=NA,
    	vertex.frame.color="transparent",
    	edge.color=sn_color,
    	edge.width=1
    	)
	legend("bottomleft", 
      	 legend=paste(levels(as.factor(color_var)), sep=""), 
      	 col = coul , 
      	 bty = "n", pch=20 , pt.cex = 1, cex = 0.5,
      	 text.col="white" , horiz = F)

   }

PrepForCytoscape <- function() {
	network_df <<- as.data.frame(as.table(fastANI_matrix))
	colnames(network_df) <<- c("from", "to", "weight")
	ref.sam.metadata.sample <<- network_df %>%
		left_join(ref.sam.metadata, by = c('from' = 'contig')) %>%
		select(-to, -weight)
	network_df <<- network_df[network_df$weight != 0, ]
}

PrepMetadata()

################# run as loop ove all matrices available

#file_list <- list.files(path = "fastANI/", pattern="Apis_rhabdovirus.*.matrix", full.names=TRUE, recursive=FALSE)

#for (file in file_list) {
#	output_name <- gsub("fastANI/|\\.matrix", "", file)
# 	ReadMatrix(file)

# 	png(paste0(output_name, ".network", ".png"), width=800, height=800, res=300, unit="px")
# 	PlotNetwork(ref.sam.metadata.sample$genus)
 #	dev.off()

 #	pdf(paste0(output_name, ".network", ".pdf"), width=8, height=8)
 #	PlotNetwork(ref.sam.metadata.sample$genus)
 #	dev.off()
#}



################ OR jsut prep input files for cytoscape
file_list <- list.files(path = "fastANI/", pattern="*.NoMmseqs.matrix", full.names=TRUE, recursive=FALSE)

for (file in file_list) {
	output_name <- gsub("fastANI/|\\.matrix", "", file)
 	ReadMatrix(file)
	PrepForCytoscape()
	write.table(network_df, paste0(output_name, ".network.cytoscape.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
 	write.table(ref.sam.metadata.sample, paste0(output_name, ".metadata.cytoscape.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
}

########## concat all networks into a snigle csv (optional)
file.remove("All.NoMmseqs.metadata.cytoscape.tsv")
file_list <- list.files(pattern="*.NoMmseqs.metadata.cytoscape.tsv", full.names=TRUE, recursive=FALSE)
### exclude other file to prevent duplicattion
file_list <- file_list[ !grepl("Dicistroviridae.NoMmseqs..*|Iflaviridae_Iflavirus.NoMmseqs..*", file_list) ]

# Initialize an empty dataframe to store concatenated data
concatenated_data <- data.frame()

# Loop through each file, read it, and concatenate it to the dataframe
for (file in file_list) {
  data <- read.table(file, header = TRUE, sep = "\t") # Adjust header argument based on your data
  concatenated_data <- rbind(concatenated_data, data)
}

# Write the concatenated data to a new CSV file
write.table(concatenated_data, "All.NoMmseqs.metadata.cytoscape.tsv", row.names = FALSE, quote = FALSE, sep = "\t")


file.remove("All.NoMmseqs.network.cytoscape.tsv")
file_list <- list.files(pattern="*.NoMmseqs.network.cytoscape.tsv", full.names=TRUE, recursive=FALSE)
### exclude other file to prevent duplicattion
file_list <- file_list[ !grepl("Dicistroviridae.NoMmseqs..*|Iflaviridae_Iflavirus.NoMmseqs..*", file_list) ]

# Initialize an empty dataframe to store concatenated data
concatenated_data <- data.frame()

# Loop through each file, read it, and concatenate it to the dataframe
for (file in file_list) {
  data <- read.table(file, header = TRUE, sep = "\t") # Adjust header argument based on your data
  concatenated_data <- rbind(concatenated_data, data)
}

# Write the concatenated data to a new CSV file
write.table(concatenated_data, "All.NoMmseqs.network.cytoscape.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

##############################################################################
############ prep the metadata for the protein blasts
file_list <- list.files(path = "iprscansplit/", pattern="*.blastp", full.names=TRUE, recursive=FALSE)

PrepForCytoscapeBlast <- function() {
	blastp.ref.sam.metadata <<- blastp %>%
		left_join(ref.sam.metadata, by = 'contig') %>%
		select(-to, -pident, -length, -mismatch, -gapopen, -qstart, -qend, -sstart, -send, -evalue, -bitscore) %>%
		distinct
}


for (file in file_list) {
	output_name <- gsub("iprscansplit/|\\.blastp", "", file)
	blastp <- read.table(file, header = TRUE, sep = "\t") # Adjust header argument based on your data
	blastp$contig <- gsub("_[0-9]+:[0-9]+-[0-9]+$", "", blastp$from)
	PrepForCytoscapeBlast()
 	write.table(blastp.ref.sam.metadata, paste0(output_name, ".blastp.metadata.cytoscape.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
}