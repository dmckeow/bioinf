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
library(qgraph)

# Import data
ref.metadata <- read_sheet("https://docs.google.com/spreadsheets/d/14BDmFgMYqHksvIRgdRG_jLISXi7MqqAoBVyAn96cm8g/edit#gid=1598046758")
sam.metadata <- read_sheet("https://docs.google.com/spreadsheets/d/1NhxQtQRc7nGmX0t0FwOZJAziMn3r0FshtpehCa-Lo5U/edit#gid=476945537")

# reduce to relevant metadata columns only
FilterMD <- function(x) {
    x %>% select(contig, sample_info1, sample_info2, sample_info3, SuperBin, RepresentativeName, contig_length, miuvig_quality, Project, Sample_metadata_code, BioRep, genus, species, collection_month, collection_year, apiary, distance, flower_genus, flower_species, num_of_bees)
}

HostFilter <- function(data){
      data %>%
      filter((genus == "Apis" | (genus == "Bombus" & species == "impatiens") | Project != "Current_study"))
}

HostFilter_nw <- function(data){
  data %>%
    filter((genus.x == "Apis" | (genus.x == "Bombus" & species.x == "impatiens") | Project.x != "Current_study") &
           (genus.y == "Apis" | (genus.y == "Bombus" & species.y == "impatiens") | Project.y != "Current_study"))
}


ref.metadata <- FilterMD(ref.metadata)
ref.metadata$tip.name <- paste(ref.metadata$RepresentativeName,ref.metadata$contig)

sam.metadata <- FilterMD(sam.metadata)
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
ref.sam.metadata$Project <<- ifelse(ref.sam.metadata$Project != "NCBI", "Current_study", ref.sam.metadata$Project)
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

############################################################################################################
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

concatenated_data <- HostFilter(concatenated_data)

# Write the concatenated data to a new CSV file
write.table(concatenated_data, "All.NoMmseqs.metadata.cytoscape.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

#################################################################
#################################################################

file.remove("All.NoMmseqs.network.cytoscape.tsv")
file_list <- list.files(pattern="*.NoMmseqs.network.cytoscape.tsv", full.names=TRUE, recursive=FALSE)
### exclude other file to prevent duplicattion
file_list <- file_list[ !grepl("Dicistroviridae.NoMmseqs..*|Iflaviridae_Iflavirus.NoMmseqs..*", file_list) ]

# Initialize an empty dataframe to store concatenated data
concatenated_data_nw <- data.frame()

# Loop through each file, read it, and concatenate it to the dataframe
for (file in file_list) {
  data <- read.table(file, header = TRUE, sep = "\t") # Adjust header argument based on your data
  concatenated_data_nw <- rbind(concatenated_data_nw, data)
}

concatenated_data_nw <- concatenated_data_nw %>%
		left_join(ref.sam.metadata, by = c('from' = 'contig')) %>%
		left_join(ref.sam.metadata, by = c('to' = 'contig'))

concatenated_data_nw <- HostFilter_nw(concatenated_data_nw) %>% select(from, to, weight)

# Write the concatenated data to a new CSV file
write.table(concatenated_data_nw, "All.NoMmseqs.network.cytoscape.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

##############################################################################
############ prep the metadata for the protein blasts
file_list <- list.files(path = "iprscansplit/", pattern="*.blastp", full.names=TRUE, recursive=FALSE)

PrepForCytoscapeBlast <- function() {
	blastp.ref.sam.metadata <<- blastp %>%
		left_join(ref.sam.metadata, by = 'contig') %>%
		HostFilter() %>%
		select(-to, -pident, -length, -mismatch, -gapopen, -qstart, -qend, -sstart, -send, -evalue, -bitscore) %>%
		distinct
}

PrepForCytoscapeBlast_nw <- function() {
	blastp.ref.sam.nw <<- blastp %>%
		left_join(ref.sam.metadata, by = c('contig_from' = 'contig')) %>%
		left_join(ref.sam.metadata, by = c('contig_to' = 'contig')) %>%
		HostFilter_nw() #%>%
		#select(from, to, pident)
}


for (file in file_list) {
	output_name <- gsub("iprscansplit/|\\.blastp", "", file)
	blastp <- read.table(file, header = TRUE, sep = "\t")
	blastp$contig <- gsub("_[0-9]+:[0-9]+-[0-9]+$", "", blastp$from)
	PrepForCytoscapeBlast()
 	write.table(blastp.ref.sam.metadata, paste0(output_name, ".blastp.metadata.cytoscape.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")

 	blastp$contig_from <- gsub("_[0-9]+:[0-9]+-[0-9]+$", "", blastp$from)
 	blastp$contig_to <- gsub("_[0-9]+:[0-9]+-[0-9]+$", "", blastp$to)
 	PrepForCytoscapeBlast_nw()
 	write.table(blastp.ref.sam.nw, paste0(output_name, ".blastp.network.cytoscape.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
}



########### plot with igraph, mostly works, but abandoned due to lack of way to label clusters



PlotNetwork2 <- function(INPUT_DATA, METADATA, color_var, color_mapping, VERTEX_SIZE, LAYOUT_AREA_FACTOR, REPULSION_FACTOR, MAX_DISPLACE_FACTOR) {
	set.seed(4)
	network <- graph_from_data_frame(INPUT_DATA, vertices = METADATA)
	isolated = which(degree(network)==0) ## get isolated vertices
	network = delete_vertices(network, isolated) ## remove the isolated verteices
	layout = layout_with_fr(network) ## get network layout without the isolated vertices
	layout = layout[-isolated,] ## get network layout without the isolated vertices

	network <- simplify(network, remove.loops = TRUE, remove.multiple = TRUE)
	

	# Filter metadata to keep only the retained vertices
    metadata_kept <- METADATA[match(V(network)$name, METADATA$from), ]

	### Extract unique categories from the COLORS
    col_categories <- unique(metadata_kept[[color_var]])

    # Generate colors for each category
    #num_col_categories <- length(col_categories)
    #palette_net <- brewer.pal(num_col_categories, "Paired")
    #palette_net <- c("#E31A1C", "#1F78B4", "#FB9A99", "#B2DF8A", "#A6CEE3", "#33A02C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
	
    # Create a mapping from category to color for vertices present in the network
    #color_mapping <- setNames(palette_net, col_categories)

    # Ensure consistent order of levels in the factor
    metadata_kept[[color_var]] <- factor(metadata_kept[[color_var]], levels = col_categories)

    # retain only color mapping in dataset
    color_mapping <- color_mapping[intersect(names(color_mapping), levels(metadata_kept[[color_var]]))]

    # Reorder color_mapping based on the levels of color_var
    color_mapping <- color_mapping[levels(metadata_kept[[color_var]])]

    # Assign colors to vertices based on their category
    V(network)$color <- color_mapping[metadata_kept[[color_var]]]

 	
e <- as_edgelist(network,names=FALSE)
#layout <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(network))

layout <- qgraph.layout.fruchtermanreingold(
	e,
	vcount=vcount(network),
  area=LAYOUT_AREA_FACTOR*(vcount(network)^2),
  repulse.rad=(vcount(network)^REPULSION_FACTOR),
  max.delta = (vcount(network))/MAX_DISPLACE_FACTOR
  ) # 8, 3.1 looked ok

  #layout <- layout_with_fr(network)
  #layout <- norm_coords(layout, ymin=-1, ymax=1, xmin=-1, xmax=1)*2


	plot.igraph(network, 
    	vertex.size=VERTEX_SIZE,
    	vertex.label=NA,
    	vertex.frame.color="black",
    	edge.width=1,
    	edge.arrow.size = 0,
    	margin = 0,
    	layout=layout,
    	#rescale=FALSE
    	)

	legend("topleft", 
      	 legend=paste(levels(as.factor(metadata_kept[[color_var]])), sep=""), 
      	 col = color_mapping, 
      	 #col = palette_net, 
      	 bty = "n", pch=20 , pt.cex = 1, cex = 0.5,
      	 text.col="black" , horiz = F)
	return(network)
    }


### make a palette manually
palette_genus_project <- c("#BB0D0E", "#FB9A99", "#B2DF8A", "#1F78B4", "#A6CEE3")
names(palette_genus_project) <- c("Apis_Current_study", "Apis_NCBI", "Other_NCBI", "Bombus_Current_study", "Bombus_NCBI")



all_RdRp.blastp <- read.csv("iprscansplit/all_RdRp.blastp", header=TRUE, sep="\t")
all_RdRp.blastp <- all_RdRp.blastp %>%
		select(from, to, pident) %>%
		filter(pident >= 95)

all_RdRp.blastp_md <- read.csv("all_RdRp.blastp.metadata.cytoscape.tsv", header=TRUE, sep="\t")

all_RdRp.blastp_md$genus_project <- paste(all_RdRp.blastp_md$genus, all_RdRp.blastp_md$Project, sep = "_")

## nrows, ncols
#par(mfrow=c(1, 2))

#VIRUS = "Iflavirus aladeformis"
#all_RdRp.blastp_md_test <- all_RdRp.blastp_md %>% filter(str_detect(RepresentativeName, VIRUS))
#all_RdRp.blastp_test <- all_RdRp.blastp %>% filter(from %in% all_RdRp.blastp_md_test$from & to %in% all_RdRp.blastp_md_test$from)
#network <- PlotNetwork2(all_RdRp.blastp_test, all_RdRp.blastp_md_test, "genus_project", palette_genus_project)

#par(mfrow=c(1, 1))

network <- PlotNetwork2(all_RdRp.blastp, all_RdRp.blastp_md, "genus_project", palette_genus_project, 1.5, 5, 2.7, 5)

png(paste0("test_network", ".png"), width=48, height=48, res=300, unit="cm")
network <- PlotNetwork2(all_RdRp.blastp, all_RdRp.blastp_md, "genus_project", palette_genus_project)
dev.off()
pdf(paste0("test_network", ".pdf"), width=7, height=6)
network <- PlotNetwork2(all_RdRp.blastp, all_RdRp.blastp_md, "genus_project", palette_genus_project)

dev.off()