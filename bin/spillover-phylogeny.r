#!/usr/bin/env Rscript

#### setwd("C:/Users/Dean Mckeown/Downloads/Spillover_FINAL/phylogeny")

library(ape)
library(ggtree)
library(ggrepel)
library(googlesheets4)
library(dplyr)
library(viridis)
library(TDbook)
library(Biostrings)
library(treeio)

#### there is somethnig wrong with ggtree - this code fixes it:
nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")

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

PrepMetaDataTree <- function() {

## reshape metadata
ref.sam.metadata <<- rbind(ref.metadata, sam.metadata)
ref.sam.metadata <<- lapply(ref.sam.metadata, as.character)
ref.sam.metadata <<- data.frame(ref.sam.metadata)
ref.sam.metadata$contig_length <<- as.numeric(ref.sam.metadata$contig_length)
ref.sam.metadata$genus <<- ifelse(ref.sam.metadata$genus != "Apis" & ref.sam.metadata$genus != "Bombus", "Other", ref.sam.metadata$genus)
ref.sam.metadata$Project <<- gsub("GenBank","NCBI",ref.sam.metadata$Project)
ref.sam.metadata$Project <<- ifelse(ref.sam.metadata$Project != "NCBI", "Current study", ref.sam.metadata$Project)

# get names of sequences from fasta
    fasta_seq <- readDNAStringSet(paste0(virus_name, ".aln"))
    deflines <- fasta_seq %>% names() %>% data.frame()
    colnames(deflines) <- c("contigAln")
    deflines$contigAln <- gsub(" .*", "", deflines$contigAln)
    deflines$contig <- gsub("^_R_", "", deflines$contigAln)

# reduce to relvant metarows only - we add it twice to get around weird first column interactiobn with ggtree
    ref.sam.metadata <<- deflines %>%
        left_join(ref.sam.metadata, by='contig')
    ref.sam.metadata <<- deflines %>%
        left_join(ref.sam.metadata, by='contig')

# Read in tree
    tree <<- read.tree(paste0(virus_name, tree_ext))

# Root tree
    #tree <- root(tree, "AF092924.1_Sacbrood_virus", resolve.root=TRUE, edgelabel=TRUE)
}


DrawTree <- function(input_tree) {
    p <<- input_tree %>% 
    ggtree(linewidth=0.9,
           layout="rectangular",
             #branch.length = "none"
    ) %<+% ref.sam.metadata +
    geom_tiplab(aes(label=tip.name), size = 3, color = "black", offset=0.01) + 
    geom_treescale(x=0, y=-5, width=0.1, color='black') +
    geom_point(aes(color=genus, shape=Project), alpha=0.7, na.rm=TRUE, size = 3) +
    scale_color_manual(values = c(Apis = 'red', Bombus = 'blue', Other = 'grey')) +
    geom_text2(aes(label=label,
             subset = !is.na(as.numeric(label)) & as.numeric(label) > 50),
             nudge_x = -0.01, nudge_y = 0.04, 
             check_overlap = TRUE, size = 3)

    ggsave(plot = p, paste0(output_name, ".pdf"), dpi=300, scale=2, units = "cm")
    ggsave(plot = p, paste0(output_name, ".png"), dpi=300, scale=2, units = "cm")

    p + geom_text(aes(label=node), hjust=-.3) # table the nodes to find them
}

CollapseTree <- function(node_to_collapse) {
    p <<- p %>% ggtree::collapse(node=node_to_collapse)
    ggsave(plot = p, paste0(output_name, ".pdf"), dpi=300, scale=2, units = "cm")
    ggsave(plot = p, paste0(output_name, ".png"), dpi=300, scale=2, units = "cm")
    p + geom_text(aes(label=node), hjust=-.3) # table the nodes to find them
}

############### TEMPLATE
#### main tree
#tree_ext <- ".aln.fasttree"
#virus_name <- "test"

#output_name <- "test"
#PrepMetaDataTree()
#DrawTree(tree)

##### subtree
#sub_tree <- tree_subset(tree, node=813, levels_back=0)
#output_name <- "test"
#DrawTree(sub_tree)

###################################################

##### Apis rhabovirus
#### main tree
tree_ext <- ".contree"
virus_name <- "Apis_rhabdovirus"

output_name <- "Apis_rhabdovirus"
PrepMetaDataTree()
DrawTree(tree)

##### Dicistroviridae
#### main tree
tree_ext <- ".aln.fasttree"
virus_name <- "Dicistroviridae"

output_name <- "Dicistroviridae"
PrepMetaDataTree()
DrawTree(tree)

##### subtree
sub_tree <- tree_subset(tree, node=813, levels_back=0)
output_name <- "Dicistroviridae_1"
DrawTree(sub_tree)

sub_tree <- tree_subset(tree, node=676, levels_back=0)
output_name <- "Dicistroviridae_2"
DrawTree(sub_tree)
CollapseTree(493)

output_name <- "Dicistroviridae_3"
PrepMetaDataTree()
DrawTree(tree)
CollapseTree(676)

##### Iflaviridae
## main tree
tree_ext <- ".aln.fasttree"
virus_name <- "Iflaviridae"

output_name <- "Iflaviridae"
PrepMetaDataTree()
DrawTree(tree)

##### subtrees
sub_tree <- tree_subset(tree, node=3214, levels_back=0)
output_name <- "Iflaviridae_1"
DrawTree(sub_tree)

sub_tree <- tree_subset(tree, node=2817, levels_back=9)
output_name <- "Iflaviridae_2"
DrawTree(sub_tree)
CollapseTree(1506)


output_name <- "Iflaviridae_3"
PrepMetaDataTree()
DrawTree(tree)
CollapseTree(2817)


##### Negevirus like
tree_ext <- ".aln.fasttree"
virus_name <- "Negevirus_like"

output_name <- "Negevirus_like"
PrepMetaDataTree()
DrawTree(tree)
##### Partiti like
tree_ext <- ".contree"
virus_name <- "Partiti_like"

output_name <- "Partiti_like"
PrepMetaDataTree()
DrawTree(tree)

##### Phasmaviridae
tree_ext <- ".contree"
virus_name <- "Phasmaviridae"

output_name <- "Phasmaviridae"
PrepMetaDataTree()
DrawTree(tree)

##### Picorna like Mayfield
tree_ext <- ".contree"
virus_name <- "Picorna_like_Mayfield"

output_name <- "Picorna_like_Mayfield"
PrepMetaDataTree()
DrawTree(tree)

##### Reo like
tree_ext <- ".contree"
virus_name <- "Reo_like"

output_name <- "Reo_like"
PrepMetaDataTree()
DrawTree(tree)

##### Sinaivirus
tree_ext <- ".contree"
virus_name <- "Sinaivirus"

output_name <- "Sinaivirus"
PrepMetaDataTree()
DrawTree(tree)

##### Virga like
tree_ext <- ".contree"
virus_name <- "Virga_like"

output_name <- "Virga_like"
PrepMetaDataTree()
DrawTree(tree)








  
  # Lable the A and B clades
  #geom_cladelabel(node=128, label="  DWV-A", color="red", barsize = 1, offset = 1) +
  #geom_cladelabel(node=111, label="  DWV-B", color="blue", barsize = 1, offset = 1) +

  # Add bootstrap values
  #geom_text2(aes(label=label,
             #subset = !is.na(as.numeric(label)) & as.numeric(label) > 50),
             #nudge_x = -0.01, nudge_y = 0.04, 
             #check_overlap = TRUE, size = 3)