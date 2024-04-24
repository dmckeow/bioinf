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
library(patchwork)


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

PrepMetaDataTreeDNA <- function() {

## reshape metadata
ref.sam.metadata <<- rbind(ref.metadata, sam.metadata)
ref.sam.metadata <<- lapply(ref.sam.metadata, as.character)
ref.sam.metadata <<- data.frame(ref.sam.metadata)
ref.sam.metadata$contig_length <<- as.numeric(ref.sam.metadata$contig_length)
ref.sam.metadata$genus <<- ifelse(is.na(ref.sam.metadata$genus) | (ref.sam.metadata$genus != "Apis" & ref.sam.metadata$genus != "Bombus"), "Other", ref.sam.metadata$genus)
ref.sam.metadata$genus <<- ifelse(is.na(ref.sam.metadata$genus), "Other", ref.sam.metadata$genus)
ref.sam.metadata$Project <<- gsub("GenBank","NCBI",ref.sam.metadata$Project)
ref.sam.metadata$Project <<- ifelse(ref.sam.metadata$Project != "NCBI", "Current study", ref.sam.metadata$Project)
ref.sam.metadata$Project <<- ifelse(is.na(ref.sam.metadata$Project), "NCBI", ref.sam.metadata$Project)

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
    phytree <<- read.tree(paste0(virus_name, tree_ext))

# Root tree
    #tree <- root(tree, "AF092924.1_Sacbrood_virus", resolve.root=TRUE, edgelabel=TRUE)
}


PrepMetaDataTreeAA <- function() {

## reshape metadata
ref.sam.metadata <<- rbind(ref.metadata, sam.metadata)
ref.sam.metadata <<- lapply(ref.sam.metadata, as.character)
ref.sam.metadata <<- data.frame(ref.sam.metadata)
ref.sam.metadata$contig_length <<- as.numeric(ref.sam.metadata$contig_length)
ref.sam.metadata$genus <<- ifelse(is.na(ref.sam.metadata$genus) | (ref.sam.metadata$genus != "Apis" & ref.sam.metadata$genus != "Bombus"), "Other", ref.sam.metadata$genus)
ref.sam.metadata$genus <<- ifelse(is.na(ref.sam.metadata$genus), "Other", ref.sam.metadata$genus)
ref.sam.metadata$Project <<- gsub("GenBank","NCBI",ref.sam.metadata$Project)
ref.sam.metadata$Project <<- ifelse(ref.sam.metadata$Project != "NCBI", "Current study", ref.sam.metadata$Project)
ref.sam.metadata$Project <<- ifelse(is.na(ref.sam.metadata$Project), "NCBI", ref.sam.metadata$Project)

# get names of sequences from fasta
    fasta_seq <- readAAStringSet(paste0(virus_name, ".aln"))
    deflines <- fasta_seq %>% names() %>% data.frame()
    colnames(deflines) <- c("prot")
    deflines$prot <- gsub(" .*", "", deflines$prot)
    deflines$contig <- gsub("^_R_", "", deflines$prot)
    deflines$contig <- gsub("_[0-9]+:[0-9]+-[0-9]+$", "", deflines$contig)

# reduce to relvant metarows only 
    ref.sam.metadata <<- deflines %>%
        left_join(ref.sam.metadata, by='contig')
    ref.sam.metadata$prot <<- gsub(":", "_", ref.sam.metadata$prot)

# Read in tree
    phytree <<- read.tree(paste0(virus_name, tree_ext))

# Root tree
    #tree <- root(tree, "AF092924.1_Sacbrood_virus", resolve.root=TRUE, edgelabel=TRUE)
}


DrawTree1 <- function(input_tree) {
    p <<- input_tree %>% 
    ggtree(linewidth=0.9,
           layout="daylight",
           branch.length="none",
           aes(color=genus, linetype=Project),
           na.rm=TRUE
    ) %<+% ref.sam.metadata +
    #geom_point(aes(color=genus, shape=Project), alpha=1.0, na.rm=TRUE, size = 2) +
    scale_color_manual(values = c(Apis = 'red', Bombus = 'blue', Other = 'grey'), na.value = "grey") +
     scale_linetype_manual(values = c(`Current study`= "solid", NCBI = "dotted", Other = "solid"), na.value = "solid")

    plot(p + geom_text(aes(label=node), hjust=-.3) +
        geom_tiplab(aes(label=tip.name), size = 3, color = "black", offset=0.01))
    plot(p)
}

CollapseTree <- function(node_to_collapse) {
    p <<- p %>% ggtree::collapse(node=node_to_collapse)
    p + geom_text(aes(label=node), hjust=-.3) # table the nodes to find them
    plot(p)
}

LabelTreeClade <- function(NODE, LABEL) {
    p <<- p + geom_cladelab(node=NODE, label=LABEL, align=TRUE, 
                  geom='label', fill='grey', hjust = 1)
    plot(p)
}

############### TEMPLATE
#### main tree
#tree_ext <- ".contree"
#virus_name <- "test"

#output_name <- "test"
#PrepMetaDataTreeDNA()
#DrawTree1(phytree)

##### subtree
#sub_tree <- tree_subset(tree, node=813, levels_back=0)
#output_name <- "test"
#DrawTree1(sub_tree)

###################################################

tree_ext <- ".contree"

##### Apis rhabovirus
#### main tree
virus_name <- "Apis_rhabdovirus"
output_name <- "Apis_rhabdovirus"
PrepMetaDataTreeDNA()
DrawTree1(phytree)

assign(output_name, p)

##### Dicistroviridae
#### main tree
virus_name <- "Dicistroviridae"
output_name <- "Dicistroviridae"
PrepMetaDataTreeDNA()
DrawTree1(phytree)

assign(output_name, p)

#### Dicistroviridae_Aparavirus
virus_name <- "Dicistroviridae_Aparavirus"
output_name <- "Dicistroviridae_Aparavirus"
PrepMetaDataTreeDNA()
DrawTree1(phytree)
assign(output_name, p)

#### Dicistroviridae_Cripavirus
virus_name <- "Dicistroviridae_Cripavirus"
output_name <- "Dicistroviridae_Cripavirus"
PrepMetaDataTreeDNA()
DrawTree1(phytree)
assign(output_name, p)

#### Dicistroviridae_Triatovirus
virus_name <- "Dicistroviridae_Triatovirus"
output_name <- "Dicistroviridae_Triatovirus"
PrepMetaDataTreeDNA()
DrawTree1(phytree)
assign(output_name, p)


##### Iflaviridae_Iflavirus_aladeformis
virus_name <- "Iflaviridae_Iflavirus_aladeformis"
output_name <- "Iflaviridae_Iflavirus_aladeformis"
PrepMetaDataTreeDNA()
DrawTree1(phytree)
assign(output_name, p)

##### Iflaviridae_Iflavirus_Bactrocera_tryoni_iflavirus_1
virus_name <- "Iflaviridae_Iflavirus_Bactrocera_tryoni_iflavirus_1"
output_name <- "Iflaviridae_Iflavirus_Bactrocera_tryoni_iflavirus_1"
PrepMetaDataTreeDNA()
DrawTree1(phytree)
assign(output_name, p)

##### Iflaviridae_Iflavirus_sacbroodi
virus_name <- "Iflaviridae_Iflavirus_sacbroodi"
output_name <- "Iflaviridae_Iflavirus_sacbroodi"
PrepMetaDataTreeDNA()
DrawTree1(phytree)
assign(output_name, p)

##### Negevirus like
virus_name <- "Negevirus_Negevirus_like"
output_name <- "Negevirus_Negevirus_like"
PrepMetaDataTreeDNA()
DrawTree1(phytree)
assign(output_name, p)

##### Partiti like
virus_name <- "Partiti_like"
output_name <- "Partiti_like"
PrepMetaDataTreeDNA()
DrawTree1(phytree)
assign(output_name, p)

##### Phasmaviridae
virus_name <- "Phasmaviridae"
output_name <- "Phasmaviridae"
PrepMetaDataTreeDNA()
DrawTree1(phytree)
assign(output_name, p)

##### Picorna like Mayfield
virus_name <- "Picorna_like_Mayfield"
output_name <- "Picorna_like_Mayfield"
PrepMetaDataTreeDNA()
DrawTree1(phytree)
assign(output_name, p)

##### Reo like
virus_name <- "Reo_like"
output_name <- "Reo_like"
PrepMetaDataTreeDNA()
DrawTree1(phytree)
assign(output_name, p)

##### Sinaivirus
virus_name <- "Sinaivirus"
output_name <- "Sinaivirus"
PrepMetaDataTreeDNA()
DrawTree1(phytree)
assign(output_name, p)

##### Virga like
virus_name <- "Virga_like"
output_name <- "Virga_like"
PrepMetaDataTreeDNA()
DrawTree1(phytree)
assign(output_name, p)


Apis_rhabdovirus + ggtitle('Apis rhabdovirus') + 
Dicistroviridae_Triatovirus + ggtitle('Dicistroviridae\nTriatovirus') + 
Dicistroviridae_Aparavirus + ggtitle('Dicistroviridae\nAparavirus') + 
Dicistroviridae_Cripavirus + ggtitle('Dicistroviridae\nCripavirus') +
Iflaviridae_Iflavirus_aladeformis + ggtitle('Iflaviridae\nIflavirus aladeformis') + 
Iflaviridae_Iflavirus_Bactrocera_tryoni_iflavirus_1 + ggtitle('Iflaviridae\nBactrocera tryoni iflavirus 1') + 
Iflaviridae_Iflavirus_sacbroodi + ggtitle('Iflaviridae\nIflavirus sacbroodi') +
Negevirus_Negevirus_like + ggtitle('Negevirus-like') +
Partiti_like + ggtitle('Partiti-like') +
Phasmaviridae + ggtitle('Phasmaviridae') +
Picorna_like_Mayfield + ggtitle('Picorna-like Mayfield') +
Reo_like + ggtitle('Reo-like') +
Sinaivirus + ggtitle('Sinaivirus') +
Virga_like + ggtitle('Virga-like') +
plot_layout(guides = 'collect', widths= c(1,1,1,1,0.5))

ggsave(plot = last_plot(), "All_fastANI_trees.pdf", dpi=300, scale=2, units = "cm")
ggsave(plot = last_plot(), "All_fastANI_trees.png", dpi=300, scale=2, units = "cm")


########################################################
########################################################
########################## protein trees 
DrawTree2 <- function(input_tree, LAYOUT, BRANCH) {
    p <<- input_tree %>% 
    ggtree(linewidth=0.9,
           layout=LAYOUT,
           branch.length=BRANCH,
           aes(color=genus, linetype=Project),
           na.rm=TRUE
    ) %<+% ref.sam.metadata +
    scale_color_manual(values = c(Apis = 'red', Bombus = 'blue', Other = 'grey'), na.value = "grey") +
     scale_linetype_manual(values = c(`Current study`= "solid", NCBI = "dotted", Other = "solid"), na.value = "solid")

    plot(p + geom_text(aes(label=node), hjust=-.3) +
        geom_tiplab(aes(label=RepresentativeName), size = 3, color = "black", offset=0.01))
    plot(p)
}


########### prep info to remove redundant sequences from trees
blastp <- read.table("all_RdRp.blastp", header = TRUE, sep = "\t")

ref.sam.metadata <- rbind(ref.metadata, sam.metadata)
ref.sam.metadata <- lapply(ref.sam.metadata, as.character)
ref.sam.metadata <- data.frame(ref.sam.metadata)
ref.sam.metadata$contig_length <- as.numeric(ref.sam.metadata$contig_length)
ref.sam.metadata$genus <- ifelse(is.na(ref.sam.metadata$genus) | (ref.sam.metadata$genus != "Apis" & ref.sam.metadata$genus != "Bombus"), "Other", ref.sam.metadata$genus)
ref.sam.metadata$genus <- ifelse(is.na(ref.sam.metadata$genus), "Other", ref.sam.metadata$genus)
ref.sam.metadata$Project <- gsub("GenBank","NCBI",ref.sam.metadata$Project)
ref.sam.metadata$Project <- ifelse(ref.sam.metadata$Project != "NCBI", "Current study", ref.sam.metadata$Project)
ref.sam.metadata$Project <- ifelse(is.na(ref.sam.metadata$Project), "NCBI", ref.sam.metadata$Project)

blastp$from_contig <- gsub("_[0-9]+:[0-9]+-[0-9]+$", "", blastp$from)
blastp$to_contig <- gsub("_[0-9]+:[0-9]+-[0-9]+$", "", blastp$to)


blastp.ref.sam.metadata <- blastp %>%
        left_join(ref.sam.metadata, by = c('from_contig' = 'contig')) %>%
        left_join(ref.sam.metadata, by = c('to_contig' = 'contig'))

# Step 1: Filter the dataframe
blastp.ref.sam.metadata <- blastp.ref.sam.metadata %>%
  filter(pident == 100, genus.x == genus.y, Project.x == Project.y)

# Step 2: Create a new column with sorted and concatenated values of 'from' and 'to'
blastp.ref.sam.metadata <- blastp.ref.sam.metadata %>%
  mutate(combination = ifelse(from < to, paste(from, to, sep = "_"), paste(to, from, sep = "_")))

# Step 3: Remove duplicates based on the 'combination' column
blastp.ref.sam.metadata <- blastp.ref.sam.metadata %>%
  distinct(combination, .keep_all = TRUE) %>%
  select(-combination)  # Remove the temporary combination column

blastp.ref.sam.metadata$from <- gsub(":", "_", blastp.ref.sam.metadata$from)
blastp.ref.sam.metadata$to <- gsub(":", "_", blastp.ref.sam.metadata$to)

###### igraph to identify representatives of connected components
library(igraph)
graph <- blastp.ref.sam.metadata %>% select(from, to) %>% graph_from_data_frame()

# Find connected components
components <- components(graph)$membership

# Find representatives for each connected component
# Initialize data frames to store representatives and non-representatives
representatives_df <- data.frame(component = numeric(), representative = character(), stringsAsFactors = FALSE)
non_representatives_df <- data.frame(component = numeric(), non_representative = character(), stringsAsFactors = FALSE)

# Find representatives for each connected component
for (comp in unique(components)) {
  nodes <- names(which(components == comp))
  representative <- nodes[1]  # Choose the first node as representative
  non_representatives <- nodes[-1]  # Exclude the representative
  representatives_df <- rbind(representatives_df, data.frame(component = comp, representative = representative))
  non_representatives_df <- rbind(non_representatives_df, data.frame(component = comp, non_representative = non_representatives))
}

non_representatives_vector <- non_representatives_df$non_representative
representatives_vector <- representatives_df$non_representative


#### setwd("C:/Users/Dean Mckeown/Downloads/Spillover_FINAL/phylogeny/iprscansplit")
tree_ext <- ".contree"

###############
virus_name <- "Pfam__Bunyavirus_RNA_dependent_RNA_polymerase__PF04196"
output_name <- "Pfam__Bunyavirus_RNA_dependent_RNA_polymerase__PF04196"

PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

###############
virus_name <- "Pfam__Mononegavirales_RNA_dependent_RNA_polymerase__PF00946"
output_name <- "Pfam__Mononegavirales_RNA_dependent_RNA_polymerase__PF00946"

PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

###############
virus_name <- "Pfam__RNA_dependent_RNA_polymerase__PF00680"
output_name <- "Pfam__RNA_dependent_RNA_polymerase__PF00680"

PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
#p + theme(plot.margin = unit(c(-1.5,-1.5,-1.5,-1.5), "cm"))
assign(output_name, p)

###############
virus_name <- "Pfam__RNA_dependent_RNA_polymerase__PF00978"
output_name <- "Pfam__RNA_dependent_RNA_polymerase__PF00978"

PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

###############
virus_name <- "Pfam__Viral_RNA_dependent_RNA_polymerase__PF00998"
output_name <- "Pfam__Viral_RNA_dependent_RNA_polymerase__PF00998"

PrepMetaDataTreeAA()
to_drop <- c("MW676134.1_2_40-305")
phytree <- drop.tip(phytree, to_drop)
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

Pfam__Bunyavirus_RNA_dependent_RNA_polymerase__PF04196 + ggtitle('Bunyavirus RdRp PF04196') +
Pfam__Mononegavirales_RNA_dependent_RNA_polymerase__PF00946 + ggtitle('Mononegavirales RdRp PF00946') +
Pfam__RNA_dependent_RNA_polymerase__PF00680 + ggtitle('RdRp PF00680') +
Pfam__RNA_dependent_RNA_polymerase__PF00978 + ggtitle('RdRp PF00978') +
Pfam__Viral_RNA_dependent_RNA_polymerase__PF00998 + ggtitle('RdRp PF00998') +
plot_layout(guides = 'collect')

ggsave(plot = last_plot(), "All_PfamAA_trees.pdf", dpi=300, scale=2, units = "cm")
ggsave(plot = last_plot(), "All_PfamAA_trees.png", dpi=300, scale=2, units = "cm")

##########################################################

virus_name <- "Apis_rhabdovirus__Pfam__Mononegavirales_RNA_dependent_RNA_polymerase__PF00946"
output_name <- "Apis_rhabdovirus__Pfam__Mononegavirales_RNA_dependent_RNA_polymerase__PF00946"
PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

virus_name <- "Dicistroviridae_Aparavirus__Pfam__RNA_dependent_RNA_polymerase__PF00680"
output_name <- "Dicistroviridae_Aparavirus__Pfam__RNA_dependent_RNA_polymerase__PF00680"
PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

virus_name <- "Dicistroviridae_Cripavirus__Pfam__RNA_dependent_RNA_polymerase__PF00680"
output_name <- "Dicistroviridae_Cripavirus__Pfam__RNA_dependent_RNA_polymerase__PF00680"
PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

virus_name <- "Dicistroviridae_Triatovirus__Pfam__RNA_dependent_RNA_polymerase__PF00680"
output_name <- "Dicistroviridae_Triatovirus__Pfam__RNA_dependent_RNA_polymerase__PF00680"
PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

virus_name <- "Iflaviridae_Iflavirus_aladeformis__Pfam__RNA_dependent_RNA_polymerase__PF00680"
output_name <- "Iflaviridae_Iflavirus_aladeformis__Pfam__RNA_dependent_RNA_polymerase__PF00680"
PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

virus_name <- "Iflaviridae_Iflavirus_Bactrocera_tryoni_iflavirus_1__Pfam__RNA_dependent_RNA_polymerase__PF00680"
output_name <- "Iflaviridae_Iflavirus_Bactrocera_tryoni_iflavirus_1__Pfam__RNA_dependent_RNA_polymerase__PF00680"
PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

virus_name <- "Iflaviridae_Iflavirus_sacbroodi__Pfam__RNA_dependent_RNA_polymerase__PF00680"
output_name <- "Iflaviridae_Iflavirus_sacbroodi__Pfam__RNA_dependent_RNA_polymerase__PF00680"
PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

virus_name <- "Negevirus_Negevirus_like__Pfam__RNA_dependent_RNA_polymerase__PF00978"
output_name <- "Negevirus_Negevirus_like__Pfam__RNA_dependent_RNA_polymerase__PF00978"
PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

virus_name <- "Partiti_like__Pfam__RNA_dependent_RNA_polymerase__PF00680"
output_name <- "Partiti_like__Pfam__RNA_dependent_RNA_polymerase__PF00680"
PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

virus_name <- "Phasmaviridae__Pfam__Bunyavirus_RNA_dependent_RNA_polymerase__PF04196"
output_name <- "Phasmaviridae__Pfam__Bunyavirus_RNA_dependent_RNA_polymerase__PF04196"
PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

virus_name <- "Picorna_like_Mayfield__Pfam__RNA_dependent_RNA_polymerase__PF00680"
output_name <- "Picorna_like_Mayfield__Pfam__RNA_dependent_RNA_polymerase__PF00680"
PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

virus_name <- "Sinaivirus__Pfam__RNA_dependent_RNA_polymerase__PF00978"
output_name <- "Sinaivirus__Pfam__RNA_dependent_RNA_polymerase__PF00978"
PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

virus_name <- "Sinaivirus__Pfam__Viral_RNA_dependent_RNA_polymerase__PF00998"
output_name <- "Sinaivirus__Pfam__Viral_RNA_dependent_RNA_polymerase__PF00998"
PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)

virus_name <- "Virga_like__Pfam__RNA_dependent_RNA_polymerase__PF00978"
output_name <- "Virga_like__Pfam__RNA_dependent_RNA_polymerase__PF00978"
PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "daylight", "none")
assign(output_name, p)


Apis_rhabdovirus__Pfam__Mononegavirales_RNA_dependent_RNA_polymerase__PF00946 +
ggtitle('Apis rhabdovirus\nMononegavirales RdRp PF00946') +
Dicistroviridae_Aparavirus__Pfam__RNA_dependent_RNA_polymerase__PF00680 +
ggtitle('Dicistroviridae\nAparavirus RdRp PF00680') +
Dicistroviridae_Cripavirus__Pfam__RNA_dependent_RNA_polymerase__PF00680 +
ggtitle('Dicistroviridae\nCripavirus RdRp PF00680') +
Dicistroviridae_Triatovirus__Pfam__RNA_dependent_RNA_polymerase__PF00680 +
ggtitle('Dicistroviridae\nTriatovirus RdRp PF00680') +
Iflaviridae_Iflavirus_aladeformis__Pfam__RNA_dependent_RNA_polymerase__PF00680 +
ggtitle('Iflaviridae\nIflavirus aladeformis\nRdRp PF00680') +
Iflaviridae_Iflavirus_Bactrocera_tryoni_iflavirus_1__Pfam__RNA_dependent_RNA_polymerase__PF00680 +
ggtitle('Iflaviridae\nBactrocera tryoni iflavirus 1\nRdRp PF00680') +
Iflaviridae_Iflavirus_sacbroodi__Pfam__RNA_dependent_RNA_polymerase__PF00680 +
ggtitle('Iflaviridae\nIflavirus sacbroodi\nRdRp PF00680') +
Negevirus_Negevirus_like__Pfam__RNA_dependent_RNA_polymerase__PF00978 +
ggtitle('Negevirus\nRdRp PF00978') +
Partiti_like__Pfam__RNA_dependent_RNA_polymerase__PF00680 +
ggtitle('Partiti-like\nRdRp PF00680') +
Phasmaviridae__Pfam__Bunyavirus_RNA_dependent_RNA_polymerase__PF04196 +
ggtitle('Phasmaviridae\nBunyavirus RdRp PF04196') +
Picorna_like_Mayfield__Pfam__RNA_dependent_RNA_polymerase__PF00680 +
ggtitle('Picorna-like Mayfield\n RdRp PF00680') +
Sinaivirus__Pfam__RNA_dependent_RNA_polymerase__PF00978 +
ggtitle('Sinaivirus\nRdRp PF00978') +
Sinaivirus__Pfam__Viral_RNA_dependent_RNA_polymerase__PF00998 +
ggtitle('Sinaivirus\nRdRp PF00998') +
Virga_like__Pfam__RNA_dependent_RNA_polymerase__PF00978 +
ggtitle('Virga-like\nRdRp PF00978') +
plot_layout(guides = 'collect')

ggsave(plot = last_plot(), "All_PfamByTaxaAA_trees.pdf", dpi=300, scale=2, units = "cm")
ggsave(plot = last_plot(), "All_PfamByTaxaAA_trees.png", dpi=300, scale=2, units = "cm")

##################### all of the RdRps
virus_name <- "all_RdRp"
output_name <- "all_RdRp"
PrepMetaDataTreeAA()
phytree <- drop.tip(phytree, non_representatives_vector)
DrawTree2(phytree, "circular", "none")

ggsave(plot = last_plot(), "All_PfamAA_trees.pdf", dpi=300, scale=2, units = "cm")
ggsave(plot = last_plot(), "All_PfamAA_trees.png", dpi=300, scale=2, units = "cm")