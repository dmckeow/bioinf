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
library(phangorn)
library(ggnewscale)




#### there is somethnig wrong with ggtree - this code fixes it:
nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")

# Import data
ref.metadata <- read_sheet("https://docs.google.com/spreadsheets/d/14BDmFgMYqHksvIRgdRG_jLISXi7MqqAoBVyAn96cm8g/edit#gid=1598046758")
sam.metadata <- read_sheet("https://docs.google.com/spreadsheets/d/1NhxQtQRc7nGmX0t0FwOZJAziMn3r0FshtpehCa-Lo5U/edit#gid=476945537")

# reduce to relevant metadata columns only
FilterMD <- function(x) {
    x %>% select(contig, sample_info1, sample_info2, sample_info3, SuperBin, RepresentativeName, contig_length, miuvig_quality, Project, Sample_metadata_code, BioRep, genus, species, collection_month, collection_year, apiary, distance, flower_genus, flower_species, num_of_bees)
}

HostFilter <- function(data) {
      data %>%
      filter(genus == "Bombus" & species != "impatiens" & Project == "Current_study")
}


TaxaNameCollapse1 <- function() {
    ref.sam.metadata$RepresentativeName <- gsub("Apis rhabdovirus .*", "Apis rhabdovirus", ref.sam.metadata$RepresentativeName)
    ref.sam.metadata$RepresentativeName <- gsub(".*cripavirus.*|.*Cripavirus.*", "Cripavirus", ref.sam.metadata$RepresentativeName)
    ref.sam.metadata$RepresentativeName <- gsub(".*partiti-like.*", "Partiti-like virus", ref.sam.metadata$RepresentativeName)
    ref.sam.metadata$RepresentativeName <- gsub("Allermuir.*|.*virga-like.*", "Virga-like virus", ref.sam.metadata$RepresentativeName)
    ref.sam.metadata$RepresentativeName <- gsub("Bactrocera tryoni.*", "Bactrocera tryoni iflavirus", ref.sam.metadata$RepresentativeName)
    ref.sam.metadata$RepresentativeName <- gsub(".*phasma.*", "Phasmaviridae", ref.sam.metadata$RepresentativeName)
    ref.sam.metadata$RepresentativeName <- gsub("Aparavirus.*", "Aparavirus", ref.sam.metadata$RepresentativeName)
    ref.sam.metadata$RepresentativeName <- gsub(".*Mayfield.*", "Picorna-like virus", ref.sam.metadata$RepresentativeName)

    ref.sam.metadata <- ref.sam.metadata %>%
    mutate(RepresentativeName = ifelse(RepresentativeName %in% c("Triatovirus himetobi", "Triatovirus hocoagulatae", "Triatovirus plastali", "Triatovirus triatomae", "Triatovirus nigereginacellulae"), "Triatovirus", RepresentativeName))
    ref.sam.metadata$RepresentativeName <- gsub("Lake Sinai .irus.*|.*sinaivirus.*", "Lake Sinai virus", ref.sam.metadata$RepresentativeName)

    ref.sam.metadata$RepresentativeName <- gsub(".*nege-like.*|Wallerfield.*", "Nege-like virus", ref.sam.metadata$RepresentativeName)
    return(ref.sam.metadata)
}

TaxaNameCollapse2 <- function() {

    ref.sam.metadata <- ref.sam.metadata %>%
    mutate(RepresentativeName = ifelse(RepresentativeName %in% c("Triatovirus himetobi", "Triatovirus hocoagulatae", "Triatovirus plastali", "Triatovirus triatomae"), "Other triatoviruses", RepresentativeName))
    ref.sam.metadata <- ref.sam.metadata %>%
    mutate(RepresentativeName = ifelse(RepresentativeName %in% c("Aparavirus cancerluti", "Aparavirus kashmirense", "Aparavirus tauraense", "Aparavirus vallesi"), "Other aparaviruses", RepresentativeName))

    ref.sam.metadata$RepresentativeName <- gsub("Bactrocera dorsalis cripavirus.*", "Bactrocera dorsalis cripavirus", ref.sam.metadata$RepresentativeName)
    
    return(ref.sam.metadata)
}

ref.sam.metadata <- ref.sam.metadata %>%
    mutate(RepresentativeName = ifelse(RepresentativeName %in% c("Triatovirus himetobi", "Triatovirus hocoagulatae", "Triatovirus plastali", "Triatovirus triatomae", "Triatovirus nigereginacellulae"), "Triatovirus", RepresentativeName))
ref.sam.metadata$RepresentativeName <- gsub("Lake Sinai .irus.*|.*sinaivirus.*", "Lake Sinai virus", ref.sam.metadata$RepresentativeName)

ref.sam.metadata$RepresentativeName <- gsub(".*nege-like.*|Wallerfield.*", "Nege-like virus", ref.sam.metadata$RepresentativeName)

ref.metadata <- FilterMD(ref.metadata)
sam.metadata <- FilterMD(sam.metadata)

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
ref.sam.metadata$Project <<- ifelse(ref.sam.metadata$Project != "NCBI", "Current_study", ref.sam.metadata$Project)
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
    ref.sam.metadata$genus_project <<- paste(ref.sam.metadata$genus, ref.sam.metadata$Project, sep = "_")


# Read in tree
    phytree <<- read.tree(paste0(virus_name, tree_ext))
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
ref.sam.metadata$Project <<- ifelse(ref.sam.metadata$Project != "NCBI", "Current_study", ref.sam.metadata$Project)
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

    ref.sam.metadata$genus_project <<- paste(ref.sam.metadata$genus, ref.sam.metadata$Project, sep = "_")

ref.sam.metadata$tip.name.prot <<- ref.sam.metadata$prot

# Read in tree
    phytree <<- read.tree(paste0(virus_name, tree_ext))

}


DrawTree <- function(input_tree, LAYOUT, BRANCH) {
    p <- input_tree %>% 
    ggtree(linewidth=0.9,
           layout=LAYOUT,
           branch.length=BRANCH,
           na.rm=TRUE,
           open.angle = 180
    ) %<+% ref.sam.metadata +
    
scale_color_manual(values = c(
    Apis_Current_study = '#BB0D0E',
    Apis_NCBI = '#FB9A99',
    Other_NCBI = "black",
    Bombus_Current_study = '#1F78B4',
    Bombus_NCBI = '#A6CEE3')) +

     geom_tippoint(aes(color = genus_project), size = 3) +
     geom_treescale() +
     scale_shape_manual(values = c(`Current_study`= 16, NCBI = 10, Other = 10), na.value = 16) +
     geom_text2(aes(label=label,  # add the bootstrap values
            subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),
            nudge_x = -0.01, nudge_y = 0.01, 
            check_overlap = TRUE, size = 3, color = "black")
     
     return(p)

    #plot(p + geom_text(aes(label=node), hjust=-.3) +
    #   geom_tiplab(aes(label=tip.name.contig), size = 3, color = "black", offset=0.01))
    #plot(p)
}

CollapseTree <- function(node_to_collapse) {
    p <<- p %>% ggtree::collapse(node=node_to_collapse)
    p + geom_text(aes(label=node), hjust=-.3) # table the nodes to find them
    plot(p)
}

LabelTreeClade <- function(NODE, LABEL) {
    p <- p + geom_cladelab(node=NODE, label=LABEL, align=TRUE, 
                  geom='label', fill='white', hjust = 1)
    plot(p)
    return(p)
}


TreeHeatmap <- function(INPUT, VARS_TO_SELECT, ROW_NAMES_COLUMN, WIDTH, OFFSET){
    # Reorder the dataframe based on the specified column for row names
    ref.sam.metadata.HM <- ref.sam.metadata[match(phytree$tip.label, ref.sam.metadata[, ROW_NAMES_COLUMN]), ]

    # Subset the dataframe based on the selected variables
    ref.sam.metadata.HM <- ref.sam.metadata.HM %>%
        select(all_of(VARS_TO_SELECT))
    
    # Set the row names
    rownames(ref.sam.metadata.HM) <- phytree$tip.label
    
    # Generate the heatmap
    p <- gheatmap(INPUT, ref.sam.metadata.HM, width = WIDTH, colnames = FALSE, offset = OFFSET)
    
    # Plot the heatmap
    plot(p)
    return(list(
        ref.sam.metadata.HM = ref.sam.metadata.HM,
        p = p
    ))
}


###################################################

tree_ext <- ".contree"

##### Apis rhabovirus
virus_name <- "Apis_rhabdovirus"
output_name <- "Apis_rhabdovirus"
PrepMetaDataTreeDNA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(contigAln.x) %>% as.character()
phytree <- drop.tip(phytree, to_drop)

phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "rectangular", "branch.length")

output <- TreeHeatmap(p, c("RepresentativeName"), "contigAln.x", 0.05, 0.05)
output$p <- output$p + scale_fill_brewer(palette = "Set3")

assign(output_name, output$p)

##### Dicistroviridae
#### main tree
virus_name <- "Dicistroviridae"
output_name <- "Dicistroviridae"
PrepMetaDataTreeDNA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(contigAln.x) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
phytree <- phangorn::midpoint(phytree)

p <- DrawTree(phytree, "rectangular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "contigAln.x", 0.05, 0.05)
output$p <- output$p + scale_fill_brewer(palette = "Set3")

assign(output_name, output$p)

#### Dicistroviridae_Aparavirus
virus_name <- "Dicistroviridae_Aparavirus"
output_name <- "Dicistroviridae_Aparavirus"
PrepMetaDataTreeDNA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(contigAln.x) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
to_drop <- c("_R_29JAN24DM1___barcode66.contigs.fasta______tig00000096")
phytree <- drop.tip(phytree, to_drop)
phytree <- phangorn::midpoint(phytree)

p <- DrawTree(phytree, "rectangular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "contigAln.x", 0.05, 0.05)
output$p <- output$p + scale_fill_brewer(palette = "Set3")

assign(output_name, output$p)

#### Dicistroviridae_Cripavirus
virus_name <- "Dicistroviridae_Cripavirus"
output_name <- "Dicistroviridae_Cripavirus"
PrepMetaDataTreeDNA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(contigAln.x) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "rectangular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "contigAln.x", 0.05, 0.05)
output$p <- output$p + scale_fill_brewer(palette = "Set3")

assign(output_name, output$p)

#### Dicistroviridae_Triatovirus
virus_name <- "Dicistroviridae_Triatovirus"
output_name <- "Dicistroviridae_Triatovirus"
PrepMetaDataTreeDNA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(contigAln.x) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "rectangular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "contigAln.x", 0.05, 0.05)
output$p <- output$p + scale_fill_brewer(palette = "Set3")

assign(output_name, output$p)


##### Iflaviridae_Iflavirus_aladeformis
virus_name <- "Iflaviridae_Iflavirus_aladeformis"
output_name <- "Iflaviridae_Iflavirus_aladeformis"
PrepMetaDataTreeDNA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(contigAln.x) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "rectangular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "contigAln.x", 0.05, 0.05)
output$p <- output$p + scale_fill_brewer(palette = "Set3")

assign(output_name, output$p)

##### Iflaviridae_Iflavirus_Bactrocera_tryoni_iflavirus_1
virus_name <- "Iflaviridae_Iflavirus_Bactrocera_tryoni_iflavirus_1"
output_name <- "Iflaviridae_Iflavirus_Bactrocera_tryoni_iflavirus_1"
PrepMetaDataTreeDNA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(contigAln.x) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "rectangular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "contigAln.x", 0.05, 0.05)
output$p <- output$p + scale_fill_brewer(palette = "Set3")

assign(output_name, output$p)

##### Iflaviridae_Iflavirus_sacbroodi
virus_name <- "Iflaviridae_Iflavirus_sacbroodi"
output_name <- "Iflaviridae_Iflavirus_sacbroodi"
PrepMetaDataTreeDNA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(contigAln.x) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
to_drop <- c("26JAN24DM1___barcode54.contigs.fasta______tig00000104")
phytree <- drop.tip(phytree, to_drop)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "rectangular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "contigAln.x", 0.05, 0.05)
output$p <- output$p + scale_fill_brewer(palette = "Set3")

assign(output_name, output$p)

##### Negevirus like
virus_name <- "Negevirus_Negevirus_like"
output_name <- "Negevirus_Negevirus_like"
PrepMetaDataTreeDNA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(contigAln.x) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "rectangular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "contigAln.x", 0.05, 0.05)
output$p <- output$p + scale_fill_brewer(palette = "Set3")

assign(output_name, output$p)

##### Partiti like
virus_name <- "Partiti_like"
output_name <- "Partiti_like"
PrepMetaDataTreeDNA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(contigAln.x) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "rectangular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "contigAln.x", 0.05, 0.05)
output$p <- output$p + scale_fill_brewer(palette = "Set3")

assign(output_name, output$p)

##### Phasmaviridae
virus_name <- "Phasmaviridae"
output_name <- "Phasmaviridae"
PrepMetaDataTreeDNA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(contigAln.x) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "rectangular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "contigAln.x", 0.05, 0.05)
output$p <- output$p + scale_fill_brewer(palette = "Set3")

assign(output_name, output$p)

##### Picorna like Mayfield
virus_name <- "Picorna_like_Mayfield"
output_name <- "Picorna_like_Mayfield"
PrepMetaDataTreeDNA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(contigAln.x) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "rectangular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "contigAln.x", 0.05, 0.05)
output$p <- output$p + scale_fill_brewer(palette = "Set3")

assign(output_name, output$p)

##### Reo like
virus_name <- "Reo_like"
output_name <- "Reo_like"
PrepMetaDataTreeDNA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(contigAln.x) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
to_drop <- ref.sam.metadata$contigAln.x[ref.sam.metadata$RepresentativeName == "Phocid orthoreovirus 1"]
phytree <- drop.tip(phytree, to_drop)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "rectangular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "contigAln.x", 0.05, 0.05)
output$p <- output$p + scale_fill_brewer(palette = "Set3")

assign(output_name, output$p)

##### Sinaivirus
virus_name <- "Sinaivirus"
output_name <- "Sinaivirus"
PrepMetaDataTreeDNA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(contigAln.x) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "rectangular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "contigAln.x", 0.05, 0.05)
output$p <- output$p + scale_fill_brewer(palette = "Set3")

assign(output_name, output$p)

##### Virga like
virus_name <- "Virga_like"
output_name <- "Virga_like"
PrepMetaDataTreeDNA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(contigAln.x) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "rectangular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "contigAln.x", 0.05, 0.05)
output$p <- output$p + scale_fill_brewer(palette = "Set3")

assign(output_name, output$p)

pdf("FigPhyNT_whole_multi.pdf", width = 9, height = 9)

Apis_rhabdovirus + ggtitle('A\nApis rhabdovirus')
Dicistroviridae_Triatovirus + ggtitle('B\nDicistroviridae Triatovirus')
Dicistroviridae_Aparavirus + ggtitle('C\nDicistroviridae Aparavirus')
Dicistroviridae_Cripavirus + ggtitle('D\nDicistroviridae Cripavirus')
Iflaviridae_Iflavirus_aladeformis + ggtitle('E\nIflaviridae Iflavirus aladeformis')
Iflaviridae_Iflavirus_Bactrocera_tryoni_iflavirus_1 + ggtitle('F\nIflaviridae Bactrocera tryoni iflavirus 1')
Iflaviridae_Iflavirus_sacbroodi + ggtitle('G\nIflaviridae Iflavirus sacbroodi')
Negevirus_Negevirus_like + ggtitle('H\nNegevirus-like')
Partiti_like + ggtitle('I\nPartiti-like')
Phasmaviridae + ggtitle('J\nPhasmaviridae')
Picorna_like_Mayfield + ggtitle('K\nPicorna-like Mayfield')
Reo_like + ggtitle('L\nReo-like')
Sinaivirus + ggtitle('M\nSinaivirus')
Virga_like + ggtitle('N\nVirga-like')

dev.off()


########################################################
########################################################
########################## protein trees 
#### setwd("C:/Users/Dean Mckeown/Downloads/Spillover_FINAL/phylogeny/iprscansplit")

########### prep info to remove redundant sequences from trees
blastp <- read.table("all_RdRp.blastp", header = TRUE, sep = "\t")

ref.sam.metadata <- rbind(ref.metadata, sam.metadata)
ref.sam.metadata <- lapply(ref.sam.metadata, as.character)
ref.sam.metadata <- data.frame(ref.sam.metadata)
ref.sam.metadata$contig_length <- as.numeric(ref.sam.metadata$contig_length)
ref.sam.metadata$genus <- ifelse(is.na(ref.sam.metadata$genus) | (ref.sam.metadata$genus != "Apis" & ref.sam.metadata$genus != "Bombus"), "Other", ref.sam.metadata$genus)
ref.sam.metadata$genus <- ifelse(is.na(ref.sam.metadata$genus), "Other", ref.sam.metadata$genus)
ref.sam.metadata$Project <- gsub("GenBank","NCBI",ref.sam.metadata$Project)
ref.sam.metadata$Project <- ifelse(ref.sam.metadata$Project != "NCBI", "Current_study", ref.sam.metadata$Project)
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

##############################################################################
##############################################################################

tree_ext <- ".contree"

###############
virus_name <- "Pfam__Bunyavirus_RNA_dependent_RNA_polymerase__PF04196"
output_name <- "Pfam__Bunyavirus_RNA_dependent_RNA_polymerase__PF04196"

PrepMetaDataTreeAA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
##phytree <- drop.tip(phytree, non_representatives_vector)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "circular", "branch.length")

output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.05, 0.01)
output$p <- output$p + scale_fill_brewer(palette = "Paired")
assign(output_name, output$p)


###############
virus_name <- "Pfam__Mononegavirales_RNA_dependent_RNA_polymerase__PF00946"
output_name <- "Pfam__Mononegavirales_RNA_dependent_RNA_polymerase__PF00946"

PrepMetaDataTreeAA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
##phytree <- drop.tip(phytree, non_representatives_vector)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "circular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.05, 0.01)
output$p <- output$p + scale_fill_brewer(palette = "Paired")
assign(output_name, output$p)

############### too many colors
virus_name <- "Pfam__RNA_dependent_RNA_polymerase__PF00680"
output_name <- "Pfam__RNA_dependent_RNA_polymerase__PF00680"

PrepMetaDataTreeAA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
phytree <- drop.tip(phytree, non_representatives_vector)
ref.sam.metadata <- TaxaNameCollapse1()

phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "circular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.05, 0.01)
output$p <- output$p + scale_fill_brewer(palette = "Paired")
assign(output_name, output$p)

###############
virus_name <- "Pfam__RNA_dependent_RNA_polymerase__PF00978"
output_name <- "Pfam__RNA_dependent_RNA_polymerase__PF00978"

PrepMetaDataTreeAA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
##phytree <- drop.tip(phytree, non_representatives_vector)
p <- DrawTree(phytree, "circular", "branch.length")
phytree <- phangorn::midpoint(phytree)
output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.05, 0.01)
output$p <- output$p + scale_fill_brewer(palette = "Paired")
assign(output_name, output$p)

###############
virus_name <- "Pfam__Viral_RNA_dependent_RNA_polymerase__PF00998"
output_name <- "Pfam__Viral_RNA_dependent_RNA_polymerase__PF00998"

PrepMetaDataTreeAA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
to_drop <- c("MW676134.1_2_40-305")
phytree <- drop.tip(phytree, to_drop)
##phytree <- drop.tip(phytree, non_representatives_vector)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "circular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
output$p <- output$p + scale_fill_brewer(palette = "Paired")
assign(output_name, output$p)


#pdf("FigPhyAA_Pfam_multi.pdf", width = 9, height = 9)

#Pfam__Bunyavirus_RNA_dependent_RNA_polymerase__PF04196 + ggtitle('A\nBunyavirus RdRp PF04196')
#Pfam__Mononegavirales_RNA_dependent_RNA_polymerase__PF00946 + ggtitle('B\nMononegavirales RdRp PF00946')
#Pfam__RNA_dependent_RNA_polymerase__PF00680 + ggtitle('C\nRdRp PF00680')
#Pfam__RNA_dependent_RNA_polymerase__PF00978 + ggtitle('D\nRdRp PF00978')
#Pfam__Viral_RNA_dependent_RNA_polymerase__PF00998 + ggtitle('E\nRdRp PF00998')

#dev.off()

##########################################################

virus_name <- "Apis_rhabdovirus__Pfam__Mononegavirales_RNA_dependent_RNA_polymerase__PF00946"
output_name <- "Apis_rhabdovirus__Pfam__Mononegavirales_RNA_dependent_RNA_polymerase__PF00946"
PrepMetaDataTreeAA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
#phytree <- drop.tip(phytree, non_representatives_vector)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "circular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
output$p <- output$p + scale_fill_brewer(palette = "Paired")
assign(output_name, output$p)

##

virus_name <- "Dicistroviridae__Pfam__RNA_dependent_RNA_polymerase__PF00680"
output_name <- "Dicistroviridae__Pfam__RNA_dependent_RNA_polymerase__PF00680"
PrepMetaDataTreeAA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
#phytree <- drop.tip(phytree, non_representatives_vector)

to_drop <- c("OR597291.1_2_1123-1601", "NC_025219.1_1_1331-1799")
phytree <- drop.tip(phytree, to_drop)

ref.sam.metadata <- TaxaNameCollapse2()
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "circular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
output$p <- output$p + scale_fill_brewer(palette = "Paired")
assign(output_name, output$p)

##
#virus_name <- "Dicistroviridae_Aparavirus__Pfam__RNA_dependent_RNA_polymerase__PF00680"
#output_name <- "Dicistroviridae_Aparavirus__Pfam__RNA_dependent_RNA_polymerase__PF00680"
#PrepMetaDataTreeAA()
#to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
#phytree <- drop.tip(phytree, to_drop)
###phytree <- drop.tip(phytree, non_representatives_vector)
#to_drop <- c("29JAN24DM1___barcode66.contigs.fasta______tig00000096_3_16-491", "29JAN24DM1___barcode66.contigs.fasta______tig00000096_1_14-127")
#phytree <- drop.tip(phytree, to_drop)
#phytree <- phangorn::midpoint(phytree)
#p <- DrawTree(phytree, "circular", "branch.length")
#output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
#output$p <- output$p + scale_fill_brewer(palette = "Paired")
#assign(output_name, output$p)
#
###
#virus_name <- "Dicistroviridae_Cripavirus__Pfam__RNA_dependent_RNA_polymerase__PF00680"
#output_name <- "Dicistroviridae_Cripavirus__Pfam__RNA_dependent_RNA_polymerase__PF00680"
#PrepMetaDataTreeAA()
#to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
#phytree <- drop.tip(phytree, to_drop)
###phytree <- drop.tip(phytree, non_representatives_vector)
#phytree <- phangorn::midpoint(phytree)
#p <- DrawTree(phytree, "circular", "branch.length")
#output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
#output$p <- output$p + scale_fill_brewer(palette = "Paired")
#assign(output_name, output$p)
#
###
#
#virus_name <- "Dicistroviridae_Triatovirus__Pfam__RNA_dependent_RNA_polymerase__PF00680"
#output_name <- "Dicistroviridae_Triatovirus__Pfam__RNA_dependent_RNA_polymerase__PF00680"
#PrepMetaDataTreeAA()
#to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
#phytree <- drop.tip(phytree, to_drop)
###phytree <- drop.tip(phytree, non_representatives_vector)
#phytree <- phangorn::midpoint(phytree)
#p <- DrawTree(phytree, "circular", "branch.length")
#output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
#output$p <- output$p + scale_fill_brewer(palette = "Paired")
#assign(output_name, output$p)

##

virus_name <- "Iflaviridae_Iflavirus__Pfam__RNA_dependent_RNA_polymerase__PF00680"
output_name <- "Iflaviridae_Iflavirus__Pfam__RNA_dependent_RNA_polymerase__PF00680"
PrepMetaDataTreeAA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
phytree <- drop.tip(phytree, non_representatives_vector)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "circular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
output$p <- output$p + scale_fill_brewer(palette = "Paired")
assign(output_name, output$p)

##

#virus_name <- "Iflaviridae_Iflavirus_aladeformis__Pfam__RNA_dependent_RNA_polymerase__PF00680"
#output_name <- "Iflaviridae_Iflavirus_aladeformis__Pfam__RNA_dependent_RNA_polymerase__PF00680"
#PrepMetaDataTreeAA()
#to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
#phytree <- drop.tip(phytree, to_drop)
##phytree <- drop.tip(phytree, non_representatives_vector)
#phytree <- phangorn::midpoint(phytree)
#p <- DrawTree(phytree, "circular", "branch.length")
#output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
#output$p <- output$p + scale_fill_brewer(palette = "Paired")
#assign(output_name, output$p)

##

#virus_name <- "Iflaviridae_Iflavirus_Bactrocera_tryoni_iflavirus_1__Pfam__RNA_dependent_RNA_polymerase__PF00680"
#output_name <- "Iflaviridae_Iflavirus_Bactrocera_tryoni_iflavirus_1__Pfam__RNA_dependent_RNA_polymerase__PF00680"
#PrepMetaDataTreeAA()
#to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
#phytree <- drop.tip(phytree, to_drop)
###phytree <- drop.tip(phytree, non_representatives_vector)
#phytree <- phangorn::midpoint(phytree)
#p <- DrawTree(phytree, "circular", "branch.length")
#output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
#output$p <- output$p + scale_fill_brewer(palette = "Paired")
#assign(output_name, output$p)
###

#virus_name <- "Iflaviridae_Iflavirus_sacbroodi__Pfam__RNA_dependent_RNA_polymerase__PF00680"
#output_name <- "Iflaviridae_Iflavirus_sacbroodi__Pfam__RNA_dependent_RNA_polymerase__PF00680"
#PrepMetaDataTreeAA()
#to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
#phytree <- drop.tip(phytree, to_drop)
###phytree <- drop.tip(phytree, non_representatives_vector)
#phytree <- phangorn::midpoint(phytree)
#p <- DrawTree(phytree, "circular", "branch.length")
#output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
#output$p <- output$p + scale_fill_brewer(palette = "Paired")
#assign(output_name, output$p)

##

virus_name <- "Negevirus_Negevirus_like__Pfam__RNA_dependent_RNA_polymerase__PF00978"
output_name <- "Negevirus_Negevirus_like__Pfam__RNA_dependent_RNA_polymerase__PF00978"
PrepMetaDataTreeAA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
#phytree <- drop.tip(phytree, non_representatives_vector)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "circular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
output$p <- output$p + scale_fill_brewer(palette = "Paired")
assign(output_name, output$p)

##

virus_name <- "Partiti_like__Pfam__RNA_dependent_RNA_polymerase__PF00680"
output_name <- "Partiti_like__Pfam__RNA_dependent_RNA_polymerase__PF00680"
PrepMetaDataTreeAA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
#phytree <- drop.tip(phytree, non_representatives_vector)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "circular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
output$p <- output$p + scale_fill_brewer(palette = "Paired")
assign(output_name, output$p)

virus_name <- "Phasmaviridae__Pfam__Bunyavirus_RNA_dependent_RNA_polymerase__PF04196"
output_name <- "Phasmaviridae__Pfam__Bunyavirus_RNA_dependent_RNA_polymerase__PF04196"
PrepMetaDataTreeAA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
#phytree <- drop.tip(phytree, non_representatives_vector)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "circular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
output$p <- output$p + scale_fill_brewer(palette = "Paired")
assign(output_name, output$p)

virus_name <- "Picorna_like_Mayfield__Pfam__RNA_dependent_RNA_polymerase__PF00680"
output_name <- "Picorna_like_Mayfield__Pfam__RNA_dependent_RNA_polymerase__PF00680"
PrepMetaDataTreeAA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
#phytree <- drop.tip(phytree, non_representatives_vector)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "circular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
output$p <- output$p + scale_fill_brewer(palette = "Paired")
assign(output_name, output$p)

virus_name <- "Sinaivirus__Pfam__RNA_dependent_RNA_polymerase__PF00978"
output_name <- "Sinaivirus__Pfam__RNA_dependent_RNA_polymerase__PF00978"
PrepMetaDataTreeAA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
#phytree <- drop.tip(phytree, non_representatives_vector)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "circular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
output$p <- output$p + scale_fill_brewer(palette = "Paired")
assign(output_name, output$p)

virus_name <- "Sinaivirus__Pfam__Viral_RNA_dependent_RNA_polymerase__PF00998"
output_name <- "Sinaivirus__Pfam__Viral_RNA_dependent_RNA_polymerase__PF00998"
PrepMetaDataTreeAA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
#phytree <- drop.tip(phytree, non_representatives_vector)
to_drop <- c("MW676134.1_2_40-305")
phytree <- drop.tip(phytree, to_drop)

ref.sam.metadata <- ref.sam.metadata %>%
    mutate(RepresentativeName = ifelse(RepresentativeName %in% c("Lake Sinai virus", "Unclassified sinaivirus"), "Unclassified lake sinaivirus", RepresentativeName))

phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "circular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
output$p <- output$p + scale_fill_brewer(palette = "Paired")
assign(output_name, output$p)

virus_name <- "Virga_like__Pfam__RNA_dependent_RNA_polymerase__PF00978"
output_name <- "Virga_like__Pfam__RNA_dependent_RNA_polymerase__PF00978"
PrepMetaDataTreeAA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
##phytree <- drop.tip(phytree, non_representatives_vector)
phytree <- phangorn::midpoint(phytree)
p <- DrawTree(phytree, "circular", "branch.length")
output <- TreeHeatmap(p, c("RepresentativeName"), "prot", 0.1, 0.01)
output$p <- output$p + scale_fill_brewer(palette = "Paired")
assign(output_name, output$p)

pdf("FigPhyAA_taxPfam_multi.pdf", width = 12, height = 12)

Apis_rhabdovirus__Pfam__Mononegavirales_RNA_dependent_RNA_polymerase__PF00946 +
ggtitle('A   Apis rhabdovirus\nMononegavirales RdRp PF00946')
Dicistroviridae_Aparavirus__Pfam__RNA_dependent_RNA_polymerase__PF00680 +
ggtitle('B   Dicistroviridae\nAparavirus RdRp PF00680')
Dicistroviridae_Cripavirus__Pfam__RNA_dependent_RNA_polymerase__PF00680 +
ggtitle('C   Dicistroviridae\nCripavirus RdRp PF00680')
Dicistroviridae_Triatovirus__Pfam__RNA_dependent_RNA_polymerase__PF00680 +
ggtitle('D   Dicistroviridae\nTriatovirus RdRp PF00680')
Iflaviridae_Iflavirus_aladeformis__Pfam__RNA_dependent_RNA_polymerase__PF00680 +
ggtitle('E   Iflaviridae\nIflavirus aladeformis\nRdRp PF00680')
Iflaviridae_Iflavirus_Bactrocera_tryoni_iflavirus_1__Pfam__RNA_dependent_RNA_polymerase__PF00680 +
ggtitle('F   Iflaviridae\nBactrocera tryoni iflavirus 1\nRdRp PF00680')
Iflaviridae_Iflavirus_sacbroodi__Pfam__RNA_dependent_RNA_polymerase__PF00680 +
ggtitle('G   Iflaviridae\nIflavirus sacbroodi\nRdRp PF00680')
Negevirus_Negevirus_like__Pfam__RNA_dependent_RNA_polymerase__PF00978 +
ggtitle('H   Negevirus\nRdRp PF00978')
Partiti_like__Pfam__RNA_dependent_RNA_polymerase__PF00680 +
ggtitle('I   Partiti-like\nRdRp PF00680')
Phasmaviridae__Pfam__Bunyavirus_RNA_dependent_RNA_polymerase__PF04196 +
ggtitle('J   Phasmaviridae\nBunyavirus RdRp PF04196')
Picorna_like_Mayfield__Pfam__RNA_dependent_RNA_polymerase__PF00680 +
ggtitle('K   Picorna-like Mayfield\n RdRp PF00680')
Sinaivirus__Pfam__RNA_dependent_RNA_polymerase__PF00978 +
ggtitle('L   Sinaivirus\nRdRp PF00978')
Sinaivirus__Pfam__Viral_RNA_dependent_RNA_polymerase__PF00998 +
ggtitle('M   Sinaivirus\nRdRp PF00998')
Virga_like__Pfam__RNA_dependent_RNA_polymerase__PF00978 +
ggtitle('N   Virga-like\nRdRp PF00978')

dev.off()

##################### all of the RdRps
virus_name <- "all_RdRp"
output_name <- "all_RdRp"
PrepMetaDataTreeAA()
to_drop <- HostFilter(ref.sam.metadata) %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)
to_drop <- c("29JAN24DM1___barcode66.contigs.fasta______tig00000096_3_16-491",
    "29JAN24DM1___barcode66.contigs.fasta______tig00000096_1_14-127",
    "29JAN24DM1___barcode66.contigs.fasta______tig00000070_1_13-70")
phytree <- drop.tip(phytree, to_drop)
phytree <- drop.tip(phytree, non_representatives_vector)

to_drop <- ref.sam.metadata[ref.sam.metadata$RepresentativeName == "Phocid orthoreovirus 1", ] %>% select(prot) %>% as.character()
phytree <- drop.tip(phytree, to_drop)


ref.sam.metadata <- TaxaNameCollapse1()


phytree <- phangorn::midpoint(phytree)

p <- DrawTree(phytree, "circular", "branch.length")
output <- TreeHeatmap(p, c("SuperBin"), "prot", 0.05, 0.2)
output$p <- output$p + scale_fill_brewer(palette = "Set3")
output$p <- output$p + new_scale_fill()

output <- TreeHeatmap(output$p, c("RepresentativeName"), "prot", 0.05, 0.01)
output$p <- output$p + scale_fill_brewer(palette = "Paired")


assign(output_name, output$p)

ggsave(plot = all_RdRp, "FigPhyAA_AllRdRp.pdf", dpi=300, scale=2, units = "cm")
ggsave(plot = all_RdRp, "FigPhyAA_AllRdRp.png", dpi=300, scale=2, units = "cm")




# 26JAN24DM1___barcode54.contigs.fasta______tig00000104_1_168−528 = Nege-like virus DONE, filter out of NT tree
# 29JAN24DM1___barcode66.contigs.fasta______tig00000096_1_14−127 = Iflavirus sacbroodi DONE, filter out of NT tree
# reference Mayfield viruses including Mayfield virus 1-like are Picorna-like, not Partiti-like DONE


########## dN/dS

#library(seqinr)


#alignment <- seqinr::read.alignment("Pfam__RNA_dependent_RNA_polymerase__PF00680.aln", format = "fasta")

#seqinr::reverse.align(
#    "Pfam__RNA_dependent_RNA_polymerase__PF00680.concat.fa", 
#    "Pfam__RNA_dependent_RNA_polymerase__PF00680.aln",
#    input.format = 'fasta',
#    Pfam__RNA_dependent_RNA_polymerase__PF00680.aln.reverse,
#    output.format = 'fasta',
#  )

#seqinr::kaks(alignment)



########## msa

#library(ggmsa)

#ggmsa("Pfam__RNA_dependent_RNA_polymerase__PF00680.aln", start = 221, end = 280, char_width = 0.5, seq_name = T)