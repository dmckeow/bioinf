#### setwd("C:/Users/Dean Mckeown/Downloads/Spillover_FINAL/phylogeny")

library(gggenes)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(googlesheets4)
library(Biostrings)



# reduce to relevant metadata columns only
FilterMD <- function(x) {
    x %>% select(contig, contig_length, sample_info1, sample_info2, sample_info3, SuperBin, RepresentativeName, contig_length, miuvig_quality, Project, Sample_metadata_code, BioRep, genus, species, collection_month, collection_year, apiary, distance, flower_genus, flower_species, num_of_bees)
}

HostFilter <- function(data){
      data %>%
      filter((genus == "Apis" | (genus == "Bombus" & species == "impatiens") | Project != "Current_study"))
}

PrepMetadata <- function() {
ref.sam.metadata <- rbind(ref.metadata, sam.metadata)
ref.sam.metadata <- lapply(ref.sam.metadata, as.character)
ref.sam.metadata <- data.frame(ref.sam.metadata)
ref.sam.metadata$contig_length <- as.numeric(ref.sam.metadata$contig_length)
ref.sam.metadata$genus <- ifelse(ref.sam.metadata$genus != "Apis" & ref.sam.metadata$genus != "Bombus", "Other", ref.sam.metadata$genus)
ref.sam.metadata$genus <- ifelse(is.na(ref.sam.metadata$genus), "Other", ref.sam.metadata$genus)
ref.sam.metadata$Project <- gsub("GenBank","NCBI",ref.sam.metadata$Project)
ref.sam.metadata$Project <- ifelse(ref.sam.metadata$Project != "NCBI", "Current_study", ref.sam.metadata$Project)
ref.sam.metadata$Project <- ifelse(is.na(ref.sam.metadata$Project), "NCBI", ref.sam.metadata$Project)
ref.sam.metadata$genus_project <- paste(ref.sam.metadata$genus, ref.sam.metadata$Project, sep = "_")

return(ref.sam.metadata)
}


# Import data
ref.metadata <- read_sheet("https://docs.google.com/spreadsheets/d/14BDmFgMYqHksvIRgdRG_jLISXi7MqqAoBVyAn96cm8g/edit#gid=1598046758")
sam.metadata <- read_sheet("https://docs.google.com/spreadsheets/d/1NhxQtQRc7nGmX0t0FwOZJAziMn3r0FshtpehCa-Lo5U/edit#gid=476945537")


pr_pr <- read.csv("all.predicted_proteins.all", header = TRUE, sep = "\t")
pr_pr <- pr_pr %>% select(molecule, ID, start, end, orientation, source)
colnames(pr_pr) <- c("molecule", "gene", "start", "end", "orientation", "source")



tr_pr_ip <- read.csv("all.translated_proteins.iprscan.all", header = TRUE, sep = "\t")


## convert protein loci to nucleotide (with the cds)
tr_pr_ip$startaa <- as.numeric(tr_pr_ip$startaa)
tr_pr_ip$start <- as.numeric(tr_pr_ip$start)
tr_pr_ip$endaa <- as.numeric(tr_pr_ip$endaa)
tr_pr_ip$end <- as.numeric(tr_pr_ip$end)

tr_pr_ip$startaa_orf <- as.numeric((tr_pr_ip$startaa*3)-2)
tr_pr_ip$endaa_orf <- as.numeric((tr_pr_ip$endaa*3))

## convert cds nucleotide loci into loci on the whole contig
tr_pr_ip$startaa_ctg <- as.numeric((tr_pr_ip$startaa_orf + tr_pr_ip$start)-1)
tr_pr_ip$endaa_ctg <- as.numeric((tr_pr_ip$endaa_orf + tr_pr_ip$start)-1)


tr_pr_ip <- tr_pr_ip %>% select(molecule, signature_desc, startaa_ctg, endaa_ctg, orientation, source)

colnames(tr_pr_ip) <- c("molecule", "gene", "start", "end", "orientation", "source")

all <- rbind(pr_pr, tr_pr_ip)


ref.metadata <- FilterMD(ref.metadata)
sam.metadata <- FilterMD(sam.metadata)

ref.sam.metadata <- PrepMetadata()

all <- all %>%
		left_join(ref.sam.metadata, by = c('molecule' = 'contig'))

all <- HostFilter(all)
all <- distinct(all)
all$start <- as.numeric(all$start)
all$end <- as.numeric(all$end)
all$contig_length <- as.numeric(all$contig_length)


########## 
# get names of sequences from fasta
aln_files <- list.files(pattern = "\\.aln$")

all_deflines <- data.frame()

for (file in aln_files) {
    fasta_seq <- readDNAStringSet(file)
    deflines <- data.frame(molecule = names(fasta_seq))
    deflines$molecule <- gsub(" .*", "", deflines$molecule)
    all_deflines <- rbind(all_deflines, deflines)
}

rownames(all_deflines) <- NULL
all_deflines <- all_deflines %>% distinct()
all_deflines$reverse <- ifelse(grepl("_R_", all_deflines$molecule), "Y", "N")
all_deflines$molecule <- gsub("_R_","", all_deflines$molecule)


# Merge df1 and df2 based on the 'contig' column
merged_df <- merge(all_deflines, all, by = "molecule", all.x = TRUE)

# Subtract 'start' from 'contig_length' to create 'start2' only for matching contigs
merged_df$start2 <- ifelse(!is.na(merged_df$start), merged_df$contig_length - merged_df$start, NA)

# Subtract 'end' from 'contig_length' to create 'end2' only for matching contigs
merged_df$end2 <- ifelse(!is.na(merged_df$end), merged_df$contig_length - merged_df$end, NA)

# Print the resulting merged data frame
print(merged_df)



VIRUS = "Iflavirus aladeformis"
LEN = 5000
input_gggenes <- all %>% filter(str_detect(RepresentativeName, VIRUS), contig_length > LEN)


dummies <- make_alignment_dummies(
  data = subset(input_gggenes, source %in% "Pfam"),
  aes(xmin = start, xmax = end, y = molecule, id = gene),
  on = "RNA_helicase"
)

ggplot(subset(input_gggenes, source %in% "Prodigal_v2.6.3"), aes(xmin = start, xmax = end, y = molecule)) +
  geom_gene_arrow() +
  geom_gene_arrow(data = subset(input_gggenes, source %in% "Pfam"), aes(xmin = start, xmax = end, y = molecule, fill = gene),
  arrowhead_width = grid::unit(0, "mm"),
  arrowhead_height = grid::unit(0, "mm"),
  arrow_body_height = grid::unit(2, "mm")) +
  scale_fill_brewer(palette = "Paired") +
  geom_blank(data = dummies)