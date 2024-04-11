##setwd("C:/Users/Dean Mckeown/Downloads/Spillover_FINAL/phylogeny")
library(ape)
library(Biostrings)
library(dplyr)

ProcessGffs <- function(file_common_name) {
GffProd <- read.gff(paste(file_common_name, ".gff", sep = ""))
GffIprs <- read.gff(paste(file_common_name, ".iprscan", sep = "")) ### need to remove the fasta part of the gff before this
FaaProd <- readAAStringSet(paste(file_common_name, ".faa", sep = ""))
FaaProdNames <- as.data.frame(names(FaaProd)) ## get deflines

### split the prodigal defline to get link between prodigal name and nucleotide coords
split_col <- strsplit(as.character(FaaProdNames[,1]), " # ")
split_col <- as.data.frame(do.call(rbind, split_col))
FaaProdLoc <- split_col %>% select(V1,V2,V3)
FaaProdLoc$V4 <- gsub("_[0-9]+$","",FaaProdLoc[, 1])
colnames(FaaProdLoc) <- c("protein","start","end","contig")


#### fix the ipr gff
GffIprs <- left_join(GffIprs, FaaProdLoc, by = c("seqid"="protein")) ## add nucleotide coords for proteins

## convert protein loci to nucleotide (with the cds)
GffIprs$start.x <- as.numeric(GffIprs$start.x)
GffIprs$start.y <- as.numeric(GffIprs$start.y)
GffIprs$end.x <- as.numeric(GffIprs$end.x)
GffIprs$end.y <- as.numeric(GffIprs$end.y)

GffIprs$start.cds <- as.numeric((GffIprs$start.x*3)-2)
GffIprs$end.cds <- as.numeric((GffIprs$end.x*3))

## convert cds nucleotide loci into loci on the whole contig
GffIprs$start <- as.numeric((GffIprs$start.cds + GffIprs$start.y)-1)
GffIprs$end <- as.numeric((GffIprs$end.cds + GffIprs$start.y)-1)

## move the protein id, protein loci, and cds loci into attributes
GffIprs$seqid <- gsub("^","aa_seqid=",GffIprs$seqid)
GffIprs$start.x <- gsub("^","aa_start=",GffIprs$start.x)
GffIprs$end.x <- gsub("^","aa_end=",GffIprs$end.x)
GffIprs$start.cds <- gsub("^","cds_start=",GffIprs$start.cds)
GffIprs$end.cds <- gsub("^","cds_end=",GffIprs$end.cds)
GffIprs$attributes <- paste(GffIprs$attributes, GffIprs$seqid, GffIprs$start.x, GffIprs$end.x, GffIprs$start.cds, GffIprs$end.cds, sep = ";")

GffIprs <- GffIprs %>% select(contig, source, type, start, end, score, strand, phase, attributes)
GffIprs <- GffIprs %>% rename(seqid = contig)

########### other format fixes
GffIprs$source <- as.character(GffIprs$source)
col <- "source"
GffIprs[[col]][is.na(GffIprs[[col]])] <- "InterProscan"
GffIprs$phase <- as.numeric(GffIprs$phase)
GffIprs[is.na(GffIprs)] <- 0

GffProd$strand <- gsub(".*",".", GffProd$strand)


concat <- rbind(GffIprs, GffProd)
write.table(concat, paste(file_common_name, ".gff", sep = ""), sep="\t", row.names=F, col.names=F, quote = F)
}

ProcessGffs("Apis_rhabdovirus.fa.predicted_proteins")
ProcessGffs("Phasmaviridae.fa.predicted_proteins")
ProcessGffs("Dicistroviridae.fa.predicted_proteins")
ProcessGffs("Picorna_like_Mayfield.fa.predicted_proteins")
ProcessGffs("Iflaviridae.fa.predicted_proteins")
ProcessGffs("Reo_like.fa.predicted_proteins")
ProcessGffs("Negevirus_like.fa.predicted_proteins")
ProcessGffs("Sinaivirus.fa.predicted_proteins")
ProcessGffs("Partiti_like.fa.predicted_proteins")
ProcessGffs("Virga_like.fa.predicted_proteins")