#!/usr/bin/env Rscript

##### set working directory - CHANGE TOPATH TO YOUR finalresult/IVEX_stramenopiles
setwd("/shared/projects/phaeoexplorer_virus/phaeoex_screen/finalresult/IVEX_stramenopiles")

##### load ggplot
library(ggplot2)
library(grid)
library(tidyverse)
library(qpdf)
library(forcats)
library(dplyr)
library(cowplot)
library(ggh4x)

####################################################

###################################################
############### heatmaps
##### all genomes, core genes
data <- read.delim("IVEX_stramenopiles.blast.reduced.quicksummary", sep = "\t", header = T)

### stop R reordering rows
data$coregene <- as.character(data$coregene)
data$coregene <- factor(data$coregene, levels=unique(data$coregene))

## set order of factor for grid facet order
data$evalue <- factor(data$evalue, levels=c("<=e-05","<=e-20","<=e-40"))


############# all genomes figure version ##########
### showing core gene count per genome

### categorise based on core gene count
####data$countgroup <- cut(data$count, breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 19, Inf), labels = c("1","2","3","4","5","6","7","8","9","10","11-15","16-19","20+"), include.lowest=TRUE, right=TRUE)


data$countgroup <- cut(data$count, breaks = c(0, 1, 2, 5, 10, 15, 19, Inf), labels = c("0","1","2-5","6-10","11-15","16-19","20+"), include.lowest=TRUE, right=TRUE)


#####data$taxaname <- paste(data$taxa3,"|",data$genome)

hm1 <- ggplot(data, aes(coregene, genome, fill=countgroup)) + geom_tile(colour="black", size=0.2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.title = element_text(size=10, hjust = 0.5)) + scale_y_discrete(expand=c(0,0)) + labs(x="NCLDV core gene",y="Stramenopile taxa", fill = "core gene count per genome") + facet_nested(taxa1 + taxa3 ~ evalue, space="free", scales="free", switch="y", strip = strip_vanilla(clip = "off")) + ggtitle("evalue category") + scale_fill_viridis_d(option="magma", na.value="grey") + theme(axis.text=element_text(size = 7), strip.text.y = element_text(angle = 0), axis.text.y=element_blank())
hm1 <- hm1 + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_line(colour = "black"), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), axis.ticks = element_blank())

ggsave(hm1, filename="IVEX_stramenopiles.blast.reduced.quicksummary.heatmap.allgenomes.pdf",height=400,width=300,units="mm",dpi=300)

############# summary per taxa1-3 figure version #############
### categorise based on core gene count
## set count category to presence/absence (R will add up all genomes per group together, so here we set the max to 1)
####data$countgroup <- cut(data$count, breaks = c(0, 1, 2, 5, 10, 15, 19, Inf), labels = c("0","1","2-5","6-10","11-15","16-19","20+"), include.lowest=TRUE, right=TRUE)

data <- read.delim("IVEX_stramenopiles.blast.reduced.quicksummary4", sep = "\t", header = T)

#### aggregate values together to value per taxa group of 1 and 3
######data$agggregatename <- paste(data$coregene,"|", data$evalue,"|", data$taxa1,"|", data$taxa2,"|", data$taxa3, "|", data$grouptotal_taxa1_taxa3)
######data <- aggregate(percent_of_grouptotal_taxa1_taxa3 ~ agggregatename, data, sum)
######data <- gsub('|', '\t', data)
### stop R reordering rows
data$coregene <- as.character(data$coregene)
data$coregene <- factor(data$coregene, levels=unique(data$coregene))

## set order of factor for grid facet order
data$evalue <- factor(data$evalue, levels=c("<=e-05","<=e-20","<=e-40"))

######data$percent_of_grouptotal_taxa1_taxa3 <- cut(data$percent_of_grouptotal_taxa1_taxa3, breaks = c(0, 20, 40, 60, 80, 100, Inf), labels = c("0","20","40","151-200","201-250","100+"), include.lowest=TRUE, right=TRUE)

data$taxa3_N <- paste(data$taxa3,"(", data$grouptotal_taxa1_taxa3,")")

hm2 <- ggplot(data, aes(coregene, taxa3_N, fill=percent_of_grouptotal_taxa1_taxa3)) + geom_tile(colour="black", size=0.2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.title = element_text(size=10, hjust = 0.5)) + scale_y_discrete(expand=c(0,0)) + labs(x="NCLDV core gene",y="Stramenopile taxa", fill = "% of genomes per group with at least 1 copy of core gene") + facet_nested(taxa1 ~ evalue, space="free", scales="free", switch="y", strip = strip_vanilla(clip = "off")) + ggtitle("evalue category") + scale_fill_viridis_c(option="magma", na.value="grey") + theme(axis.text=element_text(size = 7), strip.text.y = element_text(angle = 0))
hm2 <- hm2 + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_line(colour = "black"), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), axis.ticks = element_blank())

ggsave(hm2, filename="IVEX_stramenopiles.blast.reduced.quicksummary.heatmap.totalgenomespergroup.pdf",height=150,width=200,units="mm",dpi=300)
