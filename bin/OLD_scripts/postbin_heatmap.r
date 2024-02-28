#!/usr/bin/env Rscript

# Script use: 
## Visualising the BAM files mapped to contigs using ggplot
## BAM files must be converted to a tsv numeric file using `samtools depth`
## check the wrangle datasheets header parts for infromation for what should be contained within
## each tsv file in order for this script to work.

library(dplyr) #
library(ggplot2) #
library(RColorBrewer) #
library(svglite) #
library(reshape) #
library(forcats) #
library(plotrix) #
library(tidyr) #
library(viridis) #
library(viridisLite) #


theme_set(theme_classic())
set.seed(1234)

################################################
############## IMPORT DATASHEETS ###############
################################################
df <- read.csv("contigs-refs-VOGs.tsv", header=TRUE, sep="\t")

##sc <- read.csv("contigs-refs-VOGs.sc", header=F, sep="\t")
##sc1 <- as.numeric(sc[1]/20)

################################################
############## WRANGLE DATASHEETS ##############
################################################

###############################################################################
############# RPKM Plotting with linecolour the filter length ################# 
###############################################################################

ggplot(df, aes(VOG, fct_reorder(contig, totalVOGs), fill = bitscore)) + geom_tile(color = "black") +
labs(x="Viral Orthologous Group", y="Query") + 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face="bold")) +
theme(axis.text.y = element_text(vjust = 1, hjust=1, face="bold")) +
scale_fill_viridis(option = "plasma", na.value = "gray80") +
geom_text(aes(label = alignlength))

ggsave(plot = last_plot(),paste0("figures/VOG_summary_heatmap", format(Sys.time(), "_%Y-%m-%d"), ".pdf"), dpi=600)


