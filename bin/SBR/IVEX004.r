#!/usr/bin/env Rscript

##### set working directory - CHANGE TOPATH TO YOUR finalresult/IVEX004
setwd("/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/IVEX004")

##### load ggplot
library(ggplot2)
library(grid)
library(tidyverse)
library(qpdf)
library(forcats)
library(dplyr)
library(cowplot)

####################################################
###################################################
############### heatmaps
##### EVEs, per genome, core genes
inputs <- list.files(path = ".", pattern = "^IVEX004_final_001_.*_ncvog_count_EVEs_heatmap$")

for (f in inputs) {
  data <- read.delim(f, sep = "\t", header = T)
  TITLE <- unique(data$genome)

  ### stop R reordering rows
  data$EVE <- as.character(data$EVE)
  data$EVE <- factor(data$EVE, levels=unique(data$EVE))

  ### categorise based on core gene count
  data$countgroup <- cut(data$count, breaks = c(0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 19, Inf), labels = c("1","2","3","4","5","6","7","8","9","10","11-15","16-19","20+"), include.lowest=TRUE, right=TRUE)
  data$legendncvog <- paste(data$legend,"|",data$ncvog)
  data$fullname <- paste(data$EVE,"|",data$context,"|",data$contig_size,"|",data$EVE_size)

  hm <- ggplot(data, aes(legendncvog, fullname, fill=countgroup)) + geom_tile(colour="black", size=0.2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.title = element_text(size=10, hjust = 0.5)) + scale_y_discrete(expand=c(0,0)) + labs(x="NCLDV core gene",y="EVE | context | contig size | EVE size") + facet_grid(. ~ group, space="free_x", scales="free_x", switch="y") + ggtitle("core gene group") + scale_fill_viridis_d(option="magma") + theme(axis.text=element_text(size = 7)) + ggtitle(TITLE) + geom_text(aes(label = count), color="grey", size = 2) + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), axis.ticks = element_blank())

  ggsave(hm, file=paste0("single_",f,".pdf"), height=600,width=300,units="mm",dpi=300)
}

qpdf::pdf_combine(input = list.files(full.names=TRUE,pattern="^single_IVEX004_final_001_.*.pdf$"), output = "IVEX004_final_001_all_ncvog_count_EVEs_heatmapseries.pdf")

  junk <- list.files(path = ".", pattern = "^single_IVEX004_final_001_.*$")
  file.remove(junk)
