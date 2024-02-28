#!/usr/bin/env Rscript

##### set working directory - CHANGE TOPATH TO YOUR finalresult/IVEX004
setwd("~/projects/DWV/finalresult/rtPCR")

##### arguments needed to pass from command line:
## must be run this way: Rscript script.r args[1]
## args[1] ## input file

args <- commandArgs(trailingOnly = TRUE)

##### load packages
library(ggplot2)
library(grid)
library(qpdf)
library(forcats)
library(cowplot)
library(sinaplot)
library(ggforce)
library(tidyverse); library(dplyr)

### using both Tidyverse and dplyr may cause conlift between plyr and dplyr and changing laod order does not work, so put dplyr:: before the following functions: arrange, count, desc, failwith, id, mutate, rename, summarise, summarize

####################################################
###################################################
#### MANUAL CHANGES TO SCRIPT NEEDED DUE TO VARIATION IN INPUT = ## !!!

CAPTION = str_wrap("Figure X. (A) Scatter plot representing the rt-PCR Cq value of DWV genotypes per honeybee sample. Circle data points represent individual samples (30 bees). (B) Barplot representing percentage occurrences of DWV genotypes across honeybee samples, detected by rt-PCR. (C) Violin plot representing the rt-PCR Cq value distribution of DWV genotypes across honeybee samples. Circle data points represent individual samples (30 bees). + = mean Cq and x = median Cq.", 60) ## !!! caption text
LEGENDTITLE = "Genotype"

####### make scatterplot of Cq values by sample group and genotype

d <- read.delim(args[1], sep = "\t", header = T)
#d$n <- 1
#d$n <- as.numeric(d$n)
#d$Cq <- as.numeric(d$Cq)

#d <- d[!is.na(d$Cq),]
#d <- d %>% mutate(n = ifelse(Cq == 'NA',NA,n))

#d <- group_by(d, sample_group_name) %>% mutate(N = sum(n, na.rm = TRUE))

#d$sample_group_name_N <- paste(d$sample_group_name,"(",d$N,")")

sp <- ggplot(d, aes(sample_group_name, Cq)) + geom_point(aes(colour=Genotype), alpha = 0.75) + scale_colour_viridis_d() + theme(axis.title.x=element_text(size=12), axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5)) + ylim(0, 40)

sp <- sp + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0.5)) + labs(fill = LEGENDTITLE) + xlab("Sample group") + ylab("Cq")

sp <- sp + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"))


####### make bar chart of percentage positive for DWV and broken down by genotypes

d <- read.delim(args[1], sep = "\t", header = T)

d$n <- 1
d$n <- as.numeric(d$n)
d <- group_by(d, sample_group_name) %>% mutate(N = sum(n))
d <- aggregate(n ~ sample_group_name + Genotype + N, data = d, FUN = sum, na.rm = T)
d$percent <- (d$n/d$N)*100

d$sample_group_name_N <- paste(d$sample_group_name,"(",d$N,")")
d$Genotype <- gsub("^$", NA, d$Genotype)

#### replace percent with NA if Genotype is NA
d <- d %>% mutate(percent = ifelse(Genotype == 'NA',NA,percent))

###### make the barplot

bc <- ggplot(d, aes(sample_group_name_N, percent, fill=Genotype)) + geom_bar(stat="identity") + scale_fill_viridis_d() + theme(axis.title.x=element_text(size=12), axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5)) + ylim(0, 100)

bc <- bc + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0.5)) + labs(fill = LEGENDTITLE) + xlab("Sample group (sample size, including all bees sampled)") + ylab("Percentage")

bc <- bc + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"))



########### violin plot with mean and median #################
  d <- read.delim(args[1], sep = "\t", header = T)
  #d$n <- 1
  #d$n <- as.numeric(d$n)
  #d$Cq <- as.numeric(d$Cq)
  #d <- d %>% mutate(n = ifelse(Cq == 'NA',NA,n))
  #d$Genotype <- gsub("^$", NA, d$Genotype)

  #d <- group_by(d, sample_group_name) %>% mutate(N = sum(n, na.rm = TRUE))

#d$sample_group_name_N <- paste(d$sample_group_name,"(",d$N,")")

vp <- ggplot(d, aes(sample_group_name, Cq)) + geom_violin() +
  theme(axis.title.x=element_text(size=12), axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5)) + ylim(0, 40) +
  stat_summary(fun=mean, colour="black", geom="point", shape=3, size=2) +
  stat_summary(fun=median, colour="black", geom="point", shape=4, size=2) +
  geom_sina(aes(colour=Genotype), alpha = 0.75, position = "dodge", maxwidth=0.75) + scale_color_viridis_d()

vp <- vp + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0.5)) + labs(fill = LEGENDTITLE) + xlab("Sample group") + ylab("Cq")

vp <- vp + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"))

#### align plots togther before saving as final image
a1 <- plot_grid(sp,bc,vp + labs(caption = CAPTION), align = "v", axis = "bt", ncol = 1, labels = "AUTO")
ggsave(a1, file=paste0(args[1],"_sp_bc_vp.pdf"), height=450,width=300,units="mm",dpi=300) ## !!! image dimensions
