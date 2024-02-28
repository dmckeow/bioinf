#!/usr/bin/env Rscript

##### set working directory - CHANGE TOPATH TO YOUR finalresult/IVEX004
#setwd("~/projects/DWV/finalresult/rtPCR") ## MSI

setwd("E:/remote_backup/UMN/DWV/finalresult/readmap_Albert") ## PC

##### arguments needed to pass from command line:
## must be run this way: Rscript script.r args[1]
## args[1] ## input file

#args <- commandArgs(trailingOnly = TRUE)

##### load packages
library(ggplot2)
library(ggh4x)
library(scales)
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
d <- read.delim("allbarcodes.read.count.summary.Rformat", sep = "\t", header = T)


######################### changing data
d$reference_genome <- gsub("_"," ", d$reference_genome)
d$Condition <- gsub("Failing","Diseased", d$Condition)
d$Condition <- gsub("Package_Sample","DWV_Free", d$Condition)
d$Condition <- gsub("Apivar.*","Apivar", d$Condition)

d$Cq <- as.numeric(d$Cq)
d$vcopy <- round(d$Cq/4)*4 ## round Cq to nearest multiple of 4

###### this block multiplies the original Cq value by an order of magnitude related to the virus copy per Cq
d$vcopy <- ifelse(d$vcopy==40, d$Cq * (10 ^ 1), d$vcopy)
d$vcopy <- ifelse(d$vcopy==36, d$Cq * (10 ^ 1), d$vcopy)
d$vcopy <- ifelse(d$vcopy==32, d$Cq * (10 ^ 2), d$vcopy)
d$vcopy <- ifelse(d$vcopy==28, d$Cq * (10 ^ 3), d$vcopy)
d$vcopy <- ifelse(d$vcopy==24, d$Cq * (10 ^ 4), d$vcopy)
d$vcopy <- ifelse(d$vcopy==20, d$Cq * (10 ^ 5), d$vcopy)
d$vcopy <- ifelse(d$vcopy==16, d$Cq * (10 ^ 6), d$vcopy)
d$vcopy <- ifelse(d$vcopy==12, d$Cq * (10 ^ 7), d$vcopy)
d$vcopy <- ifelse(d$vcopy==8, d$Cq * (10 ^ 8), d$vcopy)
d$vcopy <- ifelse(d$vcopy==4, d$Cq * (10 ^ 9), d$vcopy)
d <- d %>% mutate_at(vars("vcopy"), ~replace_na(.,0))
d$vcopy <- as.numeric(d$vcopy)

########## back calculations to convert virus copy number in PCR reaction to per bee
d$vcopyperbee <- (d$vcopy*12500)/30 ## 1.2 uL tRNA to 15 mL water (*12500) with 30 bees (/30)

f1 <- filter(d, Condition == "Apivar" | Condition == "Apivar_OX" | Condition == "Apivar_THY" | Condition == "Apivar_OX_THY" | Condition == "No_Treat" | Condition == "Low_Mites" | Condition == "High_Mites" | Condition == "Diseased_Colonies" | Condition == "DWV_Free")

f1 <- aggregate(cbind(total_reads_mapped_to_genome, vcopyperbee) ~ Condition + reference_genome + PCR + title_PCR + title_mapping + title_PCR, data = f1, FUN = mean, na.rm = T)

###### whole genome figure
bc <- ggplot(f1, aes(total_reads_mapped_to_genome+1, vcopyperbee+1)) + geom_point(aes(colour=Condition, shape=PCR), alpha = 0.75, size = 3) + scale_colour_viridis_d(option = "turbo")

##bc <- bc + scale_shape_manual(values=c(25,18,17,16,15,3))

bc <- bc + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5), axis.title=element_text(size=18), axis.text = element_text(size = 10), legend.text = element_text(size = 18), legend.title = element_text(size = 20))

bc <- bc + facet_wrap(. ~ reference_genome, nrow = 3)
bc <- bc + scale_y_log10()
bc <- bc + scale_x_log10(labels=comma)


bc <- bc + theme(panel.background = element_rect(fill = "transparent", colour="black"), plot.background = element_rect(fill = "transparent"), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), panel.grid.major = element_line(size = 0.05, linetype = 'solid', colour = "grey"))

bc <- bc + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0.5)) + ylab("Average DWV copies detected by RT-PCR") + xlab("Average reads mapped to reference genomes")

ggsave(bc, file=("Albert.Fig6.pdf"), height=200,width=400,units="mm",dpi=300)

##############################################
##############################################

f2 <- filter(d, Condition == "Apivar" | Condition == "Apivar_OX" | Condition == "Apivar_THY" | Condition == "Apivar_OX_THY" | Condition == "No_Treat")
f2$Date <- gsub(".*1.*","BEFORE treatment", f2$ID)
f2$Date <- gsub(".*2.*","AFTER treatment", f2$Date)
f2$Date <- gsub(".*3.*","AFTER treatment", f2$Date)

f2$Colonies_code <- gsub("[0-9]","", f2$ID)

f2 <- aggregate(cbind(total_reads_mapped_to_genome, vcopyperbee) ~ Condition + Colonies_code + Date + reference_genome + title_mapping + title_PCR, data = f2, FUN = mean, na.rm = T)

f2$ID <- paste(f2$Colonies_code, f2$Date, sep="_")

###### whole genome figure
bc <- ggplot(f2, aes(ID, total_reads_mapped_to_genome+1)) + geom_bar(stat="identity")

bc <- bc + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5), axis.title=element_text(size=18), axis.text = element_text(size = 10), legend.text = element_text(size = 18), legend.title = element_text(size = 20), strip.text = element_text(face = "bold", size = 16))

bc <- bc + facet_nested(reference_genome ~ Colonies_code + Condition, scales="free", labeller = labeller(reference_genome = label_wrap_gen(8), Date = label_wrap_gen(8)))
bc <- bc + scale_y_log10(labels=comma)

bc <- bc + theme(panel.background = element_rect(fill = "transparent", colour="black"), plot.background = element_rect(fill = "transparent"), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), panel.grid.major = element_line(size = 0.05, linetype = 'solid', colour = "grey"))

bc <- bc + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0.5)) + xlab("Colony ID, by experimental condition") + ylab("Average reads mapped to reference genomes")

ggsave(bc, file=("Albert.Fig7.pdf"), height=300,width=400,units="mm",dpi=300)

##############################################
##############################################


############SC
sp <- ggplot(d, aes(total_reads_mapped_to_genome, vcopyperbee)) + geom_jitter(aes(colour=reference_genome, shape=PCR_genotype), position = position_jitter(width = 0.1, height = 0.1), alpha = 0.8, size = 3) + scale_colour_viridis_d(option = "turbo")

sp <- sp + theme(axis.text.x = element_text(angle = 0, hjust = 0.95, vjust = 0.5), axis.title=element_text(size=20), axis.text = element_text(size = 10), legend.text = element_text(size = 18), legend.title = element_text(size = 20))

sp <- sp + facet_nested(title_PCR + PCR ~ .)

sp <- sp + theme(panel.background = element_rect(fill = "transparent", colour="black"), plot.background = element_rect(fill = "transparent"), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), panel.grid.major = element_line(size = 0.05, linetype = 'solid', colour = "grey"))

ggsave(sp, file=("Albert.allbarcodes.read.count.summary.wholegenome.scatterplot.pdf"), height=350,width=350,units="mm",dpi=300)

####
sp <- ggplot(d, aes(total_reads_mapped_to_genomelog, vcopyperbeelog)) + geom_jitter(aes(colour=reference_genome, shape=PCR_genotype), position = position_jitter(width = 0.1, height = 0.1), alpha = 0.8, size = 3) + scale_colour_viridis_d(option = "turbo")

sp <- sp + theme(axis.text.x = element_text(angle = 0, hjust = 0.95, vjust = 0.5), axis.title=element_text(size=20), axis.text = element_text(size = 10), legend.text = element_text(size = 18), legend.title = element_text(size = 20))

sp <- sp + facet_nested(title_PCR + PCR ~ .)

sp <- sp + theme(panel.background = element_rect(fill = "transparent", colour="black"), plot.background = element_rect(fill = "transparent"), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), panel.grid.major = element_line(size = 0.05, linetype = 'solid', colour = "grey"))

ggsave(sp, file=("Albert.allbarcodes.read.count.summary.wholegenome.log.scatterplot.pdf"), height=350,width=350,units="mm",dpi=300)

####
sp <- ggplot(d, aes(total_reads_mapped_to_specific_region, vcopyperbee)) + geom_jitter(aes(colour=reference_genome, shape=PCR_genotype), position = position_jitter(width = 0.1, height = 0.1), alpha = 0.8, size = 3) + scale_colour_viridis_d(option = "turbo")

sp <- sp + theme(axis.text.x = element_text(angle = 0, hjust = 0.95, vjust = 0.5), axis.title=element_text(size=20), axis.text = element_text(size = 10), legend.text = element_text(size = 18), legend.title = element_text(size = 20))

sp <- sp + facet_nested(title_PCR + PCR ~ .)

sp <- sp + theme(panel.background = element_rect(fill = "transparent", colour="black"), plot.background = element_rect(fill = "transparent"), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), panel.grid.major = element_line(size = 0.05, linetype = 'solid', colour = "grey"))

ggsave(sp, file=("Albert.allbarcodes.read.count.summary.RDRP.scatterplot.pdf"), height=350,width=350,units="mm",dpi=300)
