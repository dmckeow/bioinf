#!/usr/bin/env Rscript

##### set working directory - CHANGE TOPATH TO YOUR finalresult/IVEX004
#setwd("~/projects/DWV/finalresult/rtPCR") ## MSI

setwd("E:/remote_backup/UMN/DWV/finalresult/rtPCR") ## PC

##### load packages
library(ggplot2)
library(scales)
library(grid)
library(qpdf)
library(forcats)
library(cowplot)
library(sinaplot)
library(ggforce)
library(ggh4x)
library(tidyverse); library(dplyr)

### using both Tidyverse and dplyr may cause conlift between plyr and dplyr and changing laod order does not work, so put dplyr:: before the following functions: arrange, count, desc, failwith, id, mutate, rename, summarise, summarize

####################################################
###################################################
#### MANUAL CHANGES TO SCRIPT NEEDED DUE TO VARIATION IN INPUT = ## !!!

############ for _bc_vcopy_time_condition.pdf ############################################

CAPTION1 = str_wrap("Figure X. Scatterplot representing the number of DWV genome copies per bee (log10+1), as detected by rt-PCR. Each data point represents a sample of a unique set of colonies. Multiple samples of a colony set were taken across multiple time points (dates 1-3). All date 0 samples were only sampled once. All negative results are shown as a zero on the Y-axis; absent data points means no sample existed for that particular combination of variables. For every sample, 30 bees were homogenised for rt-PCR. All date 1 samples were taken before application of the condition (if applicable), whilst dates 2 and/or 3 were post application. The dates correspond for the colony codes as follows. Grouping Henrys: date 1, 06/10/2021; date 2, 06/29/2021; date 3, 08/17/2021. Grouping Kokays: date 1, 06/07/2021; date2, 08/09/2021. Grouping Martins: date 1, 05/03/2021; date 2, 06/15/2021; date 3, 08/17/2021. Key indicates whether or not Oxford Nanopore sequencing was performed on each sample. The condition acronyms are as follows: OX, oxalic acid; THY, Thymol.", 80)

f <- file("Canada_Albert_150622_DWVB_rtpcr1_bc_vcopy_time_condition.txt")
writeLines(CAPTION1, f)
close(f)

X_LAB1 = "Date sampled"
Y_LAB1 = "DWV genome copies per bee (Log10+1)"
H_FIG1 = 400
W_FIG1 = 450

#########################################
####### make scatterplot of Cq values by Date, Condition, and Colony, also with genotype
d <- read.delim("Canada_Albert_150622_DWVB_rtpcr1", sep = "\t", header = T)

#### some data changes
d$Date <- gsub("^$","Date_0", d$Date)
d$Date <- gsub("Date_1","Date_1_(pre-treatment)", d$Date)
d$Sequenced <- gsub("^Y$","Yes", d$Sequenced)
d$Sequenced <- gsub("^$","No", d$Sequenced)
d$Grouping <- gsub("_"," ", d$Grouping)
d$Grouping <- gsub(" s","s", d$Grouping)
d$Condition <- gsub("_"," ", d$Condition)
d$Condition <- gsub("Apivar Both","Apivar OX THY", d$Condition)

#### Reordering group factor levels for facets
d$Condition <- factor(d$Condition, levels = c("Apivar", "No Treat", "Apivar OX", "Apivar THY", "Apivar OX THY", "Low Mites", "High Mites", "Failing Colonies", "Package Sample"))

##### CONVERT Cq to Virus genome copies per bee
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
d$vcopyperbeelog <- log10(d$vcopyperbee+1) ## log transform

##### group and aggregate by chosen variable(s) and calc % DWV positive within them
###### make binary info for number of samples (n) and positive or negative (npos)
#### these will be used to calculate group totals (N) and totals/percentages that are postiive (Npos, percentpos)
d$n <- 1
d$npos <- (d$Cq*0)+1 ## convert to binary
d <- tidyr::replace_na(d, list(npos=0)) ## convert to binary

##### save dataframe as csv file
write.csv(d,"Canada_Albert_150622_DWVB_rtpcr1.csv", row.names = T)

sp <- ggplot(d, aes(Date, vcopyperbeelog)) + geom_jitter(aes(shape = Sequenced, colour = Lineage), position = position_jitter(width = 0.2, height = 0), size = 4, alpha = 0.75) + scale_shape_manual(values=c(1, 2)) + scale_color_viridis_d(option = "turbo")

sp <- sp + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5), axis.title=element_text(size=20), axis.text = element_text(size = 14), legend.text = element_text(size = 18), legend.title = element_text(size = 20))

sp <- sp + facet_grid(Grouping ~ Condition, switch="x", labeller = labeller(Condition = label_wrap_gen(5)))

sp <- sp + theme(strip.text = element_text(face = "bold", size = 20))

sp <- sp + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0.5)) + xlab(X_LAB1) + ylab(Y_LAB1)

sp <- sp + theme(panel.background = element_rect(fill = "transparent", colour = "black"), plot.background = element_rect(fill = "transparent"), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"))

ggsave(sp, file=("Canada_Albert_150622_DWVB_rtpcr1_bc_vcopy_time_condition.pdf"), height=H_FIG1,width=W_FIG1,units="mm",dpi=300) ## !!! image dimensions









#########################################################
####################################################################
#####FINAL REPORT FIGURES
## figure 1
d <- read.delim("Canada_Albert_300622_DWV_both_rtpcr1", sep = "\t", header = T)
d$Condition <- gsub("Failing","Diseased", d$Condition)
d$Condition <- gsub("Package_Sample","DWV_Free", d$Condition)
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

f1a <- filter(d, Condition == "Diseased_Colonies" | Condition == "DWV_Free")

bp <- ggplot(f1a, aes(Colonies_code, vcopyperbee+1)) + geom_bar(position = "dodge", stat = "summary", fun = "mean") + geom_jitter(aes(color=Bins), size=3, alpha=0.75, width=0.25, height=0.25) + scale_color_viridis_d(option = "turbo")

###bp <- bp + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))
bp <- bp + scale_y_log10()

bp <- bp + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5), axis.title=element_text(size=20), axis.text = element_text(size = 14), legend.text = element_text(size = 18), legend.title = element_text(size = 20))

bp <- bp + facet_grid(. ~ Condition, switch="x", labeller = labeller(Condition = label_wrap_gen(5)), scales="free")

bp <- bp + theme(strip.text = element_text(face = "bold", size = 20))

bp <- bp + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0.5)) + xlab("Sample ID") + ylab("DWV genome copies per bee")

bp <- bp + theme(panel.background = element_rect(fill = "transparent", colour = "black"), plot.background = element_rect(fill = "transparent"), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "transparent"))

ggsave(bp, file=("Canada_Albert_300622_Fig1a.pdf"), height=200,width=200,units="mm",dpi=300) ## !!! image dimensions


#########################

f1b <- filter(d, Condition == "Low_Mites" | Condition == "High_Mites")

bp <- ggplot(f1b, aes(Colonies_code, vcopyperbee+1)) + geom_bar(position = "dodge", stat = "summary", fun = "mean") + geom_jitter(aes(color=Bins), size=3, alpha=0.75, width=0.25, height=0.25) + scale_color_viridis_d(option = "turbo")

##bp <- bp + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))
bp <- bp + scale_y_log10()

bp <- bp + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5), axis.title=element_text(size=20), axis.text = element_text(size = 14), legend.text = element_text(size = 18), legend.title = element_text(size = 20))

bp <- bp + facet_grid(. ~ Condition, switch="x", labeller = labeller(Condition = label_wrap_gen(5)), scales="free")

bp <- bp + theme(strip.text = element_text(face = "bold", size = 20))

bp <- bp + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0.5)) + xlab("Sample ID") + ylab("DWV genome copies per bee")

bp <- bp + theme(panel.background = element_rect(fill = "transparent", colour = "black"), plot.background = element_rect(fill = "transparent"), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "transparent"))

ggsave(bp, file=("Canada_Albert_300622_Fig1b.pdf"), height=200,width=200,units="mm",dpi=300) ## !!! image dimensions



#############################################################
#############################################################

d <- read.delim("E:/remote_backup/UMN/DWV/finalresult/readmap_Albert/allbarcodes.read.count.summary.Rformat", sep = "\t", header = T)
d$Condition <- gsub("Failing","Diseased", d$Condition)
d$Condition <- gsub("Package_Sample","DWV_Free", d$Condition)

d$total_reads_mapped_to_genome <- as.numeric(d$total_reads_mapped_to_genome)
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

f2a <- filter(d, Condition == "Low_Mites" | Condition == "High_Mites" | Condition == "Diseased_Colonies" | Condition == "DWV_Free")
f2a <- filter(f2a, reference_genome == "DWV_A" | reference_genome == "DWV_B")
f2a <- aggregate(cbind(total_reads_mapped_to_genome, vcopyperbee) ~ Condition + reference_genome + PCR_genotype + PCR + title_mapping + title_PCR, data = f2a, FUN = mean, na.rm = T)

bp <- ggplot(f2a, aes(total_reads_mapped_to_genome+1, vcopyperbee+1)) + geom_point(aes(color=Condition, shape=PCR), size=2, alpha=0.75) + scale_color_viridis_d(option = "turbo")

#bp <- bp + scale_y_log10(breaks = scales::trans_breaks('log10', function(x) 10^x), labels = scales::trans_format('log10', scales::math_format(10^.x)))

#bp <- bp + scale_x_log10(breaks = scales::trans_breaks('log10', function(x) 10^x), labels = scales::trans_format('log10', scales::math_format(10^.x)))

bp <- bp + scale_y_log10()
bp <- bp + scale_x_log10(labels=comma)


bp <- bp + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5), axis.title=element_text(size=14), axis.text = element_text(size = rel(1.0)), legend.text = element_text(size = 12), legend.title = element_text(size = 14))

bp <- bp + facet_nested(. ~ title_mapping + reference_genome, switch="y", labeller = labeller(Condition = label_wrap_gen(5)))

bp <- bp + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0.5)) + ylab("Average DWV copy per bee detected by PCR") + xlab("Average number of reads mapped to reference")

bp <- bp + theme(panel.background = element_rect(fill = "transparent", colour = "black"), plot.background = element_rect(fill = "transparent"), panel.grid.minor = element_line(size = 0.1, linetype = 'solid', colour = "transparent"), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "transparent"))

ggsave(bp, file=("Canada_Albert_300622_Fig2a.pdf"), height=150,width=250,units="mm",dpi=300) ## !!! image dimensions






############# figure 3

d <- read.delim("Canada_Albert_300622_DWV_both_rtpcr1", sep = "\t", header = T)
d$Condition <- gsub("Failing","Diseased", d$Condition)
d$Condition <- gsub("Package_Sample","DWV_Free", d$Condition)
d$Condition <- gsub("Apivar_Both","Apivar_OX_THY", d$Condition)

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

f3 <- filter(d, Condition == "Apivar" | Condition == "Apivar_OX" | Condition == "Apivar_THY" | Condition == "Apivar_OX_THY" | Condition == "No_Treat")
f3$Condition <- gsub("_"," ", f3$Condition)

bp <- ggplot(f3, aes(Date, vcopyperbee+1)) + geom_bar(position = "dodge", stat = "summary", fun = "mean") + geom_jitter(aes(color=Bins, shape=Lineage), size=2, alpha=0.75, width=0.25, height=0.0) + scale_color_viridis_d(option = "turbo")

##bp <- bp + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))
bp <- bp + scale_y_log10()

bp <- bp + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5), axis.title=element_text(size=20), axis.text = element_text(size = 14), legend.text = element_text(size = 18), legend.title = element_text(size = 20))

bp <- bp + facet_grid(. ~ Condition, switch="x", labeller = labeller(Condition = label_wrap_gen(5)), scales="free")

bp <- bp + theme(strip.text = element_text(face = "bold", size = 20))

bp <- bp + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0.5)) + xlab("Treatment group") + ylab("Average DWV genome copies per bee")

bp <- bp + theme(panel.background = element_rect(fill = "transparent", colour = "black"), plot.background = element_rect(fill = "transparent"), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "transparent"))

ggsave(bp, file=("Canada_Albert_300622_Fig3.pdf"), height=200,width=250,units="mm",dpi=300) ## !!! image dimensions
