#!/usr/bin/env Rscript

##### set working directory - CHANGE TOPATH TO YOUR finalresult/IVEX004
#setwd("~/projects/DWV/finalresult/rtPCR") ## MSI

setwd("E:/remote_backup/UMN/DWV/finalresult/rtPCR") ## PC

##### arguments needed to pass from command line:
## must be run this way: Rscript script.r args[1]
## args[1] ## input file

#args <- commandArgs(trailingOnly = TRUE)

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

############ for _sp_vcopy_time_condition.pdf ############################################
COLOURLEGENDTITLE1 = "Queen lineage"

CAPTION1 = str_wrap("Figure X. Scatter plots representing the number of DWV genome copies per bee, as detected by rt-PCR for (A) Deformed Wing Virus A and (B) Deformed Wing Virus B. The key shows the name of the genetic lineage of the queen (Carniolan/Italian/Unknown/Varroa resistant), and whether or not the queen had been replaced by a queen of a different lineage by the end of the experiment. Each data point represents a single colony.", 80)

f <- file("KLee_ZBQ_full_rtpcr1_sp_vcopy_lineage.txt")
writeLines(CAPTION1, f)
close(f)

X_LAB1a = "Yard"
Y_LAB1a = "Viral genome copies per bee"

H_FIG1 = 600
W_FIG1 = 300

############ for _bc ############################################

COLOURLEGENDTITLE2 = "Queen lineage"

CAPTION2 = str_wrap("Figure X. Barplot representing the percentage of colonies that were positive for DWV, as detected by rt-PCR. Each bar portion represents the DWV positive percentage per location and queen lineage. Legend indicates the genetic lineage of the colony queen (Carniolan/Italian/Unknown/Varroa resistant) and (changed) indicates that the queen had been replaced by date 4.", 80)

X_LAB2 = "Colony location, yard (total number of colonies sampled)"
Y_LAB2 = "% of colonies DWV positive"
H_FIG2 = 200
W_FIG2 = 200

############ for _vp ############################################

COLOURLEGENDTITLE3 = "Queen lineage"

CAPTION3 = str_wrap("Figure X. Violin plot representing the number of DWV genome copies per bee, as detected by rt-PCR. Legend indicates the genetic lineage of the colony queen (Carniolan/Italian/Unknown/Varroa resistant) and (changed) indicates that the queen had been replaced by date 4. Legend also indicates queen status. Black + = mean and black x = median.", 80)

X_LAB3 = "Colony location, yard"
Y_LAB3 = "Viral genome copies per bee (Log10+1)"
H_FIG3 = 300
W_FIG3 = 500

#########################################
####### make scatterplot of Cq values by Date, Condition, and Colony, also with genotype
d <- read.delim("KLee_ABC_final_rtpcr1", sep = "\t", header = T)

###### some data changes
d$original_queen_by_experiment_end <- gsub("NO_DATA","unknown", d$og_queen)
d$queen_lineage <- gsub("NO_DATA","unknown", d$queen_lineage)
d$yard <- gsub("_"," ", d$yard)

d$Cq <- ifelse(d$DWV_A_Genotype=="NO_DATA", gsub(".*",NA, d$DWV_A_Cq), d$DWV_A_Cq) ## change Cq to NA if NO genotype


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
d$vcopyperbeelog <- log10(d$vcopyperbee + 1) ## log transform

##### group and aggregate by chosen variable(s) and calc % DWV positive within them
###### make binary info for number of samples (n) and positive or negative (npos)
#### these will be used to calculate group totals (N) and totals/percentages that are postiive (Npos, percentpos)

##### save dataframe as csv file
write.csv(d,"KLee_ZBQ_full_rtpcr1_DWVA.csv", row.names = T)

sp1 <- ggplot(d, aes(yard, vcopyperbee+1)) + geom_jitter(aes(colour=queen_lineage, shape=original_queen_by_experiment_end), position = position_jitter(width = 0.2, height = 0.0), alpha = 0.9, size=6) + scale_colour_viridis_d(option="magma") + scale_y_log10()

sp1 <- sp1 + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5), axis.title=element_text(size=18), axis.text = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22))

sp1 <- sp1 + xlab(X_LAB1a) + ylab(Y_LAB1a)
sp1 <- sp1 + labs(colour = COLOURLEGENDTITLE1)

sp1 <- sp1 + theme(strip.text = element_text(face = "bold", size = 20))

sp1 <- sp1 + theme(panel.background = element_rect(fill = "grey", colour="black"), plot.background = element_rect(fill = "transparent"), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), panel.grid.major = element_line(size = 0.05, linetype = 'solid', colour = "white"))


#ggsave(sp1, file=("KLee_ZBQ_full_rtpcr1_sp_vcopy_lineage_DWVA.pdf"), height=H_FIG1,width=W_FIG1,units="mm",dpi=300)

###################################################

d$Cq <- ifelse(d$DWV_B_Genotype=="NO_DATA", gsub(".*",NA, d$DWV_B_Cq), d$DWV_B_Cq) ## change Cq to NA if NO genotype


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
d$vcopyperbeelog <- log10(d$vcopyperbee + 1) ## log transform

##### group and aggregate by chosen variable(s) and calc % DWV positive within them
###### make binary info for number of samples (n) and positive or negative (npos)
#### these will be used to calculate group totals (N) and totals/percentages that are postiive (Npos, percentpos)

##### save dataframe as csv file
write.csv(d,"KLee_ZBQ_full_rtpcr1_DWVB.csv", row.names = T)

sp2 <- ggplot(d, aes(yard, vcopyperbee+1)) + geom_jitter(aes(colour=queen_lineage, shape=original_queen_by_experiment_end), position = position_jitter(width = 0.2, height = 0.0), alpha = 0.9, size=6) + scale_colour_viridis_d(option="magma") + scale_y_log10()

sp2 <- sp2 + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5), axis.title=element_text(size=18), axis.text = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 22))

sp2 <- sp2 + xlab(X_LAB1a) + ylab(Y_LAB1a)
sp2 <- sp2 + labs(colour = COLOURLEGENDTITLE1)

sp2 <- sp2 + theme(strip.text = element_text(face = "bold", size = 20))

sp2 <- sp2 + theme(panel.background = element_rect(fill = "grey", colour="black"), plot.background = element_rect(fill = "transparent"), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), panel.grid.major = element_line(size = 0.05, linetype = 'solid', colour = "white"))


#ggsave(sp2, file=("KLee_ZBQ_full_rtpcr1_sp_vcopy_lineage_DWVB.pdf"), height=H_FIG1,width=W_FIG1,units="mm",dpi=300)


#########################
legend <- get_legend(sp1 + theme(legend.position="right"))

sp <- plot_grid(sp1 + theme(legend.position="none"),sp2 + theme(legend.position="none"), align = "h", axis = "bt", ncol = 1, labels = "AUTO", rel_heights = c(1, 1), rel_widths = c(1, 1))
sp <- plot_grid(sp, legend, align = "h", axis = "bt", ncol =  1, rel_heights = c(1, 0.15))

ggsave(sp, file=("KLee_ZBQ_full_rtpcr1_sp_vcopy_condition.pdf"), height=H_FIG1,width=W_FIG1,units="mm",dpi=300)
