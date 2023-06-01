#!/usr/bin/env Rscript

##### set working directory - CHANGE TOPATH TO YOUR finalresult/IVEX004
#setwd("~/projects/DWV/finalresult/rtPCR") ## MSI

setwd("D:/remote_backup/UMN/DWV/finalresult/rtPCR") ## PC

##### arguments needed to pass from command line:
## must be run this way: Rscript script.r args[1]
## args[1] ## input file

##args <- commandArgs(trailingOnly = TRUE)

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


#########################################
####### make scatterplot of Cq values by Date, Condition, and Colony, also with genotype
d <- read.delim("rtPCR_new_rtpcr1", sep = "\t", header = T)

d$Cq <- as.numeric(d$Cq)

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

##### save dataframe as csv file
write.csv(d,"rtPCR_new_rtpcr1.csv", row.names = T)
