#!/usr/bin/env Rscript

# Script use: 
## Visualising the BAM files mapped to contigs using ggplot
## BAM files must be converted to a tsv numeric file using `samtools depth`
## check the wrangle datasheets header parts for infromation for what should be contained within
## each tsv file in order for this script to work.

##library(dplyr) #
##require(RColorBrewer) #
##library(svglite) #
##library(reshape) #
##library(forcats) #
##library(plotrix) #
##library(tidyr) #
##library(viridis) #
##library(viridisLite) #

require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)
require(knitr)


theme_set(theme_classic())
set.seed(1234)

################################################
############## IMPORT DATASHEETS ###############
################################################
ml <- read.dta("https://stats.idre.ucla.edu/stat/data/hsbdemo.dta")

################################################
############## WRANGLE DATASHEETS ##############
################################################

###############################################################################
############# multi nomial logistic regression #################
###############################################################################




