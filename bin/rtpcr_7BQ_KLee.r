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
COLOURLEGENDTITLE1 = "Mites per 100 bees"
SHAPELEGENDTITLE1 = "Frame of bees"
SIZELEGENDTITLE1 = "Frame of bees"

CAPTION1 = str_wrap("Figure X. Scatter plots representing the number of DWV genome copies per bee, as detected by rt-PCR. The key shows the frame of bees, a measure of bee colony population (the shape of data points) and the mites per 100 bees, the infestation level of Varroa mites on adult bees (the colour of data points). Each data point represents a single colony. For each colony, rt-PCR was performed on samples collected on date 4. The dates correspond for the months of 2021 as follows: date 1, February; date 2, June; date 3, August; date 4, September. Plots are faceted by the genetic lineage of the colony queen (Carniolan/Italian/Unknown/Varroa resistant) and (changed) indicates that the queen had been replaced by date 4. Plots are also faceted by the location of colonies (yard).", 80)

f <- file("KLee_ZBQ_full_rtpcr1_sp_vcopy_time_condition.txt")
writeLines(CAPTION1, f)
close(f)

X_LAB1a = "Date"
Y_LAB1a = "Viral genome copies per bee (Log10+1)"
##X_LAB1b = "Mites per 1000 bees (Log10+1)"
##Y_LAB1b = "Viral genome copies per bee (Log10+1)"

H_FIG1 = 400
W_FIG1 = 700

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
d$og_queen <- gsub("^FALSE$","changed", d$og_queen) ##
d$og_queen <- gsub("^TRUE$","", d$og_queen) ##
d$qstatus <- gsub("^qs$|^vs$|^dl$|^lw$","normal", d$qstatus) ## lump "normal" queen activities together
d$qstatus <- gsub("^ql$","absent", d$qstatus) ##
d <- d %>% mutate_at(vars("qstatus"), ~replace_na(.,"unknown"))
d$queen_lineage_og <- paste(d$queen_lineage," (",d$og_queen,")", sep="")
d$yard <- gsub("_"," ", d$yard)
d$queen_lineage_og <- gsub("_"," ", d$queen_lineage_og)
d$queen_lineage_og <- gsub(" \\(\\)","", d$queen_lineage_og)
d$queen_lineage_og <- gsub("_"," ", d$queen_lineage_og)

d$Cq <- ifelse(d$DWV_A_Genotype=="", gsub(".*",NA, d$DWV_A_Cq), d$DWV_A_Cq) ## change Cq to NA if NO genotype


##### change to some variables to discrete variables
d$to_bin <- d$mites_per_100bees ## specify variable to bin
bin_no=4 ## set number of categories

#### binning:
d$to_bin_noz <- ifelse(d$to_bin==0, NA, d$to_bin) ## remove zeroes
MX = max(d$to_bin_noz, na.rm =T) ## get highest value
MN = min(d$to_bin_noz, na.rm =T) ## get lowest value (non-zero)
d <- d %>% mutate(binned = cut(d$to_bin_noz, breaks = seq(MN,MX, length.out = bin_no), na.rm = T))
d$binned <- gsub("^.|.$","", d$binned)
d$binned <- gsub(","," to ", d$binned)
d$binned <- ifelse(d$to_bin==0, 0, d$binned) ## add category for zeroes if needed
d <- d %>% mutate_at(vars("binned"), ~replace_na(.,"unknown")) ## maybe specify what NAs are categorised as
#### rename binned variable to something specific
d <- d %>% rename(mites_per_100bees_binned = binned)
###

##### change to some variables to discrete variables
d$to_bin <- d$fob ## specify variable to bin
bin_no=4 ## set number of categories

#### binning:
d$to_bin_noz <- ifelse(d$to_bin==0, NA, d$to_bin) ## remove zeroes
MX = max(d$to_bin_noz, na.rm =T) ## get highest value
MN = min(d$to_bin_noz, na.rm =T) ## get lowest value (non-zero)
d <- d %>% mutate(binned = cut(d$to_bin_noz, breaks = seq(MN,MX, length.out = bin_no), na.rm = T))
d$binned <- gsub("^.|.$","", d$binned)
d$binned <- gsub(","," to ", d$binned)
d$binned <- ifelse(d$to_bin==0, 0, d$binned) ## add category for zeroes if needed
d <- d %>% mutate_at(vars("binned"), ~replace_na(.,"unknown")) ## maybe specify what NAs are categorised as
#### rename binned variable to something specific
d <- d %>% rename(fob_binned = binned)
###

##### CONVERT Cq to Virus genome copies per bee
d$fob <- as.numeric(d$fob)
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

##### group and aggregate by chosen variable(s) and calc % DWV positive within them
###### make binary info for number of samples (n) and positive or negative (npos)
#### these will be used to calculate group totals (N) and totals/percentages that are postiive (Npos, percentpos)
d$n <- 1
d$npos <- (d$Cq*0)+1 ## convert to binary
d <- tidyr::replace_na(d, list(npos=0)) ## convert to binary

##### set order for key
d$fob_binned <- factor(d$fob_binned, levels = c("unknown", "0", "0.5 to 5.33", "5.33 to 10.2", "10.2 to 15"))
d$mites_per_100bees_binned <- factor(d$mites_per_100bees_binned, levels = c("unknown", "0", "0.243 to 1.75", "1.75 to 3.25", "3.25 to 4.75"))

##### save dataframe as csv file
write.csv(d,"KLee_ZBQ_full_rtpcr1.csv", row.names = T)

#sp1 <- ggplot(d, aes(Date, vcopyperbeelog)) + geom_jitter(aes(colour=mites_per_100bees_binned, size=fob_binned, shape=qstatus), position = position_jitter(width = 0.2), alpha = 0.9) + scale_colour_viridis_d(option="magma")
#sp1 <- sp1 +scale_shape_manual(values=c(1, 4, 16, 17))

sp1 <- ggplot(d, aes(queen_lineage_og, vcopyperbee+1)) + geom_jitter(aes(colour=mites_per_100bees_binned, shape=fob_binned), position = position_jitter(width = 0.2), alpha = 0.9, size=2) + scale_colour_viridis_d(option="magma")

sp1 <- sp1 +scale_shape_manual(values=c(4, 1, 16, 17, 15))

sp1 <- sp1 + facet_grid(yard ~ .)
##sp1 <- sp1 + facet_grid(yard ~ ., labeller = labeller(queen_lineage_og = label_wrap_gen(12)))

sp1 <- sp1 + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5), axis.title=element_text(size=18), axis.text = element_text(size = 16), legend.text = element_text(size = 16), legend.title = element_text(size = 18))

sp1 <- sp1 + scale_y_log10()

sp1 <- sp1 + xlab(X_LAB1a) + ylab(Y_LAB1a)
sp1 <- sp1 + labs(colour = COLOURLEGENDTITLE1, shape = SHAPELEGENDTITLE1, size = SIZELEGENDTITLE1)

sp1 <- sp1 + theme(strip.text = element_text(face = "bold", size = 20))

sp1 <- sp1 + theme(panel.background = element_rect(fill = "grey", colour="black"), plot.background = element_rect(fill = "transparent"), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), panel.grid.major = element_line(size = 0.05, linetype = 'solid', colour = "white"))

ggsave(sp1, file=("KLee_ZBQ_full_rtpcr1_sp_vcopy_time_condition.pdf"), height=H_FIG1,width=W_FIG1,units="mm",dpi=300)

###################################################

sp2 <- ggplot(d, aes(mites_per_100bees, vcopyperbee)) + geom_jitter(aes(colour=queen_lineage_og, shape=qstatus), position = position_jitter(width = 0.1, height = 0.1), alpha = 0.7, size = 2) + scale_colour_viridis_d(option = "turbo")

sp2 <- sp2 + facet_grid(Date ~ yard)

sp2 <- sp2 + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5), axis.title=element_text(size=16), axis.text = element_text(size = 10), legend.text = element_text(size = 14), legend.title = element_text(size = 16), plot.caption = element_text(size=10))

sp2 <- sp2 + xlab(X_LAB1b) + ylab(Y_LAB1b)

sp2 <- sp2 + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"))

#########################
legend <- get_legend(sp1 + theme(legend.position="right"))

sp <- plot_grid(sp1 + theme(legend.position="none"),sp2 + theme(legend.position="none"), align = "v", axis = "bt", ncol = 2, labels = "AUTO", rel_heights = c(1, 1), rel_widths = c(1, 1))
sp <- plot_grid(sp, legend, align = "v", axis = "bt", ncol = 2, rel_widths = c(1, 0.25))

sp <- sp + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0.5)) + labs(caption = CAPTION1)

ggsave(sp, file=("KLee_ZBQ_full_rtpcr1_sp_vcopy_time_condition.pdf"), height=H_FIG1,width=W_FIG1,units="mm",dpi=300)



###########################################
###### make the barplot

###### change/make N and Npos values to suit groupings of interest
d$tmp_grouping <- paste(d$yard) ## exclude colour/shape variable

d$N <- ave(d$n, d$tmp_grouping, FUN = sum, na.rm = T) ## sum total present in group
d$N <- d$N/4 ## due to data reformat (duplciating data entries for each of 4 dates, we must divide totals by 4)
d1 <- aggregate(npos ~ yard + queen_lineage_og + N + tmp_grouping, data = d, FUN = sum, na.rm = T)
d1$npos <- d1$npos/4 ## due to data reformat (duplciating data entries for each of 4 dates, we must divide totals by 4)

d1$percentpos <- (d1$npos/d1$N)*100 ## calculate percent positive
d1$percentpos <- round(d1$percentpos, 0) ## round percentage to no decimals

d1$tmp_grouping <- paste(d1$yard,"(",d1$N,")") ## regroup data to be x axis (must include colour/shape variable and exclude facet variable)

bc <- ggplot(d1, aes(tmp_grouping, percentpos, fill=queen_lineage_og)) + geom_bar(stat="identity") + scale_fill_viridis_d(option="turbo") + ylim(0, 100)

bc <- bc + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5), axis.title=element_text(size=14), axis.text = element_text(size = 10), legend.text = element_text(size = 12), legend.title = element_text(size = 14), plot.caption = element_text(size=10))

bc <- bc + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0.5)) + labs(fill = COLOURLEGENDTITLE2, caption = CAPTION2) + xlab(X_LAB2) + ylab(Y_LAB2)

bc <- bc + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"))

ggsave(bc, file=("KLee_ZBQ_full_rtpcr1_bc_percent_condition.pdf"), height=H_FIG2,width=W_FIG2,units="mm",dpi=300)

####################################################################
########### violin plot with mean and median #################
d$vcopyperbee <- gsub("^0$", NA, d$vcopyperbee)
d$vcopyperbee <- as.numeric(d$vcopyperbee)

vp <- ggplot(d, aes(yard, vcopyperbee)) + geom_violin() +
  theme(axis.title.x=element_text(size=12), axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5)) +
  stat_summary(fun=mean, colour="black", geom="point", shape=3, size=2) +
  stat_summary(fun=median, colour="black", geom="point", shape=4, size=2) +
  geom_sina(aes(colour=queen_lineage_og, shape=qstatus), alpha = 0.75, position = "dodge", maxwidth=0.75, size = 3) + scale_color_viridis_d(option="turbo")

vp <- vp + theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5), axis.title=element_text(size=14), axis.text = element_text(size = 10), legend.text = element_text(size = 12), legend.title = element_text(size = 14), plot.caption = element_text(size=10))

vp <- vp + theme(plot.caption.position = "plot", plot.caption = element_text(hjust = 0.5)) + labs(fill = COLOURLEGENDTITLE3, caption = CAPTION3) + xlab(X_LAB3) + ylab(Y_LAB3)

vp <- vp + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent"), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), panel.grid.major = element_line(size = 0.1, linetype = 'solid', colour = "grey"))

ggsave(vp, file=("KLee_ZBQ_full_rtpcr1_vp_vcopy_condition.pdf"), height=H_FIG3, width=W_FIG3,units="mm",dpi=300)
