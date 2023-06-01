#!/usr/bin/env Rscript

##### set working directory - CHANGE TOPATH TO YOUR finalresult/IVEX002
setwd("/shared/projects/phaeoexplorer_virus/phaeoex_screen/finalresult/IVEX002_SEP22")

##### load ggplot
library(ggplot2)
library(grid)
library(tidyverse)
library(qpdf)
library(forcats)
library(dplyr)
library(cowplot)

####################################################
##### scatterplot, all genomes together
data <- read.delim("IVEX002_final_004_rbitscore_scatterplot", sep = "\t", header = T)

data <- arrange(data, sort) ## sort data by sort column with numbers setting layer of ncvog category, with 1 being the bottom layer

sga <- ggplot(data, aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.2, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3))) + theme(legend.position="bottom",legend.spacing.y=unit(-4.0, "cm")) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

ggsave(sga, filename="IVEX002_final_004_rbitscore_scatterplotall.pdf",height=450,width=600,units="mm",dpi=300)

################################################
##### scatterplot of rbitscores viral vs cellular, one graph per genome
##### note that this will require editing to match number of genomes being included - e.g. in this section I plotted the genomes in groups of max 16 per page (e.g. sg1 to sg16)
## this created 4 plots (for 60 genome, the last plot having 12 genomes instead of 16), with each section below generating one
inputs <- list.files(path = ".", pattern = "^IVEX002_final_004_rbitscore_scatterplot_.*$")

######## multiple plots per page, for main figure, manually selecting which genomes to include
data_list <- lapply(inputs, read.table, sep = '\t', header = T)

#### note that only first plot must include legend
sg1 <- ggplot(data_list[[1]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3))) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18), legend.text = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore") + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg2 <- ggplot(data_list[[2]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3))) + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg3 <- ggplot(data_list[[3]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg4 <- ggplot(data_list[[4]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg5 <- ggplot(data_list[[5]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg6 <- ggplot(data_list[[6]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg7 <- ggplot(data_list[[7]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg8 <- ggplot(data_list[[8]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg9 <- ggplot(data_list[[9]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg10 <- ggplot(data_list[[10]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg11 <- ggplot(data_list[[11]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg12 <- ggplot(data_list[[12]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg13 <- ggplot(data_list[[13]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg14 <- ggplot(data_list[[14]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg15 <- ggplot(data_list[[15]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg16 <- ggplot(data_list[[16]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

legend <- get_legend(sg1 + theme(legend.position="bottom"))
msg1 <- plot_grid(sg1 + theme(legend.position="none"),sg2,sg3,sg4,sg5,sg6,sg7,sg8,sg9,sg10,sg11,sg12,sg13,sg14,sg15,sg16, align = "h", axis = "bt", ncol = 4, labels = "AUTO")
msg1 <- plot_grid(msg1 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), legend, ncol = 1, rel_heights = c(1, .25))
ggsave(msg1, filename="IVEX002_final_004_rbitscore_scatterplotgroup001.pdf",height=650,width=600,units="mm",dpi=300)

#######################

sg17 <- ggplot(data_list[[17]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg18 <- ggplot(data_list[[18]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg19 <- ggplot(data_list[[19]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg20 <- ggplot(data_list[[20]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg21 <- ggplot(data_list[[21]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg22 <- ggplot(data_list[[22]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg23 <- ggplot(data_list[[23]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg24 <- ggplot(data_list[[24]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg25 <- ggplot(data_list[[25]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg26 <- ggplot(data_list[[26]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg27 <- ggplot(data_list[[27]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg28 <- ggplot(data_list[[28]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg29 <- ggplot(data_list[[29]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg30 <- ggplot(data_list[[30]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg31 <- ggplot(data_list[[31]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg32 <- ggplot(data_list[[32]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

msg2 <- plot_grid(sg17,sg18,sg19,sg20,sg21,sg22,sg23,sg24,sg25,sg26,sg27,sg28,sg29,sg30,sg31,sg32, align = "h", axis = "bt", ncol = 4, labels = "AUTO")
msg2 <- plot_grid(msg2 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), legend, ncol = 1, rel_heights = c(1, .25))
ggsave(msg2, filename="IVEX002_final_004_rbitscore_scatterplotgroup002.pdf",height=650,width=600,units="mm",dpi=300)

#######################

sg33 <- ggplot(data_list[[33]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg34 <- ggplot(data_list[[34]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg35 <- ggplot(data_list[[35]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg36 <- ggplot(data_list[[36]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg37 <- ggplot(data_list[[37]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg38 <- ggplot(data_list[[38]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg39 <- ggplot(data_list[[39]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg40 <- ggplot(data_list[[40]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg41 <- ggplot(data_list[[41]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg42 <- ggplot(data_list[[42]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg43 <- ggplot(data_list[[43]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg44 <- ggplot(data_list[[44]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg45 <- ggplot(data_list[[45]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg46 <- ggplot(data_list[[46]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg47 <- ggplot(data_list[[47]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg48 <- ggplot(data_list[[48]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

msg3 <- plot_grid(sg33,sg34,sg35,sg36,sg37,sg38,sg39,sg40,sg41,sg42,sg43,sg44,sg45,sg46,sg47,sg48, align = "h", axis = "bt", ncol = 4, labels = "AUTO")
msg3 <- plot_grid(msg3 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), legend, ncol = 1, rel_heights = c(1, .25))
ggsave(msg3, filename="IVEX002_final_004_rbitscore_scatterplotgroup003.pdf",height=650,width=600,units="mm",dpi=300)

#######################

sg49 <- ggplot(data_list[[49]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg50 <- ggplot(data_list[[50]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg51 <- ggplot(data_list[[51]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg52 <- ggplot(data_list[[52]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg53 <- ggplot(data_list[[53]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg54 <- ggplot(data_list[[54]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg55 <- ggplot(data_list[[55]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg56 <- ggplot(data_list[[56]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg57 <- ggplot(data_list[[57]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg58 <- ggplot(data_list[[58]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg59 <- ggplot(data_list[[59]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg60 <- ggplot(data_list[[60]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg61 <- ggplot(data_list[[60]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg62 <- ggplot(data_list[[60]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg63 <- ggplot(data_list[[60]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg64 <- ggplot(data_list[[60]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")


msg4 <- plot_grid(sg49,sg50,sg51,sg52,sg53,sg54,sg55,sg56,sg57,sg58,sg59,sg60,sg61,sg62,sg63,sg64, align = "h", axis = "bt", ncol = 4, labels = "AUTO")
msg4 <- plot_grid(msg4 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), legend, ncol = 1, rel_heights = c(1, .25))
ggsave(msg4, filename="IVEX002_final_004_rbitscore_scatterplotgroup004.pdf",height=650,width=600,units="mm",dpi=300)


#######################

sg65 <- ggplot(data_list[[49]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg66 <- ggplot(data_list[[50]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg67 <- ggplot(data_list[[51]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg68 <- ggplot(data_list[[52]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg69 <- ggplot(data_list[[53]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg70 <- ggplot(data_list[[54]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg71 <- ggplot(data_list[[55]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

sg72 <- ggplot(data_list[[56]] %>% arrange(sort), aes(vrbitscore, rbitscore)) + geom_jitter(size = 0.7, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3)))  + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")

msg5 <- plot_grid(sg65,sg66,sg67,sg68,sg69,sg70,sg71,sg72, align = "h", axis = "bt", ncol = 4, labels = "AUTO")
msg5 <- plot_grid(msg5 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), legend, ncol = 1, rel_heights = c(1, .25))
ggsave(msg5, filename="IVEX002_final_004_rbitscore_scatterplotgroup005.pdf",height=650,width=600,units="mm",dpi=300)

#######################


#######################
## plot a specific set of genomes for a particular paper figure (note that the genome data does not need to loaded again - just use the existing sg dataframes)
#msgport <- plot_grid(sg42,sg43,sg44,sg46, align = "h", axis = "bt", ncol = 2, labels = "AUTO")
#msgport <- plot_grid(msgport + theme(legend.position="none", plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(size = 16), axis.text = element_text(size = 14), strip.text.x = element_text(size = 18)), legend, ncol = 1, rel_heights = c(1, .25))
#ggsave(msgport, filename="Figure_S1_Pfluv_paper.pdf",height=500,width=450,units="mm",dpi=300)

#####################################################
######## one figure per genome, 1 per page in pdf, for supplementary

CAPTION = str_wrap("Figure S1. Scatterplot of cellular and viral hits against algal proteins. Points represent the relative BLASTP scores of algal proteins aligned against cellular (y axis, NR database with virus proteins excluded) and viral (x axis, RVDB) databases. Datapoint colours indicate the categories of virus hits (NCV = any Nucleocytoviricota gene, NCVOG_groups_1-3 = conserved NCV core genes of groups 1-3), non-NCV = any gene of any non-Nucleocytoviricota virus. For reference, public genomes with known viral content are included, as indicated by PUBLIC_.", 100)

for (f in inputs) {
  data <- read.delim(f, sep = "\t", header = T)
  data <- arrange(data, sort)
  sg <- ggplot(data, aes(vrbitscore, rbitscore)) + geom_jitter(size = 1, aes(colour = ncvog), position = position_jitter(width = 0.02, height = 0.02)) + facet_wrap(~ genome) + theme_classic() + scale_colour_viridis_d(name="category") + guides(colour = guide_legend(override.aes = list(size=3))) + labs(caption=CAPTION) + theme(plot.caption = element_text(hjust = 0.5)) + theme(legend.position="bottom",legend.spacing.y=unit(-4.0, "cm")) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x ="viral rbitscore", y = "cellular rbitscore")
  ggsave(sg, file=paste0("single_",f,".pdf"), height=225, width=300, units="mm", dpi=300)
}

#### concatenate PDF files
qpdf::pdf_combine(input = list.files(full.names=TRUE,pattern="^single_IVEX002_final_004_rbitscore_scatterplot_.*.pdf$"), output = "IVEX002_final_004_rbitscore_scatterplotallseries.pdf")

##### get rid of unncesscary pdf copies
junk <- list.files(path = ".", pattern = "^single_IVEX002_final_004_rbitscore_scatterplot_.*$")
file.remove(junk)

###################################################
############### heatmaps
##### all genomes, core genes
data <- read.delim("IVEX002_final_002_ncvog_count_genomes_heatmap", sep = "\t", header = T)

### stop R reordering rows
data$genome <- as.character(data$genome)
data$genome <- factor(data$genome, levels=unique(data$genome))

### categorise based on core gene count
data$countgroup <- cut(data$count, breaks = c(0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 19, Inf), labels = c("1","2","3","4","5","6","7","8","9","10","11-15","16-19","20+"), include.lowest=TRUE, right=TRUE)
data$legendncvog <- paste(data$legend,"|",data$ncvog)

hm <- ggplot(data, aes(legendncvog, genome, fill=countgroup)) + geom_tile(colour="black", size=0.2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.title = element_text(size=10, hjust = 0.5)) + scale_y_discrete(expand=c(0,0)) + labs(x="NCLDV core gene",y="genome", fill = "core gene count") + facet_grid(. ~ group, space="free_x", scales="free_x", switch="y") + ggtitle("core gene group") + scale_fill_viridis_d(option="magma") + theme(axis.text=element_text(size = 7)) + geom_text(aes(label = count), color="grey", size = 2) + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), axis.ticks = element_blank())

ggsave(hm, filename="IVEX002_final_002_ncvog_count_genomes_heatmap.pdf",height=300,width=300,units="mm",dpi=300)

##### contigs, per genome, core genes
inputs <- list.files(path = ".", pattern = "^IVEX002_final_001_.*_ncvog_count_contigs_heatmap$")

for (f in inputs) {
  data <- read.delim(f, sep = "\t", header = T)
  #data <- arrange(data, contig)
  TITLE <- unique(data$genome)

  ### stop R reordering rows
  data$contig <- as.character(data$contig)
  data$contig <- factor(data$contig, levels=unique(data$contig))

  ### categorise based on core gene count
  data$countgroup <- cut(data$count, breaks = c(0,1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 19, Inf), labels = c("1","2","3","4","5","6","7","8","9","10","11-15","16-19","20+"), include.lowest=TRUE, right=TRUE)
  data$legendncvog <- paste(data$legend,"|",data$ncvog)

  hm <- ggplot(data, aes(legendncvog, contig, fill=countgroup)) + geom_tile(colour="black", size=0.2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.title = element_text(size=10, hjust = 0.5)) + scale_y_discrete(expand=c(0,0)) + labs(x="NCLDV core gene",y="contig") + facet_grid(. ~ group, space="free_x", scales="free_x", switch="y") + ggtitle("core gene group") + scale_fill_viridis_d(option="magma") + theme(axis.text=element_text(size = 7)) + ggtitle(TITLE) + geom_text(aes(label = count), color="grey", size = 2) + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), axis.ticks = element_blank())

  ggsave(hm, file=paste0("single_",f,".pdf"), height=600,width=300,units="mm",dpi=300)
}
qpdf::pdf_combine(input = list.files(full.names=TRUE,pattern="^single_IVEX002_final_001_.*.pdf$"), output = "IVEX002_final_001_all_ncvog_count_contigs_heatmapseries.pdf")

  junk <- list.files(path = ".", pattern = "^single_IVEX002_final_001_.*$")
  file.remove(junk)

##################################################
##### all genomes, virus taxa
data <- read.delim("IVEX002_final_003_virus_taxa_percent_count_genomes_heatmap", sep = "\t", header = T)

hm <- ggplot(data, aes(genome,taxa3,fill=percent)) + geom_tile(colour="black", size=0.2) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + scale_fill_viridis_c(option="magma") + theme(axis.text=element_text(size = 7)) + labs(x="genome",y="virus taxa") + facet_grid(taxa1 + taxa2 ~., space="free", scales="free", switch="x") + theme(strip.text.y = element_text(angle = 0))

ggsave(hm, filename="IVEX002_final_003_virus_taxa_percent_count_genomes_heatmap.pdf",height=600,width=400,units="mm",dpi=300)

#####

