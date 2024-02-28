#!/usr/bin/env Rscript

##### set working directory - CHANGE TOPATH TO YOUR finalresult/IVEX002
setwd("/shared/projects/phaeoexplorer_virus/phaeoex_screen/finalresult/IVEX003")

##### load ggplot
library(ggplot2)
library(grid)
library(qpdf)
library(forcats)
library(cowplot)
library(treeio)
library(ggtree)
library(ape)
library(phangorn)
library(tidyverse); library(dplyr)


####################################################


raxml1 <- "/shared/projects/phaeoexplorer_virus/phaeoex_screen/finalresult/IVEX003/RAxML_bipartitionsBranchLabels.DNApolB_MAR22_60genomes.phex.public.FINAL.faa_200bs.raxml"

raxml2 <- read.raxml(raxml1)

raxml3 <- as.phylo(raxml2)

raxml4 <- phangorn::midpoint(raxml3, node.labels = "support")
raxml5 <- reorder(raxml4)

raxml6 <- as.treedata(raxml5)

traxml3 <- as_tibble(raxml3)
traxml5 <- as_tibble(raxml5)

traxml5$bootstrap <- traxml3$bootstrap
raxml7 <- as.treedata(traxml5)


tr1 <- ggtree(raxml7, layout="circular") + geom_text(aes(label=bootstrap, x=branch, vjust=0.5, size=2)) + geom_tiplab(size=3) + theme_tree2(legend.position='right')

#scale_fill_continuous(low='blue', high='red')
ggsave(tr1, filename="tr1.pdf",height=200,width=200,units="mm",dpi=300)
