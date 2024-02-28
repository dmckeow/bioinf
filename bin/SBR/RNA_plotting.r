#!/usr/bin/env Rscript
##### load ggplot
library(ggplot2)
library(grid)
library(tidyverse)
library(qpdf)
library(forcats)
library(dplyr)
library(cowplot)

#### set working directory
setwd("/shared/projects/phaeoexplorer_virus/phaeoex_screen/finalresult/IVEX004")

##### HEATMAP - plot TPM data against contig/scaffold loci

data1 <- read.delim("Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig7_heatmap", sep="\t", header=T) ### contig file
data2 <- read.delim("Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig12_heatmap", sep="\t", header=T) ### contig file
data3 <- read.delim("Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig21_heatmap", sep="\t", header=T) ### contig file
data4 <- read.delim("Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig22_heatmap", sep="\t", header=T) ### contig file
data5 <- read.delim("Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig40_heatmap", sep="\t", header=T) ### contig file
data6 <- read.delim("Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig94_heatmap", sep="\t", header=T) ### contig file


  ## stop ggplot reordering rows

  data1$START <- as.character(data1$START)
  data1$START <- factor(data1$START, levels=unique(data1$START))

  data2$START <- as.character(data2$START)
  data2$START <- factor(data2$START, levels=unique(data2$START))

    data3$START <- as.character(data3$START)
    data3$START <- factor(data3$START, levels=unique(data3$START))

      data4$START <- as.character(data4$START)
      data4$START <- factor(data4$START, levels=unique(data4$START))

        data5$START <- as.character(data5$START)
        data5$START <- factor(data5$START, levels=unique(data5$START))

          data6$START <- as.character(data6$START)
          data6$START <- factor(data6$START, levels=unique(data6$START))


  ## MANUAL - change scale fill limits to those in whole genome see *_RNAseq_TPM_LOG2p1_min_ave_max file

  hm1 <- ggplot(data1, aes(CONDITION, START)) + geom_raster(aes(fill=TPM), hjust=0.5, vjust=0.5) + theme(axis.text.x=element_blank(), axis.text.y=element_text(size=0.5),plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(y="GENES", x="TPM", fill = "TPM per gene, Log2(+1)") + scale_fill_viridis_c(option="magma", breaks=c(0,1,14.5), labels=c("0","1","14.5"), limits=c(0,14.6)) + facet_grid(. ~ CONDITION, space="free_x", scales="free_x", switch="x")

  bc1 <- ggplot(data1, aes(NUMREADS, START, fill=CONDITION)) + geom_bar(stat="identity") + scale_fill_viridis_d() + theme(axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(fill = "no. of reads per gene")

  rc1 <- ggplot(data1, aes(CONTIG, START)) + geom_raster(aes(fill=REPEATS), hjust=0.5, vjust=0.5) + scale_fill_viridis_c(option="magma") + theme(axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(x="REPEATS", fill = "no. of repeats in genome")

##### store legends to add back later in different position
  l1_1 <- get_legend(hm1 + theme(legend.position="right"))
  l1_2 <- get_legend(bc1 + theme(legend.position="right"))
  l1_3 <- get_legend(rc1 + theme(legend.position="right"))
  l1 <- plot_grid(l1_1,l1_2,l1_3, ncol = 1)

###########################
hm2 <- ggplot(data2, aes(CONDITION, START)) + geom_raster(aes(fill=TPM), hjust=0.5, vjust=0.5) + theme(axis.text.x=element_blank(), axis.text.y=element_text(size=0.5),plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(y="GENES", x="TPM", fill = "TPM per gene, Log2(+1)") + scale_fill_viridis_c(option="magma", breaks=c(0,1,14.5), labels=c("0","1","14.5"), limits=c(0,14.6)) + facet_grid(. ~ CONDITION, space="free_x", scales="free_x", switch="x")

bc2 <- ggplot(data2, aes(NUMREADS, START, fill=CONDITION)) + geom_bar(stat="identity") + scale_fill_viridis_d() + theme(axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(fill = "no. of reads per gene")

rc2 <- ggplot(data2, aes(CONTIG, START)) + geom_raster(aes(fill=REPEATS), hjust=0.5, vjust=0.5) + scale_fill_viridis_c(option="magma") + theme(axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(x="REPEATS", fill = "no. of repeats in genome")

##### store legends to add back later in different position
l2_1 <- get_legend(hm2 + theme(legend.position="right"))
l2_2 <- get_legend(bc2 + theme(legend.position="right"))
l2_3 <- get_legend(rc2 + theme(legend.position="right"))
l2 <- plot_grid(l2_1,l2_2,l2_3, ncol = 1)

###########################
hm3 <- ggplot(data3, aes(CONDITION, START)) + geom_raster(aes(fill=TPM), hjust=0.5, vjust=0.5) + theme(axis.text.x=element_blank(), axis.text.y=element_text(size=0.5),plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(y="GENES", x="TPM", fill = "TPM per gene, Log2(+1)") + scale_fill_viridis_c(option="magma", breaks=c(0,1,14.5), labels=c("0","1","14.5"), limits=c(0,14.6)) + facet_grid(. ~ CONDITION, space="free_x", scales="free_x", switch="x")

bc3 <- ggplot(data3, aes(NUMREADS, START, fill=CONDITION)) + geom_bar(stat="identity") + scale_fill_viridis_d() + theme(axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(fill = "no. of reads per gene")

rc3 <- ggplot(data3, aes(CONTIG, START)) + geom_raster(aes(fill=REPEATS), hjust=0.5, vjust=0.5) + scale_fill_viridis_c(option="magma") + theme(axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(x="REPEATS", fill = "no. of repeats in genome")

##### store legends to add back later in different position
l3_1 <- get_legend(hm3 + theme(legend.position="right"))
l3_2 <- get_legend(bc3 + theme(legend.position="right"))
l3_3 <- get_legend(rc3 + theme(legend.position="right"))
l3 <- plot_grid(l3_1,l3_2,l3_3, ncol = 1)
###########################
hm4 <- ggplot(data4, aes(CONDITION, START)) + geom_raster(aes(fill=TPM), hjust=0.5, vjust=0.5) + theme(axis.text.x=element_blank(), axis.text.y=element_text(size=0.5),plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(y="GENES", x="TPM", fill = "TPM per gene, Log2(+1)") + scale_fill_viridis_c(option="magma", breaks=c(0,1,14.5), labels=c("0","1","14.5"), limits=c(0,14.6)) + facet_grid(. ~ CONDITION, space="free_x", scales="free_x", switch="x")

bc4 <- ggplot(data4, aes(NUMREADS, START, fill=CONDITION)) + geom_bar(stat="identity") + scale_fill_viridis_d() + theme(axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(fill = "no. of reads per gene")

rc4 <- ggplot(data4, aes(CONTIG, START)) + geom_raster(aes(fill=REPEATS), hjust=0.5, vjust=0.5) + scale_fill_viridis_c(option="magma") + theme(axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(x="REPEATS", fill = "no. of repeats in genome")

##### store legends to add back later in different position
l4_1 <- get_legend(hm4 + theme(legend.position="right"))
l4_2 <- get_legend(bc4 + theme(legend.position="right"))
l4_3 <- get_legend(rc4 + theme(legend.position="right"))
l4 <- plot_grid(l4_1,l4_2,l4_3, ncol = 1)

###########################
hm5 <- ggplot(data5, aes(CONDITION, START)) + geom_raster(aes(fill=TPM), hjust=0.5, vjust=0.5) + theme(axis.text.x=element_blank(), axis.text.y=element_text(size=0.5),plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(y="GENES", x="TPM", fill = "TPM per gene, Log2(+1)") + scale_fill_viridis_c(option="magma", breaks=c(0,1,14.5), labels=c("0","1","14.5"), limits=c(0,14.6)) + facet_grid(. ~ CONDITION, space="free_x", scales="free_x", switch="x")

bc5 <- ggplot(data5, aes(NUMREADS, START, fill=CONDITION)) + geom_bar(stat="identity") + scale_fill_viridis_d() + theme(axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(fill = "no. of reads per gene")

rc5 <- ggplot(data5, aes(CONTIG, START)) + geom_raster(aes(fill=REPEATS), hjust=0.5, vjust=0.5) + scale_fill_viridis_c(option="magma") + theme(axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(x="REPEATS", fill = "no. of repeats in genome")

##### store legends to add back later in different position
l5_1 <- get_legend(hm5 + theme(legend.position="right"))
l5_2 <- get_legend(bc5 + theme(legend.position="right"))
l5_3 <- get_legend(rc5 + theme(legend.position="right"))
l5 <- plot_grid(l5_1,l5_2,l5_3, ncol = 1)

###########################

hm6 <- ggplot(data6, aes(CONDITION, START)) + geom_raster(aes(fill=TPM), hjust=0.5, vjust=0.5) + theme(axis.text.x=element_blank(), axis.text.y=element_text(size=0.5),plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(y="GENES", x="TPM", fill = "TPM per gene, Log2(+1)") + scale_fill_viridis_c(option="magma", breaks=c(0,1,14.5), labels=c("0","1","14.5"), limits=c(0,14.6)) + facet_grid(. ~ CONDITION, space="free_x", scales="free_x", switch="x")

bc6 <- ggplot(data6, aes(NUMREADS, START, fill=CONDITION)) + geom_bar(stat="identity") + scale_fill_viridis_d() + theme(axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(fill = "no. of reads per gene")

rc6 <- ggplot(data6, aes(CONTIG, START)) + geom_raster(aes(fill=REPEATS), hjust=0.5, vjust=0.5) + scale_fill_viridis_c(option="magma") + theme(axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_text(size=8), plot.margin = unit(c(0, 0, 0, 0), "cm")) + labs(x="REPEATS", fill = "no. of repeats in genome")

##### store legends to add back later in different position
l6_1 <- get_legend(hm6 + theme(legend.position="right"))
l6_2 <- get_legend(bc6 + theme(legend.position="right"))
l6_3 <- get_legend(rc6 + theme(legend.position="right"))
l6 <- plot_grid(l6_1,l6_2,l6_3, ncol = 1)

###########################
  ##### MANUAL - change position see __RNAseq_TPM_LOG2p1_EVE7_coords for line positions (you cant use the loci because discrete variable)

  hm1 <- hm1 + annotate("segment", x=1, xend=1, y=2040, yend=2356, color="white", size=2)
  hm2 <- hm2 + annotate("segment", x=1, xend=1, y=15, yend=316, color="white", size=2)
  hm3 <- hm3 + annotate("segment", x=1, xend=1, y=197, yend=322, color="white", size=2)
  hm4 <- hm4 + annotate("segment", x=1, xend=1, y=593, yend=830, color="white", size=2)
  hm5 <- hm5 + annotate("segment", x=c(1,1), xend=c(1,1), y=c(179,916), yend=c(205,1069), color="white", size=2)
  hm6 <- hm6 + annotate("segment", x=1, xend=1, y=39, yend=87, color="white", size=2)

  #### combine all charts, aligned and without legends
  hm1 <- plot_grid(hm1 + theme(legend.position="none"),bc1 + theme(legend.position="none") ,rc1 + theme(legend.position="none"), align = "h",axis = "bt", ncol=3, rel_widths = c(1, 0.5, 0.25))
  hm2 <- plot_grid(hm2 + theme(legend.position="none"),bc2 + theme(legend.position="none") ,rc2 + theme(legend.position="none"), align = "h",axis = "bt", ncol=3, rel_widths = c(1, 0.5, 0.25))
  hm3 <- plot_grid(hm3 + theme(legend.position="none"),bc3 + theme(legend.position="none") ,rc3 + theme(legend.position="none"), align = "h",axis = "bt", ncol=3, rel_widths = c(1, 0.5, 0.25))
  hm4 <- plot_grid(hm4 + theme(legend.position="none"),bc4 + theme(legend.position="none") ,rc4 + theme(legend.position="none"), align = "h",axis = "bt", ncol=3, rel_widths = c(1, 0.5, 0.25))
  hm5 <- plot_grid(hm5 + theme(legend.position="none"),bc5 + theme(legend.position="none") ,rc5 + theme(legend.position="none"), align = "h",axis = "bt", ncol=3, rel_widths = c(1, 0.5, 0.25))
  hm6 <- plot_grid(hm6 + theme(legend.position="none"),bc6 + theme(legend.position="none") ,rc6 + theme(legend.position="none"), align = "h",axis = "bt", ncol=3, rel_widths = c(1, 0.5, 0.25))

#### add back in legends
hm1 <- plot_grid(hm1 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), l1 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), ncol = 2)
hm2 <- plot_grid(hm2 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), l1 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), ncol = 2)
hm3 <- plot_grid(hm3 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), l1 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), ncol = 2)
hm4 <- plot_grid(hm4 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), l1 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), ncol = 2)
hm5 <- plot_grid(hm5 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), l1 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), ncol = 2)
hm6 <- plot_grid(hm6 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), l1 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), ncol = 2)

  ##### MANUAL - add title and caption

  CAPTION = str_wrap("Figure S9. RNA expression of contig. Transcripts per million (TPM) is plotted against genes, followed left to right by the number of reads and the number of repeats in the genome. Repeats were counted by blastn of genes against all genes in genome, counting any self match below evalue 1e-40 as a repeat, up to a max number of repeats of 500. Vertical white lines mark location of virus region. Condition abbreviations: FW=freshwater, SW=seawater", 60)

  hm1 <- hm1 + ggtitle("P-fluviatile_contig7") + labs(caption = CAPTION)
  hm2 <- hm2 + ggtitle("P-fluviatile_contig12") + labs(caption = CAPTION)
  hm3 <- hm3 + ggtitle("P-fluviatile_contig21") + labs(caption = CAPTION)
  hm4 <- hm4 + ggtitle("P-fluviatile_contig22") + labs(caption = CAPTION)
  hm5 <- hm5 + ggtitle("P-fluviatile_contig40") + labs(caption = CAPTION)
  hm6 <- hm6 + ggtitle("P-fluviatile_contig94") + labs(caption = CAPTION)


  #### modify appearance properties

  hm1 <- hm1 + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), axis.ticks = element_blank(), plot.title = element_text(color="black", size=16, hjust=0.5, face="bold"), plot.caption.position = "plot", plot.caption = element_text(hjust=0.5))

  hm2 <- hm2 + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), axis.ticks = element_blank(), plot.title = element_text(color="black", size=16, hjust=0.5, face="bold"), plot.caption.position = "plot", plot.caption = element_text(hjust=0.5))

    hm3 <- hm3 + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), axis.ticks = element_blank(), plot.title = element_text(color="black", size=16, hjust=0.5, face="bold"), plot.caption.position = "plot", plot.caption = element_text(hjust=0.5))

      hm4 <- hm4 + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), axis.ticks = element_blank(), plot.title = element_text(color="black", size=16, hjust=0.5, face="bold"), plot.caption.position = "plot", plot.caption = element_text(hjust=0.5))

        hm5 <- hm5 + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), axis.ticks = element_blank(), plot.title = element_text(color="black", size=16, hjust=0.5, face="bold"), plot.caption.position = "plot", plot.caption = element_text(hjust=0.5))

          hm6 <- hm6 + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), axis.ticks = element_blank(), plot.title = element_text(color="black", size=16, hjust=0.5, face="bold"), plot.caption.position = "plot", plot.caption = element_text(hjust=0.5))

#### save output

  ggsave(hm1, filename="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig7_heatmap.pdf", height=300, width=200, units="mm", dpi=300, bg="transparent")
  ggsave(hm2, filename="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig12_heatmap.pdf", height=300, width=200, units="mm", dpi=300, bg="transparent")
  ggsave(hm3, filename="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig21_heatmap.pdf", height=300, width=200, units="mm", dpi=300, bg="transparent")
  ggsave(hm4, filename="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig22_heatmap.pdf", height=300, width=200, units="mm", dpi=300, bg="transparent")
  ggsave(hm5, filename="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig40_heatmap.pdf", height=300, width=200, units="mm", dpi=300, bg="transparent")
  ggsave(hm6, filename="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig94_heatmap.pdf", height=300, width=200, units="mm", dpi=300, bg="transparent")

##### HEATMAP - plot TPM data against gene identity within a specific viral region (still ordered by loci) only genes expressed in at least 1 condition

#####inputs <- list.files(path = ".", pattern = "^.*EVE.*_heatmap$")
inputs <- list.files(path = ".", pattern = "^.*__RNAseq_TPM_LOG2p1_EVE.*_heatmap$")

for (f in inputs) {
  data <- read.delim(f, sep = "\t", header = T)

  TITLE <- unique(data$CONTIG)

  data$FULLNAME <- paste(data$NAME,"|",data$ABBREV,"|",data$VSALLTITLES,"|",data$SALLTITLES,"|",data$PFAM)

  CAPTION = str_wrap("Figure S10. Expressed genes from viral region of contig (all genes with non-zero TPM in any condition). Transcripts per million (TPM) is plotted against genes and TPM numerical value is shown on tiles. Colour scale shows range and average TPM from whole genome. Y axis labels show: gene name | NCVOG abbreviation | viral BLASTp salltitles | cellular BLASTp salltitles | PFAM ID and description. Condition abbreviations: FW=freshwater, SW=seawater", 100)


hm <- data %>% mutate(FULLNAME = fct_reorder(FULLNAME, TPM)) %>% ggplot(aes(CONDITION, FULLNAME)) + geom_raster(aes(fill=TPM), hjust=0.5, vjust=0.5) + theme(axis.text.x=element_blank(), axis.text.y=element_text(size=5)) + labs(y="GENES", fill = "Log2(TPM+1)") + scale_fill_viridis_c(option="magma", breaks=c(0,1,14.5), labels=c("0","1","14.5"), limits=c(0,14.6)) + facet_grid(. ~ CONDITION, space="free_x", scales="free_x", switch="x") + ggtitle(TITLE) + labs(caption = CAPTION) + geom_text(aes(label = round(TPM,2)), color="grey")

#### modify appearance properties
  hm <- hm + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), axis.ticks = element_blank(), plot.title = element_text(color="black", size=16, hjust=0.5, face="bold"), plot.caption.position = "plot", plot.caption = element_text(hjust = 0.5))

  ggsave(hm, file=paste0(f,".pdf"), height=300, width=400, units="mm", dpi=300, bg = "transparent")
}

#### concatenate PDF files

qpdf::pdf_combine(input = c("Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig7_heatmap.pdf", "Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig12_heatmap.pdf", "Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig21_heatmap.pdf", "Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig22_heatmap.pdf", "Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig40_heatmap.pdf", "Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig94_heatmap.pdf"), output = "Figure_S9_Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contigs_heatmaps.pdf")

qpdf::pdf_combine(input = c("Porterinema-fluviatile__RNAseq_TPM_LOG2p1_EVE7_heatmap.pdf", "Porterinema-fluviatile__RNAseq_TPM_LOG2p1_EVE12_heatmap.pdf", "Porterinema-fluviatile__RNAseq_TPM_LOG2p1_EVE21_heatmap.pdf", "Porterinema-fluviatile__RNAseq_TPM_LOG2p1_EVE22_heatmap.pdf", "Porterinema-fluviatile__RNAseq_TPM_LOG2p1_EVE40-1_heatmap.pdf", "Porterinema-fluviatile__RNAseq_TPM_LOG2p1_EVE40-2_heatmap.pdf", "Porterinema-fluviatile__RNAseq_TPM_LOG2p1_EVE94_heatmap.pdf"), output = "Figure_S10_Porterinema-fluviatile__RNAseq_TPM_LOG2p1_EVEs_heatmaps.pdf")
