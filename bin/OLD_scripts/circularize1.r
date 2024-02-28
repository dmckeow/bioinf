#!/usr/bin/env Rscript

# Libraries
library(tidyverse) #
library(viridis) #
library(patchwork) #
library(hrbrthemes) #
library(circlize) #
library(chorddiag)  #devtools::install_github("mattflor/chorddiag") #
library(bio3d) #
library(optparse) #
library(ggplot2) #

# Define the command line options and arguments
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Path or name of input multiple alignment fasta file"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output figure file name")
  )

# Parse the command line options and arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

########################################################
########################################################

### LOAD the alignment to calculate percent ids
aln <- read.fasta(opt$input)

data <- seqidentity(aln, normalize=TRUE, similarity=FALSE)

# Set duplicates in the upper triangle to 0
data[lower.tri(data)] <- 0
## Set the diagonal (self alignments) to 0
diag(data) <- 0

data <- as.data.frame(data)

# shorten names
###colnames(data2) <- c("Africa", "East Asia", "Europe", "Latin Ame.",   "North Ame.",   "Oceania", "South Asia", "South East Asia", "Soviet Union", "West.Asia")
###rownames(data2) <- colnames(data2)

# I need a long format
data_long <- data %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname)

# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

# color palette NOTE that here is given the number of categories that must match the number of variables in the data
nvariables <- nrow(data)
mycolor <- viridis(nvariables, alpha = 1, begin = 0, end = 1, option = "D")
mycolor <- mycolor[sample(1:nvariables)]

# Base plot
plot <- chordDiagram(
  x = data_long, 
  grid.col = mycolor,
  transparency = 0.25,
  #directional = 1,
  #direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  #link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = FALSE)

# Add text and axis
plot <- plot + circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    plot <- plot + circos.text(
      x = mean(xlim), 
      y = 3.2, 
      labels = sector.index, 
      facing = "bending", 
      cex = 0.8
      )

    # Add graduation on axis
    plot <- plot + circos.axis(
      h = "top", 
      major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>10, yes = 2, no = 1)), 
      minor.ticks = 1, 
      major.tick.length = 0.5,
      labels.niceFacing = FALSE)
  }
)

ggsave(plot, paste0(opt$output, ".circosplot", format(Sys.time(), ".%Y-%m-%d"), ".png"), dpi=300)
