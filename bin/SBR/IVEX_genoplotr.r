#!/usr/bin/env Rscript

##### load packages
library(ade4)
library(grid)
library(genoPlotR)

##### required input file per molecule being aligned (df1, etc):
### name start end strand col
### feat1 2 600 -1 blue

###### load data and prepare as dna segments
### duplicate following step for each molecule, increasing output names by 1 each time
df1 <- read.delim("df1", sep = "\t", header = T)
dna_seg1 <- dna_seg(df1)

df2 <- read.delim("df2", sep = "\t", header = T)
dna_seg2 <- dna_seg(df2)

df3 <- read.delim("df3", sep = "\t", header = T)
dna_seg3 <- dna_seg(df3)

df4 <- read.delim("df4", sep = "\t", header = T)
dna_seg4 <- dna_seg(df4)

df5 <- read.delim("df5", sep = "\t", header = T)
dna_seg5 <- dna_seg(df5)

df6 <- read.delim("df6", sep = "\t", header = T)
dna_seg6 <- dna_seg(df6)

df7 <- read.delim("df7", sep = "\t", header = T)
dna_seg7 <- dna_seg(df7)


### change to include all dna_segs
dna_segs <- list(dna_seg1, dna_seg2, dna_seg3, dna_seg4, dna_seg5, dna_seg6, dna_seg7)

### set names for each molecule
names <- c("FsV-158", "EVE7", "EsV-1", "EVE12", "FsV-158", "EVE22", "EsV-1")

names(dna_segs) <- names

###### prepare comparison objects
###### comparison file input (from promer, blast, etc):
### start1 end1 start2 end2 col
### must de duplicated for each comparison, i.e. pair of molecules to be compared (adjacent; e.g. genome 1 and 2, 2 and 3, etc)

comparison1 <- read_comparison_from_tab("comparison1.tab")
comparison2 <- read_comparison_from_tab("comparison2.tab")
comparison3 <- read_comparison_from_tab("comparison3.tab")
comparison4 <- read_comparison_from_tab("comparison4.tab")
comparison5 <- read_comparison_from_tab("comparison5.tab")
comparison6 <- read_comparison_from_tab("comparison6.tab")

### change to include all comparisons
comparisons <- list(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6)

### parameters to show all feature labels
annots <- lapply(dna_segs, function(x){
mid <- middle(x)
annot <- annotation(x1=mid, text=x$name, rot=90)
})

### set color intensity scale on comparisons so higher % identity darker color (per comparison)
comparisons[[1]]$col <- apply_color_scheme(comparisons[[1]]$values, color_scheme="grey", transparency = 0.7)
comparisons[[2]]$col <- apply_color_scheme(comparisons[[2]]$values, color_scheme="grey", transparency = 0.7)
comparisons[[3]]$col <- apply_color_scheme(comparisons[[3]]$values, color_scheme="grey", transparency = 0.7)
comparisons[[4]]$col <- apply_color_scheme(comparisons[[4]]$values, color_scheme="grey", transparency = 0.7)
comparisons[[5]]$col <- apply_color_scheme(comparisons[[5]]$values, color_scheme="grey", transparency = 0.7)
comparisons[[6]]$col <- apply_color_scheme(comparisons[[6]]$values, color_scheme="grey", transparency = 0.7)

### apply reverse orientation if needed
## these values are the specific lengths of each molecule - they can be found in the last coordinate in the relevant dataframe
xlims <- list(c(1, 152847), c(314398, 1), c(1, 333659), c(320239, 1), c(1, 152847), c(227719, 1), c(1, 333659))

##### generate the plot - it will have default name - change it after
plot_gene_map(dna_segs = dna_segs, comparisons = comparisons, annotation_cex = 0.3, dna_seg_label_cex = 1, scale_cex = 0.5, annotations = annots, gene_type ="side_blocks", dna_seg_scale = TRUE, scale = FALSE, offsets = NULL, xlims = xlims)
