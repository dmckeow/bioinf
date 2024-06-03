#!/usr/bin/env Rscript

# Script use: 
## Visualising the BAM files mapped to contigs using ggplot
## BAM files must be converted to a tsv numeric file using `samtools depth`
## check the wrangle datasheets header parts for infromation for what should be contained within
## each tsv file in order for this script to work.

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggalluvial)
library(RColorBrewer)
library(svglite)
library(ggrepel)
library(treemapify)

#library(reshape) #
#library(forcats) #
#library(plotrix) #
#library(tidyr) #
#library(viridis) #
#library(viridisLite) #

theme_set(theme_classic())
set.seed(1234)

message("\n\n\nTip: If you have run bioinf-read-check.py on multiple assemblies, you can combine their tabular output to use as input for this script (bioinf-read-check.r). From a directory that is a parent to all of your outputs from bioinf-read-check.py, simply do the following in your bash terminal:\n\nfor f in $(find ASSEMBLY-CANU_* -name bioinf_readcheck_allcompare.final); do cat $f ;done | sed '/.*NCBI_taxid_all,NCBI_accession_all,matching_seq_frag,NCBI_full_taxon_lineage/d' | sed -z 's/^/bctrimmedreads.fastq,correctedReads.fasta,trimmedReads.fasta,contigs,classified_unclassified,read_name,NCBI_taxid_best,kaiju_match_score,NCBI_taxid_all,NCBI_accession_all,matching_seq_frag,NCBI_full_taxon_lineage/g' > bioinf_readcheck_allcompare.final; conda activate bioinftools; bioinf-read-check.r\n\n\n")

################################################
############## IMPORT DATASHEETS ###############
################################################
df <- read.csv("bioinf_readcheck_allcompare.final", header=TRUE)

################################################
############## WRANGLE DATASHEETS ##############
################################################

#### set a function to use for rounding decimal places
specify_decimal <- function(x, k) as.numeric(trimws(format(round(x, k), nsmall=k)))


colnames(df)[1] ="orig_read_name"
colnames(df)[2] ="canu_correctedreads"
colnames(df)[3] ="canu_trimmedreads"
colnames(df)[4] ="canu_contigs"

df$raw_bctrimmed_reads <- as.numeric(sub(".*", "1", df$orig_read_name))

df$canu_correctedreads <- sub("0", "FAILED_co", df$canu_correctedreads)
df$canu_correctedreads <- sub("1", "passed_co", df$canu_correctedreads)

df$canu_trimmedreads <- sub("0", "FAILED_tr", df$canu_trimmedreads)
df$canu_trimmedreads <- sub("1", "passed_tr", df$canu_trimmedreads)

df$canu_contigs <- sub("0", "FAILED_as", df$canu_contigs)
df$canu_contigs <- sub("1", "passed_as", df$canu_contigs)


df$NCBI_full_taxon_lineage <- gsub("[^A-Za-z0-9_;.]", "_", df$NCBI_full_taxon_lineage)
df$NCBI_full_taxon_lineage <- gsub("_+", "_", df$NCBI_full_taxon_lineage)
df$NCBI_full_taxon_lineage <- sub("$", ";", df$NCBI_full_taxon_lineage)
df$NCBI_full_taxon_lineage <- gsub("NA;", ";", df$NCBI_full_taxon_lineage)
df$NCBI_full_taxon_lineage <- gsub(";+", ";", df$NCBI_full_taxon_lineage)

df$NCBI_full_taxon_lineage <- sub("Archaea;", "AR;", df$NCBI_full_taxon_lineage)
df$NCBI_full_taxon_lineage <- sub("Bacteria;", "BA.;", df$NCBI_full_taxon_lineage)
df$NCBI_full_taxon_lineage <- sub("Eukaryota;", "EU.;", df$NCBI_full_taxon_lineage)
df$NCBI_full_taxon_lineage <- sub("Viruses;", "VI.;", df$NCBI_full_taxon_lineage)

df <- within(df, NCBI_full_taxon_lineage[NCBI_full_taxon_lineage == ';'] <- 'UN.')
df <- within(df, NCBI_full_taxon_lineage[NCBI_full_taxon_lineage == ''] <- 'UN.')
df <- within(df, NCBI_full_taxon_lineage[classified_unclassified == 'U'] <- 'UN.')
df <- within(df, classified_unclassified[NCBI_full_taxon_lineage == 'UN.'] <- 'U')

## store original taxa lineage

df$orig_NCBI_full_taxon_lineage <- df$NCBI_full_taxon_lineage

#### output table for virus data
write.csv(df, "bioinf_readcheck_data_all_byreads.csv", row.names=FALSE)

##############################
## filter out only viruses

df_v <- df %>% filter(grepl('^VI.;', NCBI_full_taxon_lineage))

### reduce taxa levels (repeat to remove next lowest level)
df_v$NCBI_full_taxon_lineage <- sub(";$", "", df_v$NCBI_full_taxon_lineage)
df_v$NCBI_full_taxon_lineage <- sub("^VI.;", "", df_v$NCBI_full_taxon_lineage)
df_v$NCBI_full_taxon_lineage <- sub("^VI.$", "UN.", df_v$NCBI_full_taxon_lineage)
df_v$NCBI_full_taxon_lineage <- sub(";[^;]*;", ";", df_v$NCBI_full_taxon_lineage)
df_v$NCBI_full_taxon_lineage <- sub(";[^;]*;", ";", df_v$NCBI_full_taxon_lineage)
df_v$NCBI_full_taxon_lineage <- sub(";[^;]*;", ";", df_v$NCBI_full_taxon_lineage)
df_v$NCBI_full_taxon_lineage <- sub("^[^;]*;", "", df_v$NCBI_full_taxon_lineage)

#### save full taxa lineages ########
df_v_tax <- subset(df_v, select = c(orig_NCBI_full_taxon_lineage, NCBI_full_taxon_lineage))
df_v_tax <- df_v_tax %>% distinct()
write.csv(df_v_tax, "bioinf_readcheck_data_viruses_taxalineages.csv", row.names=FALSE)

## add n of reads per taxa
df_v <- transform(df_v, n_per_taxa = ave(raw_bctrimmed_reads, NCBI_full_taxon_lineage, FUN = sum))

## add n of read per taxa AND assembly steps
df_v <- transform(df_v, n_per_taxa_asssemblystep = ave(raw_bctrimmed_reads, NCBI_full_taxon_lineage, canu_correctedreads, canu_trimmedreads, canu_contigs, FUN = sum))

## remove duplicate information
df_v <- subset(df_v, select = c(NCBI_full_taxon_lineage, canu_correctedreads, canu_trimmedreads, canu_contigs, n_per_taxa, n_per_taxa_asssemblystep))
df_v <- df_v %>% distinct()

## add percent of reads per taxa group AND assembly steps
df_v <- transform(df_v, percent_grouped_by_taxa_assemblysteps = ave(n_per_taxa_asssemblystep, FUN = prop.table))
df_v$percent_grouped_by_taxa_assemblysteps <- specify_decimal(100*(df_v$percent_grouped_by_taxa_assemblysteps), 2)

## add percent of reads per taxa group
df_v <- transform(df_v, sumpercentbytaxa = ave(percent_grouped_by_taxa_assemblysteps, NCBI_full_taxon_lineage, FUN = sum))

#### output table for virus data
write.csv(df_v, "bioinf_readcheck_data_viruses.csv", row.names=FALSE)

### calcuate percent of reads per assembly step, within each taxa
df_v$percent_perassemblystep_withintaxa <- (df_v$n_per_taxa_asssemblystep/df_v$n_per_taxa)*100


## remove all taxa with less than 10 reads
df_v <- as.data.frame(filter(df_v, n_per_taxa > 10))

#######################
#######################

df_a <- df

### reduce taxa levels (repeat to remove next lowest level)
df_a$NCBI_full_taxon_lineage <- sub(";[^;]*;$", ";", df_a$NCBI_full_taxon_lineage)
df_a$NCBI_full_taxon_lineage <- sub(";[^;]*;$", ";", df_a$NCBI_full_taxon_lineage)
df_a$NCBI_full_taxon_lineage <- sub(";[^;]*;$", ";", df_a$NCBI_full_taxon_lineage)
df_a$NCBI_full_taxon_lineage <- sub(";[^;]*;$", ";", df_a$NCBI_full_taxon_lineage)
df_a$NCBI_full_taxon_lineage <- sub(";[^;]*;$", ";", df_a$NCBI_full_taxon_lineage)
df_a$NCBI_full_taxon_lineage <- sub(";$", "", df_a$NCBI_full_taxon_lineage)

#### save full taxa lineages ########
df_a_tax <- subset(df_a, select = c(orig_NCBI_full_taxon_lineage, NCBI_full_taxon_lineage))
df_a_tax <- df_a_tax %>% distinct()
write.csv(df_a_tax, "bioinf_readcheck_data_all_taxalineages.csv", row.names=FALSE)

## add n of reads per taxa
df_a <- transform(df_a, n_per_taxa = ave(raw_bctrimmed_reads, NCBI_full_taxon_lineage, FUN = sum))

## add n of read per taxa AND assembly steps
df_a <- transform(df_a, n_per_taxa_asssemblystep = ave(raw_bctrimmed_reads, NCBI_full_taxon_lineage, canu_correctedreads, canu_trimmedreads, canu_contigs, FUN = sum))

## remove duplicate information
df_a <- subset(df_a, select = c(NCBI_full_taxon_lineage, canu_correctedreads, canu_trimmedreads, canu_contigs, n_per_taxa, n_per_taxa_asssemblystep))
df_a <- df_a %>% distinct()


## add percent of reads per taxa group AND assembly steps
df_a <- transform(df_a, percent_grouped_by_taxa_assemblysteps = ave(n_per_taxa_asssemblystep, FUN = prop.table))
df_a$percent_grouped_by_taxa_assemblysteps <- specify_decimal(100*(df_a$percent_grouped_by_taxa_assemblysteps), 2)

## add percent of reads per taxa group
df_a <- transform(df_a, sumpercentbytaxa = ave(percent_grouped_by_taxa_assemblysteps, NCBI_full_taxon_lineage, FUN = sum))

### calcuate percent of reads per assembly step, within each taxa
df_a$percent_perassemblystep_withintaxa <- (df_a$n_per_taxa_asssemblystep/df_a$n_per_taxa)*100

#### output table for all data
write.csv(df_a, "bioinf_readcheck_data_all.csv", row.names=FALSE)

## remove all taxa with less than 10 reads
df_a <- as.data.frame(filter(df_a, n_per_taxa > 10))

## filter taxonomic groups by percent of raw reads
df_a_lt10 <- filter(df_a, sumpercentbytaxa < 10)
df_a_lt1 <- filter(df_a, sumpercentbytaxa < 1)

###############################################################################
###############################################################################

#### check if df is right format
is_alluvia_form(as.data.frame(df_a), axes = 1:5, silent = TRUE)

num_colors <- length(unique(df_a$NCBI_full_taxon_lineage))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(num_colors)
plottitle = str_wrap("Alluvial plot showing the percent of all raw reads that passed or failed assembly correction, trimming, and contig assembly (left to right). Taxonomic classifications performed on reads by Kaiju in greedy mode against the RVDB and nr_euk databases. Total number of reads per taxa is n. Taxa with less than 10 reads are not shown. AR. = Archaea, BA. = Bacteria, EU. = Eukaryota, VI. = Viruses, UN. = Unclassified, co = correction of reads by canu, tr = trimming of reads by canu, as = assembly into contigs by canu", 75)


  ggplot(as.data.frame(df_a), aes(y = percent_grouped_by_taxa_assemblysteps, axis1 = paste(NCBI_full_taxon_lineage,"(n=", n_per_taxa,")"), axis2 = canu_correctedreads, axis3 = canu_trimmedreads, axis4 = canu_contigs)) +
  geom_alluvium(aes(fill = factor(NCBI_full_taxon_lineage)), width = 1/12, alpha = 0.7) +
  geom_stratum(aes(fill = factor(NCBI_full_taxon_lineage)), width = 1/12, color = "black", alpha = 0.7) +
  geom_label_repel(stat = "stratum", aes(label = after_stat(stratum), fill = factor(NCBI_full_taxon_lineage)), size = 4, max.overlaps = Inf, force_pull = 0.2, force = 0.2, nudge_x = -0.25, direction = "y", hjust = 1.0, segment.size = 0.1, show.legend = FALSE) +
  scale_x_discrete(limits = c("raw barcode trimmed reads", "canu read correction", "canu read trim", "canu contig assembly"), expand = expansion(mult = c(0.35, 0))) +
  scale_fill_manual(values = mycolors) +
  theme(legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=12), axis.text.x = element_text(size=14)) +
  labs(caption = plottitle) +
  theme(plot.caption = element_text(hjust=0.5, size = 16)) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(alpha = 0.7))) +
  coord_cartesian(clip = "off")

ggsave(plot = last_plot(),paste0("bioinf_readcheck_alluvial_all", format(Sys.time(), "_%Y-%m-%d"), ".pdf"), dpi=300, scale=2, units = "cm")

##############################

num_colors <- length(unique(df_v$NCBI_full_taxon_lineage))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(num_colors)
plottitle = str_wrap("Alluvial plot showing the percent of viral raw reads that passed or failed assembly correction, trimming, and contig assembly (left to right). Taxonomic classifications performed on reads by Kaiju in greedy mode against the RVDB and nr_euk databases. Percents are of raw reads that classified as viruses. Total number of reads per taxa is n. Taxa with less than 10 reads are not shown. UN. = Unclassified, co = correction of reads by canu, tr = trimming of reads by canu, as = assembly into contigs by canu", 75)


ggplot(df_v, aes(y = percent_grouped_by_taxa_assemblysteps, axis1 = paste(NCBI_full_taxon_lineage,"(n=", n_per_taxa,")"), axis2 = canu_correctedreads, axis3 = canu_trimmedreads, axis4 = canu_contigs)) +
  geom_alluvium(aes(fill = factor(NCBI_full_taxon_lineage)), width = 1/12, alpha = 0.7) +
  geom_stratum(aes(fill = factor(NCBI_full_taxon_lineage)), width = 1/12, color = "black", alpha = 0.7) +
  geom_label_repel(stat = "stratum", aes(label = after_stat(stratum), fill = factor(NCBI_full_taxon_lineage)), size = 3, max.overlaps = Inf, force_pull = 0.2, force = 0.2, nudge_x = -0.25, direction = "y", hjust = 1.0, segment.size = 0.1, show.legend = FALSE) +
  scale_x_discrete(limits = c("raw barcode trimmed reads", "canu read correction", "canu read trim", "canu contig assembly"), expand = expansion(mult = c(0.7, 0))) +
  scale_fill_manual(values = mycolors) +
  theme(legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=12), axis.text.x = element_text(size=14)) +
  labs(caption = plottitle) +
  theme(plot.caption = element_text(hjust=0.5, size = 16)) +
  guides(fill=guide_legend(ncol =1, override.aes = list(alpha = 0.7))) +
  coord_cartesian(clip = "off")

ggsave(plot = last_plot(),paste0("bioinf_readcheck_alluvial_viruses", format(Sys.time(), "_%Y-%m-%d"), ".pdf"), dpi=300, scale=2, units = "cm")

##############################

num_colors <- length(unique(df_a_lt10$NCBI_full_taxon_lineage))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(num_colors)
plottitle = str_wrap("Alluvial plot showing the percent of all raw reads that passed or failed assembly correction, trimming, and contig assembly (left to right). Taxonomic classifications performed on reads by Kaiju in greedy mode against the RVDB and nr_euk databases. Only taxa comprising 10 percent or less of reads each are shown. Total number of reads per taxa is n. Taxa with less than 10 reads are not shown. AR. = Archaea, BA. = Bacteria, EU. = Eukaryota, VI. = Viruses, UN. = Unclassified, co = correction of reads by canu, tr = trimming of reads by canu, as = assembly into contigs by canu", 75)

ggplot(as.data.frame(df_a_lt10), aes(y = percent_grouped_by_taxa_assemblysteps, axis1 = paste(NCBI_full_taxon_lineage,"(n=", n_per_taxa,")"), axis2 = canu_correctedreads, axis3 = canu_trimmedreads, axis4 = canu_contigs)) +
  geom_alluvium(aes(fill = factor(NCBI_full_taxon_lineage)), width = 1/12, alpha = 0.7) +
  geom_stratum(aes(fill = factor(NCBI_full_taxon_lineage)), width = 1/12, color = "black", alpha = 0.7) +
  geom_label_repel(stat = "stratum", aes(label = after_stat(stratum), fill = factor(NCBI_full_taxon_lineage)), size = 4, max.overlaps = Inf, force_pull = 0.2, force = 0.2, nudge_x = -0.25, direction = "y", hjust = 1.0, segment.size = 0.1, show.legend = FALSE) +
  scale_x_discrete(limits = c("raw barcode trimmed reads", "canu read correction", "canu read trim", "canu contig assembly"), expand = expansion(mult = c(0.35, 0))) +
  scale_fill_manual(values = mycolors) +
  theme(legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=12), axis.text.x = element_text(size=14)) +
  labs(caption = plottitle) +
  theme(plot.caption = element_text(hjust=0.5, size = 16)) +
  guides(fill=guide_legend(ncol =1, override.aes = list(alpha = 0.7))) +
  coord_cartesian(clip = "off")

ggsave(plot = last_plot(),paste0("bioinf_readcheck_alluvial_lt10percent", format(Sys.time(), "_%Y-%m-%d"), ".pdf"), dpi=300, scale=2, units = "cm")

#############################

num_colors <- length(unique(df_a_lt1$NCBI_full_taxon_lineage))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(num_colors)
plottitle = str_wrap("Alluvial plot showing the percent of all raw reads that passed or failed assembly correction, trimming, and contig assembly (left to right). Taxonomic classifications performed on reads by Kaiju in greedy mode against the RVDB and nr_euk databases. Only taxa comprising 1 percent or less of reads each are shown. Total number of reads per taxa is n. Taxa with less than 10 reads are not shown. AR. = Archaea, BA. = Bacteria, EU. = Eukaryota, VI. = Viruses, UN. = Unclassified, co = correction of reads by canu, tr = trimming of reads by canu, as = assembly into contigs by canu", 75)

  ggplot(as.data.frame(df_a_lt1), aes(y = percent_grouped_by_taxa_assemblysteps, axis1 = paste(NCBI_full_taxon_lineage,"(n=", n_per_taxa,")"), axis2 = canu_correctedreads, axis3 = canu_trimmedreads, axis4 = canu_contigs)) +
  geom_alluvium(aes(fill = factor(NCBI_full_taxon_lineage)), width = 1/12, alpha = 0.7) +
  geom_stratum(aes(fill = factor(NCBI_full_taxon_lineage)), width = 1/12, color = "black", alpha = 0.7) +
  geom_label_repel(stat = "stratum", aes(label = after_stat(stratum), fill = factor(NCBI_full_taxon_lineage)), size = 4, max.overlaps = Inf, force_pull = 0.2, force = 0.2, nudge_x = -0.25, direction = "y", hjust = 1.0, segment.size = 0.1, show.legend = FALSE) +
  scale_x_discrete(limits = c("raw barcode trimmed reads", "canu read correction", "canu read trim", "canu contig assembly"), expand = expansion(mult = c(0.35, 0))) +
  scale_fill_manual(values = mycolors) +
  theme(legend.key.size = unit(0.5, 'cm'), legend.text = element_text(size=12), axis.text.x = element_text(size=14)) +
  labs(caption = plottitle) +
  theme(plot.caption = element_text(hjust=0.5, size = 16)) +
  guides(fill=guide_legend(ncol =1, override.aes = list(alpha = 0.7))) +
  coord_cartesian(clip = "off")

ggsave(plot = last_plot(),paste0("bioinf_readcheck_alluvial_lt1percent", format(Sys.time(), "_%Y-%m-%d"), ".pdf"), dpi=300, scale=2, units = "cm")



###### add extra cols for treemap and barplot


#############################
############################

###### treemap showing read classifications

### for all data

num_colors <- length(unique(df_a$NCBI_full_taxon_lineage))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(num_colors)
plottitle = str_wrap("Treemap of read taxa composition. Taxonomic classifications performed on reads by Kaiju in greedy mode against the RVDB and nr_euk databases. For each taxonomic classification, percentages are of all raw reads and n is the number of total reads. Taxa with less than 10 reads are not shown. AR. = Archaea, BA. = Bacteria, EU. = Eukaryota, VI. = Viruses, UN. = Unclassified", 75)

ggplot(data=df_a[!duplicated(df_a$NCBI_full_taxon_lineage),], aes(area = sumpercentbytaxa, fill = NCBI_full_taxon_lineage, label = paste(NCBI_full_taxon_lineage, paste("n = ", n_per_taxa), paste(sumpercentbytaxa, "%", sep = " "), sep = "\n"))) + 
geom_treemap() + geom_treemap_text(colour = "white", place = "centre", size = 15, reflow = TRUE) +
scale_fill_manual(values = mycolors) +
theme(plot.caption = element_text(hjust=0.5, size = 16)) +
guides(fill=guide_legend(ncol =1, override.aes = list(alpha = 1.0))) +
coord_cartesian(clip = "off") +
labs(caption = plottitle)

ggsave(plot = last_plot(),paste0("bioinf_readcheck_treemap_all", format(Sys.time(), "_%Y-%m-%d"), ".pdf"), dpi=300, scale=2, units = "cm")

#### for all viruses

num_colors <- length(unique(df_v$NCBI_full_taxon_lineage))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(num_colors)
plottitle = str_wrap("Treemap of read viral taxa composition. Taxonomic classifications performed on reads by Kaiju in greedy mode against the RVDB and nr_euk databases. For each taxonomic classification, percentages are of viral raw reads and n is the number of total reads. Taxa with less than 10 reads are not shown. VI. = Viruses, UN. = Unclassified", 75)

ggplot(data=df_v[!duplicated(df_v$NCBI_full_taxon_lineage),], aes(area = sumpercentbytaxa, fill = NCBI_full_taxon_lineage, label = paste(NCBI_full_taxon_lineage, paste("n = ", n_per_taxa), paste(sumpercentbytaxa, "%", sep = " "), sep = "\n"))) + 
geom_treemap() + geom_treemap_text(colour = "white", place = "centre", size = 15, reflow = TRUE) +
scale_fill_manual(values = mycolors) +
theme(plot.caption = element_text(hjust=0.5, size = 16)) +
guides(fill=guide_legend(ncol =1, override.aes = list(alpha = 1.0))) +
coord_cartesian(clip = "off") +
labs(caption = plottitle)

ggsave(plot = last_plot(),paste0("bioinf_readcheck_treemap_viruses", format(Sys.time(), "_%Y-%m-%d"), ".pdf"), dpi=300, scale=2, units = "cm")

##### barplot of read classifications, plus categories

########## for all data

plottitle = str_wrap("Percentages of reads for each taxa and by assembly category. During assembly, Canu processes the reads with a correction step (co), then a trimming step (tr), and then finally contig assembly (as). For each taxa, the percentages of reads that passed each step are indicated by the legend. Taxonomic classifications performed on reads by Kaiju in greedy mode against the RVDB and nr_euk databases. The total number of reads for each taxa is n. Taxa with less than 10 reads are not shown. AR. = Archaea, BA. = Bacteria, EU. = Eukaryota, VI. = Viruses, UN. = Unclassified, co = correction of reads by canu, tr = trimming of reads by canu, as = assembly into contigs by canu", 75)

ggplot(df_a, aes(fill=paste(canu_correctedreads, canu_trimmedreads, canu_contigs, sep=" "), y=percent_perassemblystep_withintaxa, x=paste(NCBI_full_taxon_lineage,"(n=", n_per_taxa,")"))) +
geom_bar(position="stack", stat="identity") +
theme(axis.text.x = element_text(angle = 90, size = 14, vjust = 0.5)) +
theme(axis.text.y = element_text(size = 12)) +
labs(caption = plottitle) +
guides(fill=guide_legend(title="Assembly steps")) +
theme(plot.caption = element_text(hjust=0.5, size = 14, vjust = 0.5)) +
scale_fill_brewer(palette = "Spectral") +
coord_flip()

ggsave(plot = last_plot(),paste0("bioinf_readcheck_barplot_all", format(Sys.time(), "_%Y-%m-%d"), ".pdf"), dpi=300, scale=2, units = "cm")

file.remove("Rplots.pdf")