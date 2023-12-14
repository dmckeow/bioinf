#!/usr/bin/env Rscript

# Script use: 


library(dplyr)
library(tidyverse)
library(ggplot2)
#library(ggalluvial)
library(RColorBrewer)
library(svglite)
#library(ggrepel)
#library(treemapify)
library(metagenomeSeq)
#library(pheatmap)
library(googlesheets4)
library(viridis)
#library(UpSetR)
library(ggridges)
library(ggh4x)
#library(vegan)
#library(forcats)
library(hexbin)

#library(reshape) #
#library(plotrix) #
#library(tidyr) #
#library(viridisLite) #


theme_set(theme_classic())
set.seed(1234)


################################################
############## IMPORT DATASHEETS ###############
################################################

df <- read.csv("ALL.bam.reads_mapped.ALL", header=TRUE, sep="\t")

################################################
############## WRANGLE DATASHEETS ##############
################################################

#### set a function to use for rounding decimal places
specify_decimal <- function(x, k) as.numeric(trimws(format(round(x, k), nsmall=k)))

## store original taxa lineage
df$orig_kaiju_taxonomy <- df$kaiju_taxonomy

df$kaiju_taxonomy <- gsub("; ", ";", df$kaiju_taxonomy)
df$kaiju_taxonomy <- gsub(", ", ";", df$kaiju_taxonomy)
df$kaiju_taxonomy <- gsub(",", ";", df$kaiju_taxonomy)
df$kaiju_taxonomy <- gsub("[^A-Za-z0-9_;.]", "_", df$kaiju_taxonomy)
df$kaiju_taxonomy <- gsub("_+", "_", df$kaiju_taxonomy)
df$kaiju_taxonomy <- sub("$", ";", df$kaiju_taxonomy)
df$kaiju_taxonomy <- gsub("NA;", ";", df$kaiju_taxonomy)
df$kaiju_taxonomy <- gsub(";+", ";", df$kaiju_taxonomy)

df$kaiju_taxonomy <- sub("Archaea;", "AR;", df$kaiju_taxonomy)
df$kaiju_taxonomy <- sub("Bacteria;", "BA.;", df$kaiju_taxonomy)
df$kaiju_taxonomy <- sub("Eukaryota;", "EU.;", df$kaiju_taxonomy)
df$kaiju_taxonomy <- sub("Viruses;", "VI.;", df$kaiju_taxonomy)

df <- within(df, kaiju_taxonomy[kaiju_taxonomy == ';'] <- 'UN.')
df <- within(df, kaiju_taxonomy[kaiju_taxonomy == ''] <- 'UN.')
df <- within(df, kaiju_taxonomy[kaiju_taxonomy == 'NA'] <- 'UN.')
df <- df %>% mutate(kaiju_taxonomy = ifelse(is.na(kaiju_taxonomy), 'UN.', kaiju_taxonomy))
df$kaiju_taxonomy <- sub(";$", "", df$kaiju_taxonomy)
df$kaiju_taxonomy <- sub("^VI.$", "UN.", df$kaiju_taxonomy)

##############################

###### specific taxa lineage collapsing
### reduce Bacteria to Bacteria
df$kaiju_taxonomy <- gsub("^BA..*", "BA.", df$kaiju_taxonomy)
df$kaiju_taxonomy <- gsub(".*Arthropoda.*", "Arthropoda", df$kaiju_taxonomy) ## anything that is arthropod to just arthropod
df$kaiju_taxonomy <- gsub("^EU..*", "EU.", df$kaiju_taxonomy) ## anything else Eukaryotic, but not Arthropod to just Eukaryotic

###df <- df %>% filter(grepl('^VI.;', kaiju_taxonomy))

### reduce taxa levels (repeat to remove next lowest level)
df$kaiju_taxonomy <- sub(";[^;]*;", ";", df$kaiju_taxonomy)
df$kaiju_taxonomy <- sub(";[^;]*;", ";", df$kaiju_taxonomy)
df$kaiju_taxonomy <- sub(";[^;]*;", ";", df$kaiju_taxonomy)
###df$kaiju_taxonomy <- sub("^[^;]*;", "", df$kaiju_taxonomy)

####### make higher taxa groupings
df$kaiju_taxonomy_higher <- df$kaiju_taxonomy
df$kaiju_taxonomy_higher <- gsub("^BA.$", "non_viral", df$kaiju_taxonomy_higher)
df$kaiju_taxonomy_higher <- gsub("^EU.$", "non_viral", df$kaiju_taxonomy_higher)
df$kaiju_taxonomy_higher <- gsub("^Arthropoda$", "non_viral", df$kaiju_taxonomy_higher)
df$kaiju_taxonomy_higher <- gsub("^UN.$", "non_viral", df$kaiju_taxonomy_higher)
df$kaiju_taxonomy_higher <- gsub("^VI.;", "", df$kaiju_taxonomy_higher)
df$kaiju_taxonomy_higher <- gsub(";.*", "", df$kaiju_taxonomy_higher)

df$kaiju_taxonomy <- gsub("^VI.;", "", df$kaiju_taxonomy)
df$kaiju_taxonomy <- sub("^[^;]*;", "", df$kaiju_taxonomy)

########### other format fixes
df$specific_read_source <- gsub(".bctrimmedreads.fastq.gz", "", df$specific_read_source)
df$specific_contig_source <- gsub(".contigs.fasta.*", "", df$specific_contig)


## add n of read per taxa PER sample
df <- transform(df, n_per_taxa_per_sample = ave(num_specific_read_source_reads_mapped_to_specific_contig, kaiju_taxonomy, specific_read_source, FUN = sum))

## remove duplicate information (for reads)
df_r <- subset(df, select = c(n_per_taxa_per_sample, specific_read_source, kaiju_taxonomy))
df_r <- df_r %>% distinct()


######################################################
############# contig composition

df_con <- df %>% select(-specific_read_source, -num_specific_read_source_reads_mapped_to_specific_contig, -total_num_reads_mapped_to_all_contigs_from_specific_read_source, -num_reads_in_specific_read_source, -n_per_taxa_per_sample) ## remove taxa we don't want in plot
df_con <- df_con %>% distinct()

## add n of contig per taxa PER sample
df_con$one <- as.numeric(gsub(".*", "1", df_con$specific_contig_source))

df_con <- transform(df_con, n_contigs_per_taxa_per_sample = ave(one, kaiju_taxonomy, specific_contig_source, FUN = sum))
df_con$specific_contig_source <- gsub("^", "X", df_con$specific_contig_source)

df_con_counts <- subset(df_con, select = c(n_contigs_per_taxa_per_sample, specific_contig_source, kaiju_taxonomy))

df_con_counts <- df_con_counts %>% distinct()


######################################################
######################################################

####################
#### get metadata for spillover and merge with other info
dfmd <- read_sheet("https://docs.google.com/spreadsheets/d/1yDaJm30o-FcLEIP2iyX3JHZAzWVvFCbxjxdfXqnR73w/edit?pli=1#gid=0")
dfmd <- lapply(dfmd, as.character)
dfmd <- data.frame(dfmd)
dfmd$specific_read_source <- gsub("^", "X", dfmd$specific_read_source)
dfmd <- dfmd %>% select(specific_read_source, genus, species, flower_genus, flower_species, collection_month, collection_year, apiary, distance)
dfmd <- dfmd %>% column_to_rownames("specific_read_source")
dfmd <- data.frame(dfmd)
###########################
###########################

###############################
#### average by genus, species, collection_month, collection_year, apiary, distance
df_r_av <- data.frame(df_r)
df_r_av$specific_read_source <- gsub("^", "X", df_r_av$specific_read_source)
df_r_av <- merge(df_r_av, dfmd, by.x="specific_read_source", by.y=0)

## need to remove totals of 0, because they are not appropriate to include when averaging per taxa?
###df_r_av <- subset(df_r_av, df_r_av$n_per_taxa_per_sample > 0)

df_r_av <- transform(df_r_av, n_per_taxa_per_sample = ave(n_per_taxa_per_sample, kaiju_taxonomy, genus, species, collection_month, collection_year, apiary, distance, FUN = ave))

df_r_av <- df_r_av %>% group_by(n_per_taxa_per_sample, kaiju_taxonomy, genus, species, collection_month, collection_year, apiary, distance) %>% summarise(specific_read_source = paste(specific_read_source, collapse = ";"))
df_r_av <- df_r_av %>% ungroup() %>% select(n_per_taxa_per_sample, specific_read_source, kaiju_taxonomy) ## remove taxa we don't want in plot
df_r_av <- data.frame(df_r_av)


###############################
###############################

###############################
### df_r for per sample, df_r_av for average by sample groups

colnames(df_r) <- c("count", "sample", "OTU")
count_matrix <- df_r %>% pivot_wider(names_from = sample, values_from = count)
##### OR ####
#colnames(df_r_av) <- c("count", "sample", "OTU")
#count_matrix <- df_r_av %>% pivot_wider(names_from = sample, values_from = count)

count_matrix <- as.data.frame(count_matrix)
count_matrix[is.na(count_matrix)] <- 0
write.csv(count_matrix, "count_matrix.csv", row.names = FALSE)
count_matrix_obj = loadMeta("count_matrix.csv", sep = ",")

taxa = subset(count_matrix, select = c(OTU, OTU))
taxa <- as.data.frame(taxa)
write.csv(taxa, "taxa.csv", row.names = FALSE)

taxa = read.delim("taxa.csv", sep=",", stringsAsFactors = FALSE, row.names = 1)

taxa = AnnotatedDataFrame(taxa)
obj = newMRexperiment(count_matrix_obj$counts, featureData = taxa)
obj_CSS = cumNorm(obj, p=cumNormStat(obj))
taxa_obj_CSS = data.frame(MRcounts(obj_CSS, norm=TRUE, log=TRUE))
taxa_obj_CSS <- t(taxa_obj_CSS)


########## ggplot

############################
### identify which samples had contigs for a taxa
###########################
###########################

### pivot table back to long format
taxa_obj_CSS <- taxa_obj_CSS %>% as.data.frame()
taxa_obj_CSS_HM <- taxa_obj_CSS %>% rownames_to_column("samples") %>% pivot_longer(-c(samples), names_to = "taxa", values_to = "counts")

taxa_obj_CSS_HM <- merge(taxa_obj_CSS_HM, dfmd, by.x="samples", by.y=0)

##### add back all sample info and metadata
df_taxa <- df %>% select(kaiju_taxonomy, kaiju_taxonomy_higher)
df_taxa <- df_taxa %>% distinct()
taxa_obj_CSS_HM <- merge(taxa_obj_CSS_HM, df_taxa, by.x="taxa", by.y="kaiju_taxonomy")

taxa_obj_CSS_HM$distance <- as.character(taxa_obj_CSS_HM$distance)

### rename columns
names(taxa_obj_CSS_HM)[names(taxa_obj_CSS_HM) == 'samples'] <- 'specific_read_source'
names(taxa_obj_CSS_HM)[names(taxa_obj_CSS_HM) == 'taxa'] <- 'kaiju_taxonomy'

###################################################

#### add contig presence or absence per sample by taxa
taxa_obj_CSS_HM_c <- merge(df_con_counts, taxa_obj_CSS_HM, by.x=c("specific_contig_source","kaiju_taxonomy"), by.y=c("specific_read_source","kaiju_taxonomy"), all.y=TRUE)
taxa_obj_CSS_HM_c$n_contigs_per_taxa_per_sample <- gsub("[0-9]+", "1", taxa_obj_CSS_HM_c$n_contigs_per_taxa_per_sample)


###################################################
###################################################

###plot_coverage(df_merged, percent_of_contig_len_with_over_5X_coverage)
taxa_obj_CSS_HM_c %>% ggplot(aes(x=specific_contig_source, y=kaiju_taxonomy, fill=counts)) + 
  geom_tile() +
  scale_fill_viridis(option="inferno") +
  facet_nested(kaiju_taxonomy_higher + kaiju_taxonomy ~ genus, switch = "both", space = "free", scales = "free") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y=element_blank()) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1)) +
theme(strip.background.y = element_blank()) +
theme(panel.border = element_blank()) +
theme(axis.title = element_blank()) +
theme(panel.spacing.x=unit(0.2, "lines")) +
theme(panel.spacing.y=unit(0.1, "lines")) +
theme(strip.placement = "outside") +
theme(ggh4x.facet.nestline = element_line(colour = "black")) +
geom_point(aes(size = n_contigs_per_taxa_per_sample), color = "grey", show.legend = FALSE)


###############################
###############################

############# breadth of coverage

## import data from bash script

con <- read.csv("tmp.ALLcontig", header=TRUE, sep="\t")
sam_con <- read.csv("tmp.ALLstats.list.count0", header=TRUE, sep="\t")
bas_1 <- read.csv("tmp.ALLstats.list.count1", header=TRUE, sep="\t")
bas_2 <- read.csv("tmp.ALLstats.list.count2", header=TRUE, sep="\t")
bas_3 <- read.csv("tmp.ALLstats.list.count3", header=TRUE, sep="\t")


df_merged <- merge(con, sam_con, by.x="contig", by.y="contig")
df_merged <- merge(bas_1, df_merged, by.x=c("reads","contig"), by.y=c("reads","contig"), all.y=TRUE)
df_merged <- merge(bas_2, df_merged, by.x=c("reads","contig"), by.y=c("reads","contig"), all.y=TRUE)
df_merged <- merge(bas_3, df_merged, by.x=c("reads","contig"), by.y=c("reads","contig"), all.y=TRUE)

## add average coverage depth per base for each contig (total coveage per base for contig divided by contig length), then times 100 to get the percent of the contig length in bases that passed the coverage threshold in each category - e.g. 50 % of contig had at least 5X coverage
df_merged <- transform(df_merged, percent_of_contig_len_with_over_5X_coverage = total_bases_over_5_coverage / contig_length * 100)
df_merged <- transform(df_merged, percent_of_contig_len_with_over_10X_coverage = total_bases_over_10_coverage / contig_length * 100)
df_merged <- transform(df_merged, percent_of_contig_len_with_over_100X_coverage = total_bases_over_100_coverage / contig_length * 100)

## get the contig classifications info
df_r <- subset(df, select = c(specific_contig, kaiju_taxonomy, kaiju_taxonomy_higher)) ## select only some columns
df_r <- df_r %>% distinct() ## remove duplicate information

########### other format fixes
df_merged$specific_read_source <- gsub(".bctrimmedreads.fastq.gz", "", df_merged$reads)
df_merged$specific_read_source <- gsub("^", "X", df_merged$specific_read_source)

## merge meta data and data
df_merged <- merge(df_merged, df_r, by.x="contig", by.y="specific_contig") # with the contig classifications and the coverage info
df_merged <- merge(df_merged, dfmd, by.x="specific_read_source", by.y=0) # with the sampel metadata

### Reshape data by gathering data to make coverage a category

df_merged_cat <- gather(df_merged, key = "cov_category", value = "percent_of_contig_len_with_over_nX_coverage", percent_of_contig_len_with_over_5X_coverage, percent_of_contig_len_with_over_10X_coverage, percent_of_contig_len_with_over_100X_coverage)

df_merged_cat$cov_category <- gsub("percent_of_contig_len_with_over_5X_coverage", "> 5x coverage", df_merged_cat$cov_category)
df_merged_cat$cov_category <- gsub("percent_of_contig_len_with_over_10X_coverage", "> 10x coverage", df_merged_cat$cov_category)
df_merged_cat$cov_category <- gsub("percent_of_contig_len_with_over_100X_coverage", "> 100x coverage", df_merged_cat$cov_category)


# Arrange the dataframe by the 'group' column and 'Name' column
df_merged_cat <- df_merged_cat %>% arrange(kaiju_taxonomy_higher, kaiju_taxonomy)
# Set the order of the 'Name' column based on the alphabetical order of 'group'
df_merged_cat$kaiju_taxonomy <- factor(df_merged_cat$kaiju_taxonomy, levels = unique(df_merged_cat$kaiju_taxonomy))

### scatter plot showing the coverage breadth by threshold Apis vs Bombus

#plot_coverage <- function(data, xaxis){
#data %>% filter(!grepl('non_viral', kaiju_taxonomy_higher)) %>% ggplot(aes(x={{xaxis}}, y=kaiju_taxonomy, size=contig_length)) +
#geom_point(alpha=0.5) + 
#facet_grid(kaiju_taxonomy_higher ~ genus, space = "free", scales = "free")
#}

###plot_coverage(df_merged, percent_of_contig_len_with_over_5X_coverage)
###plot_coverage(df_merged, percent_of_contig_len_with_over_10X_coverage)
###plot_coverage(df_merged, percent_of_contig_len_with_over_100X_coverage)

########## as categoryised coverage thresholds SCATTER
df_merged_cat %>% 
filter(!grepl('non_viral', kaiju_taxonomy_higher)) %>% 
#filter(!grepl('> 5x coverage', cov_category)) %>% 
ggplot(aes(x=percent_of_contig_len_with_over_nX_coverage, y=kaiju_taxonomy)) +
#geom_point(size = 3, shape = "|") + 
geom_jitter(aes(color = contig_length), alpha = 1, shape = "|", size = 5, width = 0.2, height = 0.0) + 
#facet_nested(kaiju_taxonomy_higher + kaiju_taxonomy ~ fct_rev(cov_category) + genus, switch = "y") +
facet_nested(kaiju_taxonomy_higher ~ fct_rev(cov_category) + genus, switch = "y", scales = "free", space = "free") +
theme(strip.text.y.left = element_text(angle = 0, hjust = 1)) +
theme(strip.background = element_blank()) +
theme(panel.border = element_blank()) +
theme(axis.title.y = element_blank()) +
theme(panel.grid.major = element_blank()) +
theme(panel.spacing.x=unit(0.2, "lines")) +
theme(panel.spacing.y=unit(0.2, "lines")) +
theme(panel.background = element_rect(fill = NA, color = "black")) +
theme(strip.placement = "outside") +
theme(ggh4x.facet.nestline = element_line(colour = "black")) +
scale_color_binned(type="viridis", option="inferno", n.breaks = 6)





########## as categoryised coverage thresholds RIDGES
#df_merged_cat %>% filter(!grepl('non_viral', kaiju_taxonomy_higher)) %>% ggplot(aes(x=percent_of_contig_len_with_over_nX_coverage, y=genus, fill=genus, height = stat(density))) +
#geom_density_ridges(alpha=0.7, stat = "density", scale = 5) + 
#geom_density_ridges(stat = "binline", bins = 100, scale = 5, draw_baseline = FALSE) +
#facet_nested(kaiju_taxonomy_higher + kaiju_taxonomy ~ fct_rev(cov_category), switch = "y") +
#theme(strip.text.y.left = element_text(angle = 0, hjust = 1)) +
#theme(strip.background = element_blank()) +
#theme(panel.border = element_blank()) +
#theme(axis.text.y=element_blank(), axis.title.y = element_blank()) +
#theme(panel.grid.major = element_blank()) +
#guides(color = guide_legend(override.aes=list(shape = 15, size = 8), reverse=TRUE)) +
#theme(panel.spacing.x=unit(0, "lines")) +
#theme(panel.spacing.y=unit(0.2, "lines")) +
#theme(panel.background = element_rect(fill = NA, color = "black")) +
#theme(strip.placement = "outside") +
#theme(ggh4x.facet.nestline = element_line(colour = "black"))




######### UpSet Plots

#### contig classification by various factors
#df_r$specific_contig_source <- gsub(".contigs.fasta.*", "", df_r$specific_contig)
#df_r$specific_contig_source <- gsub("^", "X", df_r$specific_contig_source)
#df_r$one <- as.numeric(gsub(".*", "1", df_r$specific_contig_source))

#dfm_up <- merge(dfmd, df_r, by.x=0, by.y="specific_contig_source") # with the sampel metadata
#dfm_up <- dfm_up %>% filter(!grepl('non_viral', kaiju_taxonomy_higher))
#dfm_up <- dfm_up %>% rename(specific_contig_source = Row.names)

#dfm_upM <- dfm_up %>% pivot_wider(id_cols=specific_contig, names_from = kaiju_taxonomy, values_from = one) ## by sample 

#dfm_upM <- dfm_upM %>% as.data.frame()
#dfm_upM[is.na(dfm_upM)] <- 0

#dfm_upM$specific_contig <- gsub(".contigs.fasta.*", "", dfm_upM$specific_contig)

#dfm_upM <- dfm_upM %>% group_by(specific_contig) %>% summarise(across(everything(), ~ ifelse(sum(.) > 0, 1, 0)))

#dfm_upM <- dfm_upM %>% rename(specific_contig_source = specific_contig)
#dfm_upM$specific_contig_source <- gsub("^", "X", dfm_upM$specific_contig_source)


#dfm_upM <- merge(dfm_upM, dfm_up, by.x="specific_contig_source", by.y="specific_contig_source") # with the sampel metadata

### remove redundant data FOR contig-only plot
#dfm_upM <- dfm_upM %>% select(-one)

#dfm_upM <- dfm_upM %>% distinct() ## remove duplicate information

#upset(dfm_upM, order.by = "freq")

### viral classification vs genus & species (Apis mellifera, Bombus impatiens, etc), flower_genus (many), apiary (Golf, Bee_Vet)
######## genus, species, flower_genus, flower_species, collection_month, apiary, distance
###dfm_upM <- dfm_upM %>% mutate(genus_species = paste(genus, species, sep = "_"))
#dfm_upM$flower_genus <- gsub("^nd$", "unknown_flower_genus", dfm_upM$flower_genus)
#dfm_upM$distance <- gsub("^nd$", "unknown_distance", dfm_upM$distance)

###dfm_upM <- dfm_upM %>% pivot_wider(names_from = genus, values_from = genus, values_fn = length, values_fill = 0)

#dfm_upM <- dfm_upM %>% select(-specific_contig_source)

#dfm_upM <- dfm_upM %>% relocate(specific_contig)
#dfm_upM <- dfm_upM %>% as.data.frame()

#upset(dfm_upM, order.by = "freq", nsets=30)


########################################################
########################################################

###### get data into a phyloseq object
## otu_table - any numeric matrix
# taxa_obj_CSS
## sample_data - any dataframe - row names must match otu_table sample names
# dfmd
## tax_table - any character matrix - row names must match otu_table OTU names
# df_taxa

library("phyloseq")
library("microViz")
library("patchwork")
library("ggside")



### reshape data for PCA
df_taxa <- df_taxa %>% filter(kaiju_taxonomy != "BA." & kaiju_taxonomy != "UN." & kaiju_taxonomy != "Arthropoda" & kaiju_taxonomy != "EU.") ## remove taxa we don't want in plot
df_taxa <- df_taxa %>% filter(kaiju_taxonomy_higher != "Bromoviridae" & kaiju_taxonomy_higher != "Secoviridae" & kaiju_taxonomy_higher != "Tymoviridae" & kaiju_taxonomy_higher != "Virgaviridae"  & kaiju_taxonomy_higher != "Alphaflexiviridae") ## remove taxa we don't want in plot

df_taxa$kaiju_taxonomy_lower <- df_taxa$kaiju_taxonomy

df_taxa <- df_taxa %>% column_to_rownames("kaiju_taxonomy")
df_taxa <- as.matrix(df_taxa, mode = "character")

taxa_obj_CSS_PCA <- taxa_obj_CSS %>% select(-BA., -UN., -Arthropoda, -EU.) ## remove taxa we don't want in plot
taxa_obj_CSS_PCA <- taxa_obj_CSS_PCA[apply(taxa_obj_CSS_PCA[,-1], 1, function(x) !all(x==0)),] ## remove rows that total to 0 (except first row)
taxa_obj_CSS_PCA <- t(taxa_obj_CSS_PCA)

taxa_obj_CSS_PCA <- merge(df_taxa, taxa_obj_CSS_PCA, by.x=0, by.y=0) ## keep only taxa kept in the df_taxa
taxa_obj_CSS_PCA <- taxa_obj_CSS_PCA %>% select(-kaiju_taxonomy_higher, -kaiju_taxonomy_lower) ## remove taxa we don't want in plot
taxa_obj_CSS_PCA <- taxa_obj_CSS_PCA %>% column_to_rownames("Row.names")
taxa_obj_CSS_PCA <- t(taxa_obj_CSS_PCA)
taxa_obj_CSS_PCA <- taxa_obj_CSS_PCA[apply(taxa_obj_CSS_PCA[,-1], 1, function(x) !all(x==0)),] ## remove rows that total to 0 (except first row)

taxa_obj_CSS_PCA <- t(taxa_obj_CSS_PCA)
taxa_obj_CSS_PCA <- data.matrix(taxa_obj_CSS_PCA)

OTU = otu_table(taxa_obj_CSS_PCA, taxa_are_rows=TRUE)
TAX = tax_table(df_taxa)
sampledata = sample_data(dfmd)

physeq = phyloseq(OTU, TAX, sampledata, remove_undetected = TRUE)

physeq <- phyloseq_validate(physeq)



### calc for PCoA
PCoA_physeq <- physeq %>%
  tax_transform("identity", rank = "kaiju_taxonomy_lower") %>% # don't transform!
  dist_calc("bray") %>%
  ord_calc("PCoA")

### plot PCoA iris
PCoA_iris <- PCoA_physeq %>% ord_plot_iris(tax_level = "kaiju_taxonomy_higher", ord_plot = "none", anno_colour = "genus")
  

### plot PCoA with marginal plots
PCoA_marg <- PCoA_physeq %>%
  ord_plot(color = "genus", size = 2, alpha = 0.7, shape = "apiary") +
  ggside::geom_xsidedensity(aes(fill = distance), alpha = 0.5) +
  ggside::geom_ysidedensity(aes(fill = distance), alpha = 0.5) +
  ggside::theme_ggside_void()

PCoA_marg / PCoA_iris

###################################################################################
################## contig distribution size

df_con_dfmd <- merge(df_con, dfmd, by.x="specific_contig_source", by.y=0)


########## as categoryised coverage thresholds RIDGES
#df_con_dfmd %>%
#filter(!grepl('non_viral', kaiju_taxonomy_higher)) %>%
#filter(!grepl('Bromoviridae', kaiju_taxonomy_higher)) %>%
#filter(!grepl('Secoviridae', kaiju_taxonomy_higher)) %>%
#filter(!grepl('Tymoviridae', kaiju_taxonomy_higher)) %>%
#filter(!grepl('Virgaviridae', kaiju_taxonomy_higher)) %>%
#filter(!grepl('Alphaflexiviridae', kaiju_taxonomy_higher)) %>%
##ggplot(aes(x=specific_contig_length, y=kaiju_taxonomy, fill=genus, height = stat(density))) +
#ggplot(aes(x=specific_contig_length, y=kaiju_taxonomy, fill=genus)) +
#geom_density_ridges(alpha=0.5, jittered_points = TRUE, position = position_points_jitter(width = 0.2, height = 0),
#    point_shape = '|', point_size = 5, point_alpha = 1) + 
##geom_density_ridges(stat = "binline", bins = 100, scale = 5, draw_baseline = FALSE) +
##facet_nested(kaiju_taxonomy_higher ~ ., switch = "y", scales = "free", space = "free") +
#facet_nested(kaiju_taxonomy_higher ~ genus, switch = "y", scales = "free_y") +
#theme(strip.text.y.left = element_text(angle = 0, hjust = 1)) +
#theme(strip.background = element_blank()) +
#theme(panel.border = element_blank()) +
#theme(axis.title.y = element_blank()) +
#theme(panel.grid.major = element_blank()) +
#guides(color = guide_legend(override.aes=list(shape = 15, size = 8), reverse=TRUE)) +
#theme(panel.spacing.x=unit(0, "lines")) +
#theme(panel.spacing.y=unit(0.2, "lines")) +
#theme(panel.background = element_rect(fill = NA, color = "black")) +
#theme(strip.placement = "outside") +
#theme(ggh4x.facet.nestline = element_line(colour = "black"))

############ BOXPLOT plot for contig length distribution

df_con_dfmd %>%
filter(!grepl('non_viral', kaiju_taxonomy_higher)) %>%
filter(!grepl('Bromoviridae', kaiju_taxonomy_higher)) %>%
filter(!grepl('Secoviridae', kaiju_taxonomy_higher)) %>%
filter(!grepl('Tymoviridae', kaiju_taxonomy_higher)) %>%
filter(!grepl('Virgaviridae', kaiju_taxonomy_higher)) %>%
filter(!grepl('Alphaflexiviridae', kaiju_taxonomy_higher)) %>%
ggplot(aes(x=specific_contig_length, y=genus, fill=genus)) +
geom_boxplot(outlier.shape = NA) + 
facet_nested(kaiju_taxonomy_higher + kaiju_taxonomy ~ ., switch = "y", space = "free_y") +
theme(strip.text.y.left = element_text(angle = 0, hjust = 1)) +
theme(strip.background = element_blank()) +
theme(panel.border = element_blank()) +
theme(axis.title.y = element_blank()) +
theme(panel.grid.major = element_blank()) +
guides(color = guide_legend(override.aes=list(shape = 15, size = 8), reverse=TRUE)) +
theme(panel.spacing.x=unit(0, "lines")) +
theme(panel.spacing.y=unit(0.2, "lines")) +
theme(panel.background = element_rect(fill = NA, color = "black")) +
theme(strip.placement = "outside") +
theme(ggh4x.facet.nestline = element_line(colour = "black")) +
scale_x_continuous(limits = c(0, 15000), breaks = seq(0, 15000, 1000)) +
geom_jitter(color="black", size=0.4, alpha=0.9) +
theme(axis.text.y=element_blank()) +
theme(panel.grid.major.x = element_line(color = "grey", linetype = 1))


###################################################################################
###################################################################################


################# CONTIG coverage across length visualisation vs REF genomes

map_win <- read.csv("ALL.mapping.ref.meanwindowdepth", header=TRUE, sep="\t")
map_ref <- read.csv("ALL.mapping.ref.idxstats", header=TRUE, sep="\t")

map <- merge(map_win, map_ref, by.x=c("specific_read_source", "reference"), by.y=c("specific_read_source", "reference"))

### add metadata
map$specific_read_source <- gsub("^", "X", map$specific_read_source)
map$specific_read_source <- gsub(".bctrimmedreads.fastq.gz.*", "", map$specific_read_source)
map <- merge(map, dfmd, by.x="specific_read_source", by.y=0, all.x=TRUE)

#### calc RPKM
map <- map %>% mutate(RPKM=paste0(avg_coverage_by_window/((ref_length/1e3)*(NumReadsInSpecificReadSource/1e6))))
map$RPKM <- as.numeric(map$RPKM)
map$avg_coverage_by_window <- as.numeric(map$avg_coverage_by_window)
map <- map %>% mutate(RPKM=(round(RPKM, 1)))
map <- map %>% mutate(avg_coverage_by_window=(round(avg_coverage_by_window, 1)))

map <- map %>% mutate(avg_coverage_by_window = log10(avg_coverage_by_window + 1))
map <- map %>% mutate(RPKM_log = log10(RPKM + 1))



#####################################################

#### average by genus, species, collection_month, collection_year, apiary, distance

## need to remove totals of 0, because they are not appropriate to include when averaging per taxa?
###df_r_av <- subset(df_r_av, df_r_av$n_per_taxa_per_sample > 0)

map <- transform(map, avg_coverage_by_window_AVG = ave(avg_coverage_by_window, ref_start, reference, genus, species, collection_month, collection_year, apiary, distance, FUN = ave))

map <- transform(map, RPKM_AVG = ave(RPKM, ref_start, reference, genus, species, collection_month, collection_year, apiary, distance, FUN = ave))

map <- map %>% group_by(ref_start, reference, genus, species, collection_month, collection_year, apiary, distance) %>% mutate(AVG_grouping = paste(specific_read_source, collapse = ";")) %>% ungroup()
map <- data.frame(map)

map <- map %>% mutate(RPKM_AVG=(round(RPKM_AVG, 1)))
map <- map %>% mutate(avg_coverage_by_window_AVG=(round(avg_coverage_by_window_AVG, 1)))

map <- map %>% mutate(avg_coverage_by_window_AVG_log = log10(avg_coverage_by_window_AVG + 1))
map <- map %>% mutate(RPKM_AVG_log = log10(RPKM_AVG + 1))

map <- map %>% mutate(RPKM_AVG_log=(round(RPKM_AVG_log, 1)))
map <- map %>% mutate(avg_coverage_by_window_AVG_log=(round(avg_coverage_by_window_AVG_log, 1)))

map$seqrun <- gsub("__barcode.*", "", map$specific_read_source)

#####################################################
#####################################################

########## as categoryised coverage thresholds
map %>%
filter(grepl('NC_002066.1___Sacbrood_virus', reference)) %>%
distinct(ref_start, reference, genus, species, collection_month, collection_year, apiary, distance, .keep_all = TRUE) %>% 
ggplot() +
geom_line(aes(x=ref_start, y=RPKM_AVG_log, group=AVG_grouping, color=genus))
#geom_area(aes(x=ref_start, y=RPKM_log, group=genus, fill=genus), position = 'stack')


###for jsut as a heatmap
map %>% 
#filter(grepl('LSV', reference)) %>%
distinct(ref_start, reference, genus, species, collection_month, collection_year, apiary, distance, .keep_all = TRUE) %>% 
ggplot(aes(x=ref_start, y=specific_read_source, fill=RPKM_log)) + 
  geom_tile() +
  scale_fill_viridis(option="inferno") +
  facet_nested(genus + seqrun ~ reference, switch = "both", space = "free", scales = "free") +
  theme(axis.text.y=element_blank()) +
  theme(strip.text.y.left = element_text(angle = 0, hjust = 1)) +
theme(strip.background.y = element_blank()) +
theme(panel.border = element_blank()) +
theme(axis.title = element_blank()) +
theme(panel.spacing.x=unit(0.2, "lines")) +
theme(panel.spacing.y=unit(0.1, "lines")) +
theme(strip.placement = "outside") +
theme(ggh4x.facet.nestline = element_line(colour = "black"))