#!/usr/bin/env Rscript

# Script use: 
#### setwd("C:/Users/Dean Mckeown/Downloads/Spillover_FINAL/RProject_Spillover_FINAL")

library(here)
library(tidyverse)
library(RColorBrewer)
library(svglite)
library(viridis)
library(ggh4x)
library(vegan)
library(metagenomeSeq)
library(cowplot)
library(ggbeeswarm)

theme_set(theme_bw())
set.seed(1234)


################################################
############## IMPORT DATASHEETS ###############
################################################

###### Info for mapping reads vs contigs
ALLSamplesIdxstats <- read.csv(here("ALLSamples.idxstats"), header=TRUE, sep="\t")
ReadTotalPerSample <- read.csv(here("read_total_per_sample"), header=TRUE, sep="\t")
AllSamtoolsDepthCovMappingByWindow <- read.csv(here("AllSamtoolsDepthCovMappingByWindow.tsv"), header=TRUE, sep="\t")

SelfSamplesIdxstats <- read.csv(here("SelfSamples.idxstats"), header=TRUE, sep="\t")
SelfSamtoolsDepthCovMappingByWindow <- read.csv(here("SelfSamtoolsDepthCovMappingByWindow.tsv"), header=TRUE, sep="\t")

### Info for Each CONTIG
dfmd <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/12Z5JHjL_PFjUF3MLHlOUOjBD8viFRN7FUVJhgLQQCyQ/edit#gid=564525482")
BinList <- read.csv(here("bin_list"), header=TRUE, sep="\t")
BinMmseq <- read.csv(here("bin_mmseq"), header=TRUE, sep="\t")
SeqsRemovedByClustering <- read.csv(here("SeqsRemovedByClustering"), header=TRUE, sep="\t")
RVDBAllContigsDmndBlastx <- read.csv(here("RVDBv26_AllContigs.dmnd.blastx"), header=TRUE, sep="\t")
VOGAllContigsDmndBlastx <- read.csv(here("vog_AllContigs.dmnd.blastx"), header=TRUE, sep="\t")
AllContigsKaijuNrEuk <- read.csv(here("all-contigs.fa.kaiju.nr_euk.names"), header=TRUE, sep="\t")
AllContigsKaijuRVDB <- read.csv(here("all-contigs.fa.kaiju.rvdb.names"), header=TRUE, sep="\t")
CVQualitySummary <- read.csv(here("quality_summary.tsv"), header=TRUE, sep="\t")
AllContigLengths <- read.csv(here("AllContigLengths.tsv"), header=TRUE, sep="\t")

##### UNUSED DATA:
################# CONTIG coverage across length visualisation vs REF genomes
##df <- read.csv(here("ALL.bam.reads_mapped.ALL"), header=TRUE, sep="\t")
##map_win <- read.csv(here("ALL.mapping.ref.meanwindowdepth"), header=TRUE, sep="\t")
##map_ref <- read.csv(here("ALL.mapping.ref.idxstats"), header=TRUE, sep="\t")
##map_ref_tot <- read.csv(here("ALL.mapping.ref.totalmapped"), header=TRUE, sep="\t")
############### checkv outputs
##CVCompleteGenomes <- read.csv(here("complete_genomes.tsv"), header=TRUE, sep="\t")
##CVCompleteness <- read.csv(here("completeness.tsv"), header=TRUE, sep="\t")
##CVContamination <- read.csv(here("contamination.tsv"), header=TRUE, sep="\t")
##sam_con <- read.csv(here("tmp.ALLstats.list.count0"), header=TRUE, sep="\t")



################################################
############## WRANGLE DATASHEETS ##############
################################################

############### YOU CAN SKIP THIS PART IF YOU ALREADY CURATED dfContigs
###### Clean up kaiju taxonomy lineages
## store original taxa lineage
##Tidy_Taxonomy <- function(dfTidyTaxonomy){
##    dfTidyTaxonomy$orig_kaiju_taxonomy <- dfTidyTaxonomy$kaiju_taxonomy
##
##    dfTidyTaxonomy$kaiju_taxonomy <- gsub("; ", ";", dfTidyTaxonomy$kaiju_taxonomy)
##    dfTidyTaxonomy$kaiju_taxonomy <- gsub(", ", ";", dfTidyTaxonomy$kaiju_taxonomy)
##    dfTidyTaxonomy$kaiju_taxonomy <- gsub(",", ";", dfTidyTaxonomy$kaiju_taxonomy)
##    dfTidyTaxonomy$kaiju_taxonomy <- gsub("[^A-Za-z0-9_;.]", "_", dfTidyTaxonomy$kaiju_taxonomy)
##    dfTidyTaxonomy$kaiju_taxonomy <- gsub("_+", "_", dfTidyTaxonomy$kaiju_taxonomy)
##    dfTidyTaxonomy$kaiju_taxonomy <- sub("$", ";", dfTidyTaxonomy$kaiju_taxonomy)
##    dfTidyTaxonomy$kaiju_taxonomy <- gsub("NA;", ";", dfTidyTaxonomy$kaiju_taxonomy)
##    dfTidyTaxonomy$kaiju_taxonomy <- gsub(";+", ";", dfTidyTaxonomy$kaiju_taxonomy)
##
##    dfTidyTaxonomy$kaiju_taxonomy <- sub("Archaea;", "AR;", dfTidyTaxonomy$kaiju_taxonomy)
##    dfTidyTaxonomy$kaiju_taxonomy <- sub("Bacteria;", "BA.;", dfTidyTaxonomy$kaiju_taxonomy)
##    dfTidyTaxonomy$kaiju_taxonomy <- sub("Eukaryota;", "EU.;", dfTidyTaxonomy$kaiju_taxonomy)
##    dfTidyTaxonomy$kaiju_taxonomy <- sub("Viruses;", "VI.;", dfTidyTaxonomy$kaiju_taxonomy)
##
##    dfTidyTaxonomy <- within(dfTidyTaxonomy, kaiju_taxonomy[kaiju_taxonomy == ';'] <- 'UN.')
##    dfTidyTaxonomy <- within(dfTidyTaxonomy, kaiju_taxonomy[kaiju_taxonomy == ''] <- 'UN.')
##    dfTidyTaxonomy <- within(dfTidyTaxonomy, kaiju_taxonomy[kaiju_taxonomy == 'NA'] <- 'UN.')
##    dfTidyTaxonomy <- dfTidyTaxonomy %>% mutate(kaiju_taxonomy = ifelse(is.na(kaiju_taxonomy), 'UN.', kaiju_taxonomy))
##    dfTidyTaxonomy$kaiju_taxonomy <- sub(";$", "", dfTidyTaxonomy$kaiju_taxonomy)
##    dfTidyTaxonomy$kaiju_taxonomy <- sub("^VI.$", "UN.", dfTidyTaxonomy$kaiju_taxonomy)
##
##    ##############################
##
##    ###### specific taxa lineage collapsing
##    ### reduce Bacteria to Bacteria
##    dfTidyTaxonomy$kaiju_taxonomy <- gsub("^BA..*", "BA.", dfTidyTaxonomy$kaiju_taxonomy)
##    dfTidyTaxonomy$kaiju_taxonomy <- gsub(".*Arthropoda.*", "Arthropoda", dfTidyTaxonomy$kaiju_taxonomy) ## anything that is arthropod to just arthropod
##    dfTidyTaxonomy$kaiju_taxonomy <- gsub("^EU..*", "EU.", dfTidyTaxonomy$kaiju_taxonomy) ## anything else Eukaryotic, but not Arthropod to just Eukaryotic
##
##    ###dfTidyTaxonomy <- dfTidyTaxonomy %>% filter(grepl('^VI.;', kaiju_taxonomy))
##
##    ### reduce taxa levels (repeat to remove next lowest level)
##    dfTidyTaxonomy$kaiju_taxonomy <- sub(";[^;]*;", ";", dfTidyTaxonomy$kaiju_taxonomy)
##    dfTidyTaxonomy$kaiju_taxonomy <- sub(";[^;]*;", ";", dfTidyTaxonomy$kaiju_taxonomy)
##    dfTidyTaxonomy$kaiju_taxonomy <- sub(";[^;]*;", ";", dfTidyTaxonomy$kaiju_taxonomy)
##    ###dfTidyTaxonomy$kaiju_taxonomy <- sub("^[^;]*;", "", dfTidyTaxonomy$kaiju_taxonomy)
##
##    ####### make higher taxa groupings
##    dfTidyTaxonomy$kaiju_taxonomy_higher <- dfTidyTaxonomy$kaiju_taxonomy
##    dfTidyTaxonomy$kaiju_taxonomy_higher <- gsub("^BA.$", "non_viral", dfTidyTaxonomy$kaiju_taxonomy_higher)
##    dfTidyTaxonomy$kaiju_taxonomy_higher <- gsub("^EU.$", "non_viral", dfTidyTaxonomy$kaiju_taxonomy_higher)
##    dfTidyTaxonomy$kaiju_taxonomy_higher <- gsub("^Arthropoda$", "non_viral", dfTidyTaxonomy$kaiju_taxonomy_higher)
##    dfTidyTaxonomy$kaiju_taxonomy_higher <- gsub("^UN.$", "non_viral", dfTidyTaxonomy$kaiju_taxonomy_higher)
##    dfTidyTaxonomy$kaiju_taxonomy_higher <- gsub("^VI.;", "", dfTidyTaxonomy$kaiju_taxonomy_higher)
##    dfTidyTaxonomy$kaiju_taxonomy_higher <- gsub(";.*", "", dfTidyTaxonomy$kaiju_taxonomy_higher)
##
###    dfTidyTaxonomy$kaiju_taxonomy <- gsub("^VI.;", "", dfTidyTaxonomy$kaiju_taxonomy)
###    dfTidyTaxonomy$kaiju_taxonomy <- sub("^[^;]*;", "", dfTidyTaxonomy$kaiju_taxonomy)
###}
##
###AllContigsKaijuNrEuk <- AllContigsKaijuNrEuk %>% Tidy_Taxonomy() %>% data.frame()
###AllContigsKaijuRVDB <- Tidy_Taxonomy(AllContigsKaijuRVDB)



### Bring together All info for CONTIGS
SeqsRemovedByClustering$SeqSample <- gsub(".contigs.fasta.*", "", SeqsRemovedByClustering$contig)
dfContigs <- SeqsRemovedByClustering %>% 
      left_join(CVQualitySummary, by=c('contig'='contig_id')) %>% 
      left_join(BinList, by=c('contig'='ContigNameDerep')) %>% 
      left_join(BinMmseq, by=c('contig'='Contig')) %>% 
      left_join(AllContigLengths, by=c('contig'='Contig')) %>%
      left_join(RVDBAllContigsDmndBlastx, by=c('contig'='qseqid')) %>% ## .x
      left_join(VOGAllContigsDmndBlastx, by=c('contig'='qseqid')) %>% ## .y
      left_join(AllContigsKaijuNrEuk, by=c('contig'='kaiju_specific_contig')) %>% ## .x
      left_join(AllContigsKaijuRVDB, by=c('contig'='kaiju_specific_contig')) %>% ## .y
      left_join(dfmd, by=c('SeqSample'='specific_read_source'))
dfContigs <- lapply(dfContigs, as.character)
dfContigs <- data.frame(dfContigs)

#### manual assessment of contig info notes:
## NotAssessed: WasContigRemovedByLengthCutoff = YES
## CellularHighConfidence: warnings ~ "no viral genes" && sseqid.x = "NA" && sseqid.y = "NA" && SuperBin.y = "non_viral"
## CellularMedConfidence: warnings ~ "no viral genes" && sseqid.y = "NA" && orig_RepresentativeName.x !~ "Virus" && orig_RepresentativeName.y !~ "Virus"
## CellularLowConfidence: warnings ~ "no viral genes" && orig_RepresentativeName.y !~ "Virus" && checkv_quality = Not-determined
## CellularLowConfidence: checkv_quality = Not-determined && orig_RepresentativeName.y !~ "Virus"

## ViralLowConfidence: checkv_quality = Not-determined && orig_RepresentativeName.y ~ "Virus"
## ViralLowConfidence: warnings ~ "no viral genes" && orig_RepresentativeName.y ~ "Virus"
## ViralLowConfidence: checkv_quality = Low-quality|Medium-quality && orig_RepresentativeName.y !~ "Virus"
## ViralMedConfidence: checkv_quality = Low-quality|Medium-quality && orig_RepresentativeName.y ~ "Virus"
## ViralMedConfidence: evalue.y => 1E-80 =< 1E-20 && orig_RepresentativeName.y ~ "Virus"
## ViralMedConfidence: warnings ~ "no viral genes" && evalue.y =< 1E-80 && orig_RepresentativeName.y ~ "Virus"
## ViralHighConfidence: checkv_quality = "High-quality|Complete"

#### then create column copy of "AnvioBin" as "ManualAssessmentBin", then check each category (Viral/Cellular,Low/Med/High) and remove contaminants from viral bins by deleted the "ManualAssessmentBin" entry from each contaminant contig. Assess the similarity search identities and scores. Look for viral hits that are weak and common cellular homologs such as gag pol or reverse transcriptase, transposases, etc

### Then group sequences from across the bins by "SuperBins" in taxonomic groups according to interest
### Then add another column called "RepresentativeName" and give contigs clean names
###### finally, save the new dfContigs with the three new columns in googlesheets
##### duplicate the sheet with only needed info for reimport (don't keep info that is already available in R)

#### after manually inspecting dfContigs and adding columns, we import it back in
dfContigsManual <- googlesheets4::read_sheet('https://docs.google.com/spreadsheets/d/1b9m_CXN4CNJ63ukONzhT19SiRKKCBnzf5Y1g3SK5Dak/edit#gid=953170762')
dfContigsManual <- lapply(dfContigsManual, as.character)
dfContigsManual <- data.frame(dfContigsManual)

dfContigsManual <- dfContigsManual %>% 
            filter(grepl('NO', WasContigRemovedByLengthCutoff)) %>%
            filter(grepl('NO', WasContigRemovedByDerep)) %>%
            select(contig, SuperBin, RepresentativeName)

####### add the manual assessment info
dfContigs <- dfContigs %>% 
      left_join(dfContigsManual, by='contig')

## add n of contigs per taxa PER sample (using manually curated category RepresentativeName as taxa)
dfContigs <- transform(
      dfContigs, 
      NumContigsPerTaxaPerSample = ave(
      contig, RepresentativeName, Sample_metadata_code, 
      FUN = function(x) sum(!is.na(x))
      ))

dfContigs$NumContigsPerTaxaPerSample <- ifelse(grepl("^\\d+\\.?\\d*$", dfContigs$NumContigsPerTaxaPerSample), as.numeric(dfContigs$NumContigsPerTaxaPerSample), NA)

#### average by genus, species, collection_month, collection_year, apiary, distance (biological replicates)
dfContigs <- transform(
      dfContigs, NumContigsPerTaxaPerSampleAveragedBioReps = ave(
            NumContigsPerTaxaPerSample, 
            RepresentativeName,
            BioRep, 
            FUN = ave
            ))

write.csv(dfContigs, "dfContigs.csv", row.names = FALSE)

###### Info for mapping reads vs contigs
ReadTotalPerSample$specific_read_source <- gsub(".bctrimmedreads.fastq.gz", "", ReadTotalPerSample$specific_read_source)

dfMapping <- ALLSamplesIdxstats %>%
      left_join(ReadTotalPerSample, by=c('ReadSource'='specific_read_source')) %>%
      left_join(dfContigsManual, by=c('Contig'='contig')) %>%
      left_join(dfmd, by=c('ReadSource'='specific_read_source')) %>%
      filter(!grepl('ExcludedByLengthCutoff', SuperBin))
dfMapping <- lapply(dfMapping, as.character)
dfMapping <- data.frame(dfMapping)
dfMapping$NumReadsMapped <- as.numeric(dfMapping$NumReadsMapped)


## n of read mapped per taxa PER sample
dfMapping <- transform(
      dfMapping, NumReadsMappedPerTaxaPerSample = as.numeric(ave(
            NumReadsMapped, 
            RepresentativeName, 
            Sample_metadata_code, 
            FUN = sum
            )))

#### average by genus, species, collection_month, collection_year, apiary, distance (biological replicates)
dfMapping <- transform(
      dfMapping, NumReadsMappedPerTaxaPerSampleAveragedBioReps = as.numeric(ave(
            NumReadsMappedPerTaxaPerSample, 
            RepresentativeName,
            BioRep, 
            FUN = ave
            )))
dfMapping$num_reads_in_specific_read_source <- as.numeric(dfMapping$num_reads_in_specific_read_source)
dfMapping <- transform(
      dfMapping, 
      NumReadsInSampleLibrary = as.numeric(ave(
      num_reads_in_specific_read_source, Sample_metadata_code, 
      FUN = sum)
      ))

########## other metadata filters

## remove Bombus species that are not impatiens
HostFilter <- function(data){
      data %>%
      filter(genus == "Apis" | (genus == "Bombus" & species == "impatiens"))
}

dfContigs <- dfContigs %>% 
      HostFilter()

dfMapping <- dfMapping %>% 
      HostFilter()
#############################################
#############################################
#############################################

####### To plot average by biological replicas, make the following changes to the code below:
## replace "Sample_metadata_code" with "BioRep"
## replace "NumReadsMappedPerTaxaPerSample" with "NumReadsMappedPerTaxaPerSampleAveragedBioReps"

dfMappingUniqBySampleTaxa <- dfMapping %>% subset(select = c(NumReadsMappedPerTaxaPerSample, Sample_metadata_code, RepresentativeName))
dfMappingUniqBySampleTaxa <- dfMappingUniqBySampleTaxa %>% distinct()

colnames(dfMappingUniqBySampleTaxa) <- c("count", "sample", "OTU")
count_matrix <- dfMappingUniqBySampleTaxa %>% pivot_wider(names_from = sample, values_from = count)

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
obj = filterData(obj, present = 1)
obj_CSS = cumNorm(obj, p=cumNormStat(obj))
taxa_obj_CSS = data.frame(MRcounts(obj_CSS, norm=TRUE, log=FALSE))

taxa_obj_CSS <- t(taxa_obj_CSS)


########## ggplot

############################
### identify which samples had contigs for a taxa
###########################
###########################


RenameTaxa <- function(x) {
  x <- gsub("([0-9])", " \\1", x)  # Inserting a space before every number
  toupper(gsub("\\b([[:alnum:]-]{1})[[:alnum:]-]*\\s*", "\\1", x))
}

### pivot table back to long format
taxa_obj_CSS <- taxa_obj_CSS %>% as.data.frame()
taxa_obj_CSS_HM <- taxa_obj_CSS %>% rownames_to_column("samples") %>% pivot_longer(-c(samples), names_to = "taxa", values_to = "counts")
##### add back all sample info and metadata

taxa_obj_CSS_HM <- taxa_obj_CSS_HM %>% left_join(dfmd, by=c('samples'='Sample_metadata_code'))

df_taxa <- dfMapping %>%
            select(SuperBin, RepresentativeName) %>% 
            distinct()
taxa_obj_CSS_HM <- taxa_obj_CSS_HM %>% left_join(df_taxa, by=c('taxa'='RepresentativeName'))
taxa_obj_CSS_HM$distance <- as.character(taxa_obj_CSS_HM$distance)

### rename columns
names(taxa_obj_CSS_HM)[names(taxa_obj_CSS_HM) == 'samples'] <- 'Sample_metadata_code'
names(taxa_obj_CSS_HM)[names(taxa_obj_CSS_HM) == 'taxa'] <- 'RepresentativeName'
names(taxa_obj_CSS_HM)[names(taxa_obj_CSS_HM) == 'counts'] <- 'CSSCountReadsMapped'

###################################################
#### add back in the raw reads counts mapped

### rename columns
names(dfMappingUniqBySampleTaxa)[names(dfMappingUniqBySampleTaxa) == 'sample'] <- 'Sample_metadata_code'
names(dfMappingUniqBySampleTaxa)[names(dfMappingUniqBySampleTaxa) == 'OTU'] <- 'RepresentativeName'
names(dfMappingUniqBySampleTaxa)[names(dfMappingUniqBySampleTaxa) == 'count'] <- 'RawCountReadsMapped'

taxa_obj_CSS_HM <- taxa_obj_CSS_HM %>% left_join(dfMappingUniqBySampleTaxa, by=c('RepresentativeName', 'Sample_metadata_code'))
taxa_obj_CSS_HM$CSSCountReadsMapped[taxa_obj_CSS_HM$CSSCountReadsMapped == 0] <- NA

###################################################
###################################################

######### CHOOSE A FILTER

## filter retain only plants viruses
TaxaFilter <- function(data){
      data %>%
      filter(!grepl('Cellular', SuperBin)) %>%
      filter(grepl('flexiviridae|Bromoviridae|Secoviridae', SuperBin)) %>%
      mutate(SuperBin=gsub("; $", "", SuperBin)) %>%
      mutate(SuperBin=gsub(".*; ", "", SuperBin)) %>%
      mutate(SuperBin=gsub(";$", "", SuperBin))
}

#### filter retain only bee/insect viruses
TaxaFilter <- function(data){
      data %>%
      filter(!grepl('Cellular', SuperBin)) %>%
      filter(grepl('Chronic bee|Negevirus|Reovirales|Virgaviridae|Sinhaliviridae|Phasmaviridae|Rhabdoviridae|Permutotetraviridae|Partitiviridae|Picornavirales|Dicistroviridae|Iflaviridae', SuperBin)) %>%
      filter(!grepl('Secoviridae', SuperBin)) %>%
      filter(!grepl('Turnip', RepresentativeName)) %>%
      filter(!grepl('clover', RepresentativeName)) %>%
      filter(!grepl('Amaranthus', RepresentativeName)) %>%
      filter(!grepl('sativa', RepresentativeName)) %>%
      filter(!grepl('Phocid orthoreovirus', RepresentativeName)) %>%
      mutate(SuperBin=gsub("; $", "", SuperBin)) %>%
      mutate(SuperBin=gsub(".*; ", "", SuperBin)) %>%
      mutate(SuperBin=gsub(";$", "", SuperBin))
}

taxa_obj_CSS_HM$RepresentativeNameABB <- RenameTaxa(taxa_obj_CSS_HM$RepresentativeName)

taxa_obj_CSS_HM_p <- taxa_obj_CSS_HM %>%
      TaxaFilter() %>%
ggplot(aes(x=Sample_metadata_code, y=RepresentativeNameABB, fill=CSSCountReadsMapped)) + 
  geom_tile() +
  scale_fill_viridis(option="inferno", limits = c(1, 2000), oob = scales::oob_squish) +
  facet_nested(SuperBin ~ genus, space = "free", scales = "free") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y=element_text(face="bold", size = 12)) +
  theme(strip.text.y.right = element_text(angle = 0, hjust = 1, face="bold", size = 14)) +
  theme(strip.text.x.top = element_text(face="bold", size = 14)) +
theme(axis.title.y = element_blank()) +
theme(strip.placement = "outside") +
theme(ggh4x.facet.nestline = element_line(colour = "black")) +
labs(fill="# CSS norm. reads")

ggsave(plot=taxa_obj_CSS_HM_p, paste0("Figure_S5", ".pdf"), dpi=300, width = 36, height = 24, units = "cm")
ggsave(plot=taxa_obj_CSS_HM_p, paste0("Figure_S5", ".png"), dpi=300, width = 36, height = 24, units = "cm")

###############################
###############################


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
library("ggtext")

Paired_pal <- RColorBrewer::brewer.pal(12, "Paired")[1:12]


df_taxaPCA <- df_taxa %>% TaxaFilter()

df_taxaPCA <- as.matrix(df_taxaPCA, mode = "character")
df_taxaPCA <- as.data.frame(df_taxaPCA)

taxa_obj_CSS_PCA <- taxa_obj_CSS[apply(taxa_obj_CSS[,-1], 1, function(x) !all(x==0)),] ## remove rows that total to 0 (except first row)
taxa_obj_CSS_PCA <- taxa_obj_CSS_PCA %>% t() %>% as.data.frame()
taxa_obj_CSS_PCA <- taxa_obj_CSS_PCA %>% rownames_to_column("RepresentativeName")

taxa_obj_CSS_PCA <- left_join(df_taxaPCA, taxa_obj_CSS_PCA, by="RepresentativeName") ## keep only taxa kept in the df_taxa

taxa_obj_CSS_PCA <- taxa_obj_CSS_PCA %>% select(-SuperBin) ## remove cols we don't want in plot
taxa_obj_CSS_PCA <- taxa_obj_CSS_PCA %>% column_to_rownames("RepresentativeName")
taxa_obj_CSS_PCA <- t(taxa_obj_CSS_PCA)
taxa_obj_CSS_PCA <- taxa_obj_CSS_PCA[apply(taxa_obj_CSS_PCA[,-1], 1, function(x) !all(x==0)),] ## remove rows that total to 0 (except first row)

taxa_obj_CSS_PCA <- t(taxa_obj_CSS_PCA)
taxa_obj_CSS_PCA <- data.matrix(taxa_obj_CSS_PCA)

OTU = otu_table(taxa_obj_CSS_PCA, taxa_are_rows=TRUE)


rownames(df_taxaPCA) <- df_taxaPCA$RepresentativeName
df_taxaPCA <- as.matrix(df_taxaPCA, mode = "character")
TAX = tax_table(df_taxaPCA)

sampledata = dfmd %>% select(-real_analyses_name, -specific_read_source, -seqrun) %>% distinct() %>% column_to_rownames("Sample_metadata_code") %>% sample_data()

physeq = phyloseq(OTU, TAX, sampledata, remove_undetected = TRUE)
physeq <- phyloseq_validate(physeq)


physeq@sam_data[["Overall_BBHB_ratio"]] <- as.character(physeq@sam_data[["Overall_BBHB_ratio"]])
### calc for PCoA
PCoA_physeq <- physeq %>%
  tax_transform("identity", rank = "RepresentativeName") %>% # don't transform!
  dist_calc("bray") %>%
  ord_calc("PCoA")

### plot PCoA iris
PCoA_physeq %>% 
      ord_plot_iris(tax_level = "RepresentativeName", ord_plot = "above", anno_colour = "genus", n_taxa = 20)

ggsave(plot=last_plot(), paste0("FigR2", ".pdf"), dpi=300, scale=2, units = "cm")
ggsave(plot=last_plot(), paste0("FigR2", ".png"), dpi=300, scale=2, units = "cm")

### plot various PCoA plots
PlotOrdUnconsrained <- function() {
      PCoA_physeq %>%
  ord_plot(color = "genus",
      alpha = 0.8,
      size = 3,
      ) +
  ggside::theme_ggside_void() +
  theme(aspect.ratio = 1) +
  scale_color_manual(values = c("Apis" = Paired_pal[6], "Bombus" = Paired_pal[2]))
}

PlotOrdUnconsrainedOthers <- function(VAR) {
      PCoA_physeq %>%
  ord_plot(shape = "genus", color = VAR, alpha = 0.8, size = 2) +
  ggside::theme_ggside_void() +
  theme(aspect.ratio = 1) +
  guides(fill = guide_legend(reverse = FALSE)) +
  scale_shape_manual(values = c("Apis" = 17, "Bombus" = 16))
}

PlotOrdUnconsrainedOthers2 <- function(VAR) {
      PCoA_physeq %>%
  ord_plot(color = "genus", shape = VAR, alpha = 0.8, size = 2) +
  ggside::theme_ggside_void() +
  theme(aspect.ratio = 1) +
  guides(fill = guide_legend(reverse = FALSE)) +
  scale_color_manual(values = c("Apis" = Paired_pal[6], "Bombus" = Paired_pal[2]))
}

Ord_main <- PlotOrdUnconsrained()


################ variations of normal unconstrained ordination plot with various factors labelled
Ord_main_yearlab <- PlotOrdUnconsrainedOthers("collection_year") + 
  scale_color_manual(name = "Collection year", values = c(
  "2021" = "#fdbf6f",  # light orange
  "2022" = "#ff7f00",  # dark orange
  "2023" = "#b15928"  # dark brown
))

Ord_main_monthlab <- PlotOrdUnconsrainedOthers("collection_month") + 
      
  scale_color_manual(name = "Collection month", values = c(
  "May" = "#b2df8a",  # light green
  "June" = "#33a02c",  # dark green
  "September" = "#fdbf6f",  # light orange
  "October" = "#ff7f00",  # dark orange
  "November" = "#cab2d6",  # light purple
  "July" = "#ffff99", # light yellow
  "August" = "#b15928"  # dark brown
))

Ord_main_dist <- PlotOrdUnconsrainedOthers("distance") + 
      
scale_color_manual(name = "Distance in m from colony", values = c(
  "Colony" = "#fdbf6f",  # light orange
  "100" = "#ff7f00",  # dark orange
  "500" = "#cab2d6",  # light purple
  "1500" = "#6a3d9a" # dark purple
), na.value = "black")

Ord_main_apiary <- PlotOrdUnconsrainedOthers("apiary") + 
      
scale_color_manual(name = "Apiary", values = c(
  "Crosby" = "#33a02c",  # dark green
  "Golf" = "#ff7f00",  # dark orange
  "Vet" = "#6a3d9a" # dark purple
), na.value = "black")

Ord_main_flowers <- PlotOrdUnconsrainedOthers("flower_genus") + 
      
scale_color_manual(name = "Genus of flower collected from", values = c(
  "Solidago" = "#a6cee3",   # light blue
  "Bombus" = "#1f78b4", # dark blue
  "Agastache" = "#b2df8a",  # light green
  "Lotus" = "#33a02c",  # dark green
  "Monarda" = "#fb9a99",  # light red
  "Apis" = "#e31a1c",  # dark red
  "Trifolium" = "#fdbf6f",  # light orange
  "Calamintha" = "#ff7f00",  # dark orange
  "Cirsium" = "#cab2d6",  # light purple
  "Dalea" = "#6a3d9a", # dark purple
  "Eutrochium" = "#ffff99", # light yellow
  "Chamaecrista" = "#b15928"  # dark brown
), na.value = "black")


Ord_main_flower_shape <- PlotOrdUnconsrainedOthers2("flower_shape")
Ord_main_flower_shape_cat <- PlotOrdUnconsrainedOthers2("flower_shape_cat")
Ord_main_flower_depth_cat <- PlotOrdUnconsrainedOthers2("flower_depth_cat")
  

cowplot::plot_grid(Ord_main, Ord_main_yearlab, Ord_main_monthlab, Ord_main_dist, Ord_main_apiary, labels = c('A','B','C','D','E'))

ggsave(plot=last_plot(), paste0("FigR4", ".pdf"), dpi=300, scale=2, units = "cm")
ggsave(plot=last_plot(), paste0("FigR4", ".png"), dpi=300, scale=2, units = "cm")



df <- physeq@tax_table %>% 
      as.data.frame() %>%
      select(RepresentativeName) %>%
      RenameTaxa()
df2 <- physeq@tax_table %>% 
      as.data.frame() %>%
      select(RepresentativeName)

# Extract the comma-separated values
values <- strsplit(df, ",")[[1]]

# Remove the C(' and ')
values <- gsub("C\\(|'|\\)", "", values)

# Print each value on a new line
cat(values, sep = "\n")




 PCA_main_taxa <- physeq %>%
      tax_transform("clr", rank = "RepresentativeName") %>%
      ord_calc(method="PCA"
            ) %>%
      ord_plot(
            plot_taxa = TRUE,
            color = "genus",
            size = 2, alpha = 0.8,
            tax_lab_style = tax_lab_style(colour = "grey30", type = "text", fontface = "bold", max_angle = 90, size = 5),
            tax_vec_style_all  = vec_tax_all(colour = "grey30"),
            constraint_lab_style = constraint_lab_style(colour = "black", type = "text", fontface = "bold", max_angle = 90, size = 3),
            constraint_vec_style  = vec_constraint(colour = "black"),
            taxon_renamer = RenameTaxa
            ) +
      theme(aspect.ratio = 1) +
      scale_color_manual(values = c("Apis" = Paired_pal[6], "Bombus" = Paired_pal[2]))

PCA_main_taxa
ggsave(plot=last_plot(), paste0("Figure_2", ".pdf"), dpi=300, width=24, height=24, units = "cm")
ggsave(plot=last_plot(), paste0("Figure_2", ".png"), dpi=300, width=24, height=24, units = "cm")

#### Redundancy plot to test what variables contributed to communtiy similarity

RDS_year <- physeq %>%
      ps_mutate(
            y_2021 = as.numeric(collection_year == "2021"),
            y_2022 = as.numeric(collection_year == "2022"),
            y_2023 = as.numeric(collection_year == "2023")
            ) %>%
      tax_transform("clr", rank = "RepresentativeName") %>%
      ord_calc(method="RDA", constraints = c(
            "y_2021",
            "y_2022",
            "y_2023"
            )) %>%
      ord_plot(plot_taxa = 1:20, shape = "genus", size = 3, alpha = 0.8, color = "collection_year", 
            tax_lab_style = tax_lab_style(colour = "grey30", type = "text", fontface = "bold", max_angle = 90, size = 5),
            tax_vec_style_all  = vec_tax_all(colour = "grey30"),
            constraint_lab_style = constraint_lab_style(colour = "black", type = "text", fontface = "bold", max_angle = 90, size = 3), constraint_vec_style  = vec_constraint(colour = "black"),
            taxon_renamer = RenameTaxa) +
      theme(aspect.ratio = 1) +
scale_color_manual(values = c("2021" = Paired_pal[3], "2022" = Paired_pal[4], "2023" = Paired_pal[5]))

UnconstrainedTaxaPCAFiltered <- function(INPUT, FILTER_VAR, VAR, PLOT_TAXA, COLOR_VAR, SHAPE_VAR) {
      INPUT %>%
      ps_filter({{ FILTER_VAR }} == VAR) %>%
      tax_transform("clr", rank = "RepresentativeName") %>%
      ord_calc(method="PCA"
            ) %>%
      ord_plot(
            plot_taxa = PLOT_TAXA,
            color = COLOR_VAR,
            shape = SHAPE_VAR,
            size = 3, alpha = 0.8,
            tax_lab_style = tax_lab_style(colour = "grey30", type = "text", fontface = "bold", max_angle = 45, size = 5),
            tax_vec_style_all  = vec_tax_all(colour = "black"),
            taxon_renamer = RenameTaxa) +
      theme(aspect.ratio = 1)
}

Unconstrained_PCA_2021 <- UnconstrainedTaxaPCAFiltered(physeq, collection_year, "2021", FALSE, "genus", 16)
Unconstrained_PCA_2021 <- Unconstrained_PCA_2021 + scale_color_manual(values = c("Apis" = Paired_pal[6], "Bombus" = Paired_pal[2]))

Unconstrained_PCA_2022 <- UnconstrainedTaxaPCAFiltered(physeq, collection_year, "2022", FALSE, "genus", 16)
Unconstrained_PCA_2022 <- Unconstrained_PCA_2022 + scale_color_manual(values = c("Apis" = Paired_pal[6], "Bombus" = Paired_pal[2]))

Unconstrained_PCA_2023 <- UnconstrainedTaxaPCAFiltered(physeq, collection_year, "2023", FALSE, "genus", 16)
Unconstrained_PCA_2023 <- Unconstrained_PCA_2023 + scale_color_manual(values = c("Apis" = Paired_pal[6], "Bombus" = Paired_pal[2]))

Unconstrained_PCA_May <- UnconstrainedTaxaPCAFiltered(physeq, collection_month, "May", FALSE, "genus", 16)
Unconstrained_PCA_May <- Unconstrained_PCA_May + scale_color_manual(values = c("Apis" = Paired_pal[6], "Bombus" = Paired_pal[2]))

Unconstrained_PCA_June <- UnconstrainedTaxaPCAFiltered(physeq, collection_month, "June", FALSE, "genus", 16)
Unconstrained_PCA_June <- Unconstrained_PCA_June + scale_color_manual(values = c("Apis" = Paired_pal[6], "Bombus" = Paired_pal[2]))

Unconstrained_PCA_July <- UnconstrainedTaxaPCAFiltered(physeq, collection_month, "July", FALSE, "genus", 16)
Unconstrained_PCA_July <- Unconstrained_PCA_July + scale_color_manual(values = c("Apis" = Paired_pal[6], "Bombus" = Paired_pal[2]))

Unconstrained_PCA_August <- UnconstrainedTaxaPCAFiltered(physeq, collection_month, "August", FALSE, "genus", 16)
Unconstrained_PCA_August <- Unconstrained_PCA_August + scale_color_manual(values = c("Apis" = Paired_pal[6], "Bombus" = Paired_pal[2]))

Unconstrained_PCA_September <- UnconstrainedTaxaPCAFiltered(physeq, collection_month, "September", FALSE, "genus", 16)
Unconstrained_PCA_September <- Unconstrained_PCA_September + scale_color_manual(values = c("Apis" = Paired_pal[6], "Bombus" = Paired_pal[2]))

Unconstrained_PCA_October <- UnconstrainedTaxaPCAFiltered(physeq, collection_month, "October", FALSE, "genus", 16)
Unconstrained_PCA_October <- Unconstrained_PCA_October + scale_color_manual(values = c("Apis" = Paired_pal[6], "Bombus" = Paired_pal[2]))

Unconstrained_PCA_November <- UnconstrainedTaxaPCAFiltered(physeq, collection_month, "November", FALSE, "genus", 16)
Unconstrained_PCA_November <- Unconstrained_PCA_November + scale_color_manual(values = c("Apis" = Paired_pal[6], "Bombus" = Paired_pal[2]))

design <- "
  AAA#BBB
  AAA#BBB
  AAA#BBB
  AAA#BBB
  CDE#FGH
"

Ord_main + PCA_main_taxa +
Unconstrained_PCA_2021 + Unconstrained_PCA_2022 + Unconstrained_PCA_2023 +
Unconstrained_PCA_July + Unconstrained_PCA_August + Unconstrained_PCA_September +
plot_annotation(tag_levels = 'A') +
plot_layout(guides = 'collect', design = design)


ggsave(plot=last_plot(), paste0("FigR16", ".pdf"), dpi=300, width=36, height=24, units = "cm")
ggsave(plot=last_plot(), paste0("FigR16", ".png"), dpi=300, width=36, height=24, units = "cm")


#####
# plot for transplantation Bombus
physeq_t_Bombus <- physeq %>%
      ps_filter(genus == "Bombus") %>%
      ps_mutate(
            FilterBypass = gsub(".*", "1", collection_year),
            ColonyTransDate = gsub("^(?!.*date[0-9]+$).*$", "Not_transplanted", sample_info3, perl = TRUE),
            ColonyTransDate = ifelse(is.na(ColonyTransDate), "Not_transplanted", ColonyTransDate),
            ColonyTransExp = gsub("^(?!col34$|colCR$).*$", "Not_transplanted", sample_info1, perl = TRUE),
            ColonyTransExp_year = paste0(ColonyTransExp, collection_year),
            CompGroup = paste0(collection_year, collection_month, apiary),
            flower_or_colony = ifelse(grepl("Sentinel|Colonies", Project), "from colony", "from flower"),
            CompGroupForC = paste0(collection_year, flower_or_colony),
            ColonyTransDateCompGroup = paste0(CompGroupForC, " ", apiary, " ", ColonyTransDate)
            )
Unconstrained_PCA_Bombus_t <- UnconstrainedTaxaPCAFiltered(physeq_t_Bombus, FilterBypass, "1", 1:20, "ColonyTransDate", "ColonyTransExp_year")

colors_viridis <- viridis(7, option = "turbo")

# Create the color palette for dates 0 to 6 and add black for date 7
colors_manual <- c(colors_viridis, "grey")

# Plot with scale_colour_manual()
Unconstrained_PCA_Bombus_t <- Unconstrained_PCA_Bombus_t +
      scale_colour_manual(values = colors_manual) +
      guides(color = guide_legend(override.aes = list(size = 6))) +
      theme(legend.text=element_text(size=14))

ggsave(plot=Unconstrained_PCA_Bombus_t, paste0("FigR19", ".pdf"), dpi=300, width=36, height=24, units = "cm")
ggsave(plot=Unconstrained_PCA_Bombus_t, paste0("FigR19", ".png"), dpi=300, width=36, height=24, units = "cm")

Unconstrained_PCA_Bombus_t__TransCompoplot <- cowplot::plot_grid(Unconstrained_PCA_Bombus_t, TransCompoplot, labels = c('A','B'), ncol = 1)



ggsave(plot=Unconstrained_PCA_Bombus_t__TransCompoplot, paste0("Figure_S6", ".pdf"), dpi=300, width=36, height=28, units = "cm")
ggsave(plot=Unconstrained_PCA_Bombus_t__TransCompoplot, paste0("Figure_S6", ".png"), dpi=300, width=36, height=28, units = "cm")


#####

Unconstrained_PCA_2021 <- UnconstrainedTaxaPCAFiltered(physeq, collection_year, "2021", 1:20, "genus", 16) + 
guides(color = guide_legend(override.aes = list(size = 6))) + theme(legend.text=element_text(size=14))

Unconstrained_PCA_2022 <- UnconstrainedTaxaPCAFiltered(physeq, collection_year, "2022", 1:20, "genus", 16) +
guides(color = guide_legend(override.aes = list(size = 6))) + theme(legend.text=element_text(size=14))

Unconstrained_PCA_2023 <- UnconstrainedTaxaPCAFiltered(physeq, collection_year, "2023", 1:20, "genus", 16) +
guides(color = guide_legend(override.aes = list(size = 6))) + theme(legend.text=element_text(size=14))

Unconstrained_PCA_May <- UnconstrainedTaxaPCAFiltered(physeq, collection_month, "May", 1:20, "genus", 16) +
guides(color = guide_legend(override.aes = list(size = 6))) + theme(legend.text=element_text(size=14))

Unconstrained_PCA_June <- UnconstrainedTaxaPCAFiltered(physeq, collection_month, "June", 1:20, "genus", 16) +
guides(color = guide_legend(override.aes = list(size = 6))) + theme(legend.text=element_text(size=14))

Unconstrained_PCA_July <- UnconstrainedTaxaPCAFiltered(physeq, collection_month, "July", 1:20, "genus", 16) +
guides(color = guide_legend(override.aes = list(size = 6))) + theme(legend.text=element_text(size=14))

Unconstrained_PCA_August <- UnconstrainedTaxaPCAFiltered(physeq, collection_month, "August", 1:20, "genus", 16) +
guides(color = guide_legend(override.aes = list(size = 6))) + theme(legend.text=element_text(size=14))

Unconstrained_PCA_September <- UnconstrainedTaxaPCAFiltered(physeq, collection_month, "September", 1:20, "genus", 16) +
guides(color = guide_legend(override.aes = list(size = 6))) + theme(legend.text=element_text(size=14))

Unconstrained_PCA_October <- UnconstrainedTaxaPCAFiltered(physeq, collection_month, "October", 1:20, "genus", 16) +
guides(color = guide_legend(override.aes = list(size = 6))) + theme(legend.text=element_text(size=14))

Unconstrained_PCA_November <- UnconstrainedTaxaPCAFiltered(physeq, collection_month, "November", 1:20, "genus", 16) +
guides(color = guide_legend(override.aes = list(size = 6))) + theme(legend.text=element_text(size=14))


cowplot::plot_grid(
                  Unconstrained_PCA_2021,
                  Unconstrained_PCA_2022,
                  Unconstrained_PCA_2023,
                  #Unconstrained_PCA_May,
                  #Unconstrained_PCA_June,
                  Unconstrained_PCA_July,
                  Unconstrained_PCA_August,
                  Unconstrained_PCA_September,
                  #Unconstrained_PCA_October,
                  #Unconstrained_PCA_November,
                  labels = c('A','B','C','D','E','F'),
                  ncol = 3)

ggsave(plot=last_plot(), paste0("FigR17", ".pdf"), dpi=300, width=36, height=24, units = "cm")
ggsave(plot=last_plot(), paste0("FigR17", ".png"), dpi=300, width=36, height=24, units = "cm")

#####

RDS_month <- physeq %>%
      ps_mutate(
            May = as.numeric(collection_month == "May"),
            Jun = as.numeric(collection_month == "June"),
            Jul = as.numeric(collection_month == "July"),
            Aug = as.numeric(collection_month == "August"),
            Sep = as.numeric(collection_month == "September"),
            Oct = as.numeric(collection_month == "October"),
            Nov = as.numeric(collection_month == "November")
            ) %>%
      tax_transform("clr", rank = "RepresentativeName") %>%
      ord_calc(method="RDA", constraints = c(
            "May",
            "Jun",
            "Jul",
            "Aug",
            "Sep",
            "Oct",
            "Nov"
            )) %>%
      ord_plot(plot_taxa = 1:5, color = "collection_month", size = 3, alpha = 0.8, shape = "genus", 
            tax_lab_style = tax_lab_style(colour = "grey30", type = "text", fontface = "bold", max_angle = 90, size = 5),
            tax_vec_style_all  = vec_tax_all(colour = "grey30"),
            constraint_lab_style = constraint_lab_style(colour = "black", type = "text", fontface = "bold", max_angle = 90, size = 3), constraint_vec_style  = vec_constraint(colour = "black"),
            taxon_renamer = RenameTaxa) +
      theme(aspect.ratio = 1) +
scale_color_manual(values = c("May" = Paired_pal[3], "June" = Paired_pal[4], "July" = Paired_pal[5], "August" = Paired_pal[7], "September" = Paired_pal[8], "October" = Paired_pal[9], "November" = Paired_pal[10])) +
scale_shape_manual(values = c("Apis" = 17, "Bombus" = 16))

RDS_apiary <- physeq %>%
      ps_mutate(
            s_Crosby = as.numeric(apiary == "Crosby"),
            s_Golf = as.numeric(apiary == "Golf"),
            s_Vet = as.numeric(apiary == "Vet"),
            flower_or_colony = ifelse(grepl("Sentinel|Colonies", Project), "from colony", "from flower")
            ) %>%
      tax_transform("clr", rank = "RepresentativeName") %>%
      ord_calc(method="RDA", constraints = c(
            "s_Crosby",
            "s_Golf",
            "s_Vet"
            )) %>%
      ord_plot(plot_taxa = 1:5, shape = "genus", color = "apiary", size = 3, alpha = 0.8,
            tax_lab_style = tax_lab_style(colour = "grey30", type = "text", fontface = "bold", max_angle = 90, size = 5),
            tax_vec_style_all  = vec_tax_all(colour = "grey30"),
            constraint_lab_style = constraint_lab_style(colour = "black", type = "text", fontface = "bold", max_angle = 90, size = 3), constraint_vec_style  = vec_constraint(colour = "black"),
            taxon_renamer = RenameTaxa) +
      theme(aspect.ratio = 1) +
scale_color_manual(values = c("Crosby" = Paired_pal[3], "Golf" = Paired_pal[4], "Vet" = Paired_pal[5])) +
scale_shape_manual(values = c("Apis" = 17, "Bombus" = 16)) +
guides(color = guide_legend(override.aes = list(size = 6))) +
      theme(legend.text=element_text(size=14))

RDS_distance <- physeq %>%
      ps_mutate(
            d_100 = as.numeric(distance == "100"),
            d_500 = as.numeric(distance == "500"),
            d_1500 = as.numeric(distance == "1500"),
            d_Colony = as.numeric(distance == "Colony")
            ) %>%
      tax_transform("clr", rank = "RepresentativeName") %>%
      ord_calc(method="RDA", constraints = c(
            "d_100",
            "d_500",
            "d_1500",
            "d_Colony"
            )) %>%
      ord_plot(plot_taxa = 1:5, shape = "genus", size = 3, alpha = 0.8, color = "distance",
            tax_lab_style = tax_lab_style(colour = "grey30", type = "text", fontface = "bold", max_angle = 90, size = 5),
            tax_vec_style_all  = vec_tax_all(colour = "grey30"),
            constraint_lab_style = constraint_lab_style(colour = "black", type = "text", fontface = "bold", max_angle = 90, size = 3), constraint_vec_style  = vec_constraint(colour = "black"),
            taxon_renamer = RenameTaxa) +
      theme(aspect.ratio = 1) +
scale_color_manual(values = c("100" = Paired_pal[3], "500" = Paired_pal[4], "1500" = Paired_pal[5], "Colony" = Paired_pal[7])) +
scale_shape_manual(values = c("Apis" = 17, "Bombus" = 16)) +
guides(color = guide_legend(override.aes = list(size = 6))) +
      theme(legend.text=element_text(size=14))



valid_genera <- c("Solidago", "Agastache", "Lotus", "Monarda", "Trifolium", "Calamintha", "Cirsium", "Dalea", "Eutrochium", "Chamaecrista")

RDS_flower <- physeq %>%
      ps_filter(collection_year != "2021") %>%
      ps_filter(!grepl("Sentinel|Colonies", Project)) %>%
      ps_mutate(
            Solidago = ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Solidago")),
            Agastache = ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Agastache")),
            Lotus = ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Lotus")),
            Monarda = ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Monarda")),
            Trifolium = ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Trifolium")),
            Calamintha = ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Calamintha")),
            Cirsium = ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Cirsium")),
            Dalea = ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Dalea")),
            Eutrochium = ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Eutrochium")),
            Chamaecrista = ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Chamaecrista")),
            flower_genus = gsub(paste0("^((?!(?:", paste(valid_genera, collapse="|"), ")).)*$"), "Other", flower_genus, perl=TRUE)
            ) %>%
      tax_transform("clr", rank = "RepresentativeName") %>%
      ord_calc(method="RDA", constraints = c(
            "Solidago",
            "Agastache",
            "Lotus",
            "Monarda",
            "Trifolium",
            "Calamintha",
            "Cirsium",
            "Dalea",
            "Eutrochium",
            "Chamaecrista"
            )) %>%
      ord_plot(plot_taxa = 1:10, color = "flower_genus", shape = "genus", size = 3, alpha = 0.8,
            tax_lab_style = tax_lab_style(colour = "grey30", type = "text", fontface = "bold", max_angle = 90, size = 4),
            tax_vec_style_all  = vec_tax_all(colour = "red"),
            constraint_lab_style = constraint_lab_style(colour = "black", type = "text", fontface = "bold", max_angle = 90, size = 3),
            constraint_vec_style  = vec_constraint(colour = "black"),
            taxon_renamer = RenameTaxa) +
      theme(aspect.ratio = 1) +
      scale_color_brewer(palette = "Paired") +
      scale_shape_manual(values = c("Apis" = 17, "Bombus" = 16)) +
      guides(color = guide_legend(override.aes = list(size = 6))) +
      theme(legend.text=element_text(size=14))


physeq_flower_chars <- physeq %>%
      ps_filter(collection_year != "2021") %>%
      ps_filter(flower_shape != "") %>%
      ps_filter(!grepl("Sentinel|Colonies", Project)) %>%
      ps_mutate(
            composite = ifelse(is.na(flower_shape), 0, as.numeric(flower_shape == "composite")),
            cup = ifelse(is.na(flower_shape), 0, as.numeric(flower_shape == "cup")),
            pea = ifelse(is.na(flower_shape), 0, as.numeric(flower_shape == "pea")),
            tube = ifelse(is.na(flower_shape), 0, as.numeric(flower_shape == "tube")),
            simple = ifelse(is.na(flower_shape_cat), 0, as.numeric(flower_shape_cat == "simple")),
            complex = ifelse(is.na(flower_shape_cat), 0, as.numeric(flower_shape_cat == "complex")),
            shallow = ifelse(is.na(flower_depth_cat), 0, as.numeric(flower_depth_cat == "shallow")),
            medium = ifelse(is.na(flower_depth_cat), 0, as.numeric(flower_depth_cat == "medium")),
            deep = ifelse(is.na(flower_depth_cat), 0, as.numeric(flower_depth_cat == "deep"))
            ) 

RDS_flower_chars_shape <- physeq_flower_chars %>%
      tax_transform("clr", rank = "RepresentativeName") %>%
      ord_calc(method="RDA", constraints = c(
            "composite",
            "cup",
            "pea",
            "tube"
            )) %>%
      ord_plot(plot_taxa = 1:10, color = "genus", shape = "flower_shape", size = 3, alpha = 0.8,
            tax_lab_style = tax_lab_style(colour = "grey30", type = "text", fontface = "bold", max_angle = 90, size = 4),
            tax_vec_style_all  = vec_tax_all(colour = "red"),
            constraint_lab_style = constraint_lab_style(colour = "black", type = "text", fontface = "bold", max_angle = 90, size = 3),
            constraint_vec_style  = vec_constraint(colour = "black"),
            taxon_renamer = RenameTaxa) +
      theme(aspect.ratio = 1) +
      scale_color_manual(values = c("Apis" = Paired_pal[6], "Bombus" = Paired_pal[2])) +
      guides(color = guide_legend(override.aes = list(size = 6))) +
      theme(legend.text=element_text(size=14))


RDS_flower_chars_shape_cat <- physeq_flower_chars %>%
      tax_transform("clr", rank = "RepresentativeName") %>%
      ord_calc(method="RDA", constraints = c(
            "simple",
            "complex"
            )) %>%
      ord_plot(plot_taxa = 1:10, color = "genus", shape = "flower_shape_cat", size = 3, alpha = 0.8,
            tax_lab_style = tax_lab_style(colour = "grey30", type = "text", fontface = "bold", max_angle = 90, size = 4),
            tax_vec_style_all  = vec_tax_all(colour = "red"),
            constraint_lab_style = constraint_lab_style(colour = "black", type = "text", fontface = "bold", max_angle = 90, size = 3),
            constraint_vec_style  = vec_constraint(colour = "black"),
            taxon_renamer = RenameTaxa) +
      theme(aspect.ratio = 1) +
      scale_color_manual(values = c("Apis" = Paired_pal[6], "Bombus" = Paired_pal[2])) +
      guides(color = guide_legend(override.aes = list(size = 6))) +
      theme(legend.text=element_text(size=14))


RDS_flower_chars_depth_cat <- physeq_flower_chars %>%
      tax_transform("clr", rank = "RepresentativeName") %>%
      ord_calc(method="RDA", constraints = c(
            "shallow",
            "medium",
            "deep"
            )) %>%
      ord_plot(plot_taxa = 1:10, color = "genus", shape = "flower_depth_cat", size = 3, alpha = 0.8,
            tax_lab_style = tax_lab_style(colour = "grey30", type = "text", fontface = "bold", max_angle = 90, size = 4),
            tax_vec_style_all  = vec_tax_all(colour = "red"),
            constraint_lab_style = constraint_lab_style(colour = "black", type = "text", fontface = "bold", max_angle = 90, size = 3),
            constraint_vec_style  = vec_constraint(colour = "black"),
            taxon_renamer = RenameTaxa) +
      theme(aspect.ratio = 1) +
      scale_color_manual(values = c("Apis" = Paired_pal[6], "Bombus" = Paired_pal[2])) +
      guides(color = guide_legend(override.aes = list(size = 6))) +
      theme(legend.text=element_text(size=14))


RDS_flower_chars_proportion_of_BBs_to_HBs <- physeq_flower_chars %>%
      tax_transform("clr", rank = "RepresentativeName") %>%
      ord_calc(method="RDA", constraints = c(
            "proportion_of_BBs_to_HBs"
            )) %>%
      ord_plot(plot_taxa = 1:10, color = "proportion_of_BBs_to_HBs", shape = "genus", size = 3, alpha = 0.8,
            tax_lab_style = tax_lab_style(colour = "grey30", type = "text", fontface = "bold", max_angle = 90, size = 4),
            tax_vec_style_all  = vec_tax_all(colour = "red"),
            constraint_lab_style = constraint_lab_style(colour = "black", type = "text", fontface = "bold", max_angle = 90, size = 3),
            constraint_vec_style  = vec_constraint(colour = "black"),
            taxon_renamer = RenameTaxa) +
      theme(aspect.ratio = 1) +
      scale_color_viridis_c(option="plasma") +
      guides(color = guide_legend(override.aes = list(size = 6))) +
      theme(legend.text=element_text(size=14))



cowplot::plot_grid(RDS_year, RDS_month, labels = c('A','B'))

ggsave(plot=last_plot(), paste0("FigR14", ".pdf"), dpi=300, width = 48, height = 18, units = "cm")
ggsave(plot=last_plot(), paste0("FigR14", ".png"), dpi=300, width = 48, height = 18, units = "cm")

RDS_apiary_RDS_distance_RDS_flower <- cowplot::plot_grid(RDS_apiary, RDS_distance, RDS_flower, labels = c('A','B','C'))

ggsave(plot=RDS_apiary_RDS_distance_RDS_flower, paste0("FigR15", ".pdf"), dpi=300, width = 36, height = 24, units = "cm")
ggsave(plot=RDS_apiary_RDS_distance_RDS_flower, paste0("FigR15", ".png"), dpi=300, width = 36, height = 24, units = "cm")

RDS_flower_chars <- cowplot::plot_grid(RDS_flower_chars_shape, RDS_flower_chars_shape_cat, RDS_flower_chars_depth_cat, RDS_flower_chars_proportion_of_BBs_to_HBs, labels = c('A','B','C','D'))

RDS_flower_chars <- RDS_flower_chars +
      labs(caption = "Figure R20. PCA with redundancy analyses constraints of bee collection\n flower shape (A), flower shape cat (B), flower depth cat (C), and proportion of bumblebees to\n honey bees at site (D). Taxa influencing point distribution are labelled, as well as variables.\nGO = Ganda orthophasmavirus, BVV4 = Bombus-associated virus Vir4, AHV1 =\n Allermuir Hill virus 1, AHNV = Andrena haemorrhoa nege-like virus, BTI1 =\nBactrocera tryoni iflavirus 1, HPVO231 = Hymenopteran phasma-related virus\n OKIAV231, MV1 = Mayfield virus 1, BVR1 = Bombus-associated virus Reo1, HPV27 =\n Hubei picorna-like virus 27, VVAPV1 = Vespa velutina associated permutotetra-like\n virus 1, ARV = Agassiz Rock virus, CM = Cripavirus mortiferum, CR = Cripavirus\n ropadi, AA = Aparavirus apisacutum, AR1 = Apis rhabdovirus 1, HPV34 = Hubei\n partiti-like virus 34, VVPV2 = Vespa velutina partiti-like virus 2, LSV2 = Lake\n Sinai virus 2, HPVO233 = Hymenopteran phasma-related virus OKIAV233, US =\n Unclassified sinaivirus, LSV = Lake Sinai virus, LSV1 = Lake Sinai virus 1") +
            theme(plot.caption = element_text(hjust = 0.5, size = 14))

ggsave(plot=RDS_flower_chars, paste0("FigR20", ".pdf"), dpi=300, width = 38, height = 36, units = "cm")
ggsave(plot=RDS_flower_chars, paste0("FigR20", ".png"), dpi=300, width = 38, height = 36, units = "cm")


##############################################################################
### Permanova via adonsi2

dfmd_PAN <- dfmd[dfmd$Sample_metadata_code %in% rownames(t(taxa_obj_CSS_PCA)), ] %>%
      select(-seqrun, -real_analyses_name, -specific_read_source) %>%
      distinct()

## collapse the flowers
dfmd_PAN <- dfmd_PAN %>%
      mutate(
            flower_genus = gsub(paste0("^((?!(?:", paste(valid_genera, collapse="|"), ")).)*$"), "Other", flower_genus, perl=TRUE),
            flower_genus = ifelse(is.na(flower_genus), "Other", flower_genus),
            flower_genus = ifelse(collection_year == 2021, "flowers2021", flower_genus),
            flower_shape = ifelse(is.na(flower_shape), "Other", flower_shape),
            flower_shape_cat = ifelse(is.na(flower_shape_cat), "Other", flower_shape_cat),
            flower_depth_cat = ifelse(is.na(flower_depth_cat), "Other", flower_depth_cat),
            proportion_of_BBs_to_HBs = ifelse(is.na(proportion_of_BBs_to_HBs), "Other", proportion_of_BBs_to_HBs)
            )

permanova <- adonis2(t(taxa_obj_CSS_PCA) ~ genus * collection_year * collection_month * apiary * distance * flower_genus, data = dfmd_PAN, method="bray", na.rm = TRUE)

write.table(permanova, "permanova.txt", quote = FALSE, row.names = TRUE, col.names = TRUE)


permanova_flower_chars <- adonis2(t(taxa_obj_CSS_PCA) ~ genus * flower_shape * flower_shape_cat * flower_depth_cat * proportion_of_BBs_to_HBs, data = dfmd_PAN, method="bray", na.rm = TRUE)

write.table(permanova_flower_chars, "permanova_flower_chars.txt", quote = FALSE, row.names = TRUE, col.names = TRUE)

#### composition plot


#########or plto as faceted single figure

myPal <- tax_palette(
  data = physeq, rank = "SuperBin", n = 20, pal = "brewerPlus"
)

topPal <- tax_palette(physeq, pal = "brewerPlus", rank = "RepresentativeName", n = 20, add = NA)
names(topPal) <- sort(names(topPal))
topPal["Other"] <- "black"

CompoPlot <- function(INPUT, GROUP_BY, FILTER_VAR, TAXLEV, MERGE, PALETTE) {
      INPUT %>%
      ps_filter(genus == FILTER_VAR) %>%
      ps_mutate(
            season = recode(collection_month, May = "Spring", June = "Summer", July = "Summer", August = "Summer", September = "Autumn", October = "Autumn", November = "Autumn"),
            ColonyTransDate = paste(sample_info1, "_", sample_info3),
            ColonyTransDate = gsub("^(?!.*date[0-9]+$).*$", "Not_transplanted", ColonyTransDate, perl = TRUE),
            ColonyTransDate = ifelse(is.na(ColonyTransDate), "Not_transplanted", ColonyTransDate),
            CompGroup = paste0(collection_year, apiary),
            flower_or_colony = ifelse(grepl("Sentinel|Colonies", Project), "from colony", "from flower"),
            CompGroupForC = paste0(collection_year, flower_or_colony),
            ColonyTransDateCompGroup = paste0(CompGroupForC, " ", apiary, " ", ColonyTransDate)
            ) %>%
      phyloseq::merge_samples(MERGE) %>%
  comp_barplot(
    tax_level = TAXLEV, n_taxa = 20,
    bar_width = 0.7,
    merge_other = FALSE, bar_outline_colour = NA,
    tax_transform_for_plot = "compositional",
    group_by = GROUP_BY,
    palette = PALETTE,
    sample_order = "bray",
    #sample_order = "bray", # for trans plot only
    tax_order = names(PALETTE),
    other_name = "Other"
  )
}

CompoPlotRename <- function(INPUT, GROUP_BY, FILTER_VAR, TAXLEV, MERGE, PALETTE) {
      INPUT %>%
      ps_filter(genus == FILTER_VAR) %>%
      ps_mutate(
            season = recode(collection_month, May = "Spring", June = "Summer", July = "Summer", August = "Summer", September = "Autumn", October = "Autumn", November = "Autumn"),
            ColonyTransDate = paste(sample_info1, "_", sample_info3),
            ColonyTransDate = gsub("^(?!.*date[0-9]+$).*$", "Not_transplanted", ColonyTransDate, perl = TRUE),
            ColonyTransDate = ifelse(is.na(ColonyTransDate), "Not_transplanted", ColonyTransDate),
            CompGroup = paste0(collection_year, apiary),
            flower_or_colony = ifelse(grepl("Sentinel|Colonies", Project), "from colony", "from flower"),
            CompGroupForC = paste0(collection_year, flower_or_colony),
            ColonyTransDateCompGroup = paste0(CompGroupForC, " ", apiary, " ", ColonyTransDate)
            ) %>%
      phyloseq::merge_samples(MERGE) %>%
  comp_barplot(
    tax_level = TAXLEV, n_taxa = 20,
    bar_width = 0.7,
    merge_other = FALSE, bar_outline_colour = NA,
    tax_transform_for_plot = "compositional",
    group_by = GROUP_BY,
    palette = PALETTE,
    sample_order = "bray",
    #sample_order = "bray", # for trans plot only
    tax_order = names(PALETTE),
    other_name = "Other",
    taxon_renamer = RenameTaxa
  )
}

#### transplatnt colony samples
CompoTrans <- CompoPlotRename(physeq, "genus", "Bombus", "RepresentativeName", "ColonyTransDateCompGroup", topPal)
patchTrans <- patchwork::wrap_plots(CompoTrans, guides = 'collect', ncol = 1)

TransCompoplot <- patchTrans + coord_flip() +
      guides(fill=guide_legend(ncol =1)) +
      guides(color = guide_legend(override.aes = list(size = 6))) +
      theme(legend.text=element_text(size=14))

ggsave(plot=TransCompoplot, paste0("FigR18", ".pdf"), dpi=300, height = 30, width = 36, units = "cm")
ggsave(plot=TransCompoplot, paste0("FigR18", ".png"), dpi=300, height = 30, width = 36, units = "cm")



###########


CompoApis <- CompoPlot(physeq, "collection_year", "Apis", "SuperBin", "CompGroup", myPal)
CompoBombus <- CompoPlot(physeq, "collection_year", "Bombus", "SuperBin", "CompGroup", myPal)

CompoApisFC <- CompoPlot(physeq, "collection_year", "Apis", "SuperBin", "CompGroupForC", myPal)
CompoBombusFC <- CompoPlot(physeq, "collection_year", "Bombus", "SuperBin", "CompGroupForC", myPal)

CompoPlots <- c(CompoApis, CompoBombus, CompoApisFC, CompoBombusFC)

patch <- patchwork::wrap_plots(CompoPlots, nrow = 4, guides = 'collect', heights = c(4, 4, 1, 1)) &
      plot_annotation(tag_levels = 'A')

patch & coord_flip() & guides(fill=guide_legend(ncol =1))

ggsave(plot=last_plot(), paste0("FigR5", ".pdf"), dpi=300, scale=2, units = "cm")
ggsave(plot=last_plot(), paste0("FigR5", ".png"), dpi=300, scale=2, units = "cm")


CompoApis <- CompoPlot(physeq, "collection_year", "Apis", "RepresentativeName", "CompGroup", topPal)
CompoBombus <- CompoPlot(physeq, "collection_year", "Bombus", "RepresentativeName", "CompGroup", topPal)

CompoApisFC <- CompoPlot(physeq, "collection_year", "Apis", "RepresentativeName", "CompGroupForC", topPal)
CompoBombusFC <- CompoPlot(physeq, "collection_year", "Bombus", "RepresentativeName", "CompGroupForC", topPal)

CompoPlots <- c(CompoApis, CompoBombus, CompoApisFC, CompoBombusFC)


patch_CompoPlots <- patchwork::wrap_plots(CompoPlots, nrow = 4, guides = 'collect', heights = c(4, 4, 1, 1)) &
      plot_annotation(tag_levels = 'A')

patch_CompoPlots & coord_flip() & guides(fill=guide_legend(ncol =1))

ggsave(plot=last_plot(), paste0("FigR6.1", ".pdf"), dpi=300, scale=2, units = "cm")
ggsave(plot=last_plot(), paste0("FigR6.1", ".png"), dpi=300, scale=2, units = "cm")

###
CompoPlots <- c(CompoApis, CompoBombus)



patch_CompoPlots <- patchwork::wrap_plots(CompoPlots, byrow = FALSE, nrow = 3, guides = 'collect', heights = c(1, 1)) &
      ##plot_annotation(tag_levels = list(c('B', 'C', 'D', 'E', 'F', 'G')))
      plot_annotation(tag_levels = 'A')
patch_CompoPlots <- patch_CompoPlots & 
      coord_flip() &
      guides(fill=guide_legend(ncol =1, title="")) &
      theme(axis.text = element_text(size = 14), legend.text = element_text(size = 16))

patch_CompoPlots


ggsave(plot=last_plot(), paste0("Figure_3", ".pdf"), dpi=300, height=18, width = 36, units = "cm")
ggsave(plot=last_plot(), paste0("Figure_3", ".png"), dpi=300, height=18, width = 36, units = "cm")

####### composition heatmap
  htmp <- physeq %>%
  ps_mutate(genus = as.character(genus)) %>%
  tax_transform("log2", add = 1, chain = TRUE) %>%
  comp_heatmap(
    taxa = tax_top(physeq, n = 30), grid_col = NA, name = "Log2p",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    colors = heat_palette(),
    row_names_side = "left", row_dend_side = "right", sample_side = "bottom",
    sample_anno = sampleAnnotation(
      genus = anno_sample_cat(
        var = "genus", col = c(Apis = "red", Bombus = "blue"),
        box_col = NA, legend_title = "Host Genus", size = grid::unit(4, "mm")
      )
    )
  )

DrawCompHeatmap <- function() {
      ComplexHeatmap::draw(
  object = htmp, annotation_legend_list = attr(htmp, "AnnoLegends"),
  merge_legends = TRUE
)
}


######### correlations heatmap

# set up the data with numerical variables and filter to top taxa
ReShapeCorrel <- function() {
      physeq %>%
  ps_mutate(
      Apis = if_else(genus == "Apis", true = 1, false = 0),
      Bombus = if_else(genus == "Bombus", true = 1, false = 0),
    Colony = if_else(grepl("Sentinel|Colonies", Project), true = 1, false = 0),
    d_100 = if_else(distance == "100", true = 1, false = 0),
    d_500 = if_else(distance == "500", true = 1, false = 0),
    d_1500 = if_else(distance == "1500", true = 1, false = 0),
    y_2021 = if_else(collection_year == "2021", true = 1, false = 0),
    y_2022 = if_else(collection_year == "2022", true = 1, false = 0),
    y_2023 = if_else(collection_year == "2023", true = 1, false = 0),
      May = ifelse(grepl("May", collection_month), 1, 0),
      June = ifelse(grepl("June", collection_month), 1, 0),
      July = ifelse(grepl("July", collection_month), 1, 0),
      August = ifelse(grepl("August", collection_month), 1, 0),
      September = ifelse(grepl("September", collection_month), 1, 0),
      October = ifelse(grepl("October", collection_month), 1, 0),
      November = ifelse(grepl("November", collection_month), 1, 0),
      #Summer = ifelse(grepl("June|July|August", collection_month), 1, 0),
      #Autumn = ifelse(grepl("September|October|November", collection_month), 1, 0),
      Crosby = if_else(apiary == "Crosby", true = 1, false = 0),
      Golf = if_else(apiary == "Golf", true = 1, false = 0),
      Vet = if_else(apiary == "Vet", true = 1, false = 0),
Agastache= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Agastache")),
Allium= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Allium")),
Arctium= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Arctium")),
Asclepias= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Asclepias")),
Calamintha= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Calamintha")),
Centaurea= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Centaurea")),
Chamaecrista= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Chamaecrista")),
Cirsium= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Cirsium")),
Dalea= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Dalea")),
Eutrochium= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Eutrochium")),
Helianthus= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Helianthus")),
Heliopsis= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Heliopsis")),
Liatris= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Liatris")),
Linaria= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Linaria")),
Lotus= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Lotus")),
Medicago= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Medicago")),
Melilotus= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Melilotus")),
Monarda= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Monarda")),
Nepeta= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Nepeta")),
Persicaria= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Persicaria")),
Pycnanthemum= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Pycnanthemum")),
Ratibida= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Ratibida")),
Rudbeckia= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Rudbeckia")),
Securigera= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Securigera")),
Silphium= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Silphium")),
Solidago= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Solidago")),
Sonchus= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Sonchus")),
Symphyotrichum= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Symphyotrichum")),
Trifolium= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Trifolium")),
Verbena= ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Verbena"))
  )
}

#  draw the heatmap
DrawCorHeatmap <- function(INPUT_PHYSEQ, VARS) {
      INPUT_PHYSEQ %>%
      cor_heatmap(
      vars = VARS,
      grid_col = NA,
      tax_anno = taxAnnotation(
      Prev. = anno_tax_prev()),
      var_anno = varAnnotation(
            Abundance = anno_var_hist(fun = sum, size = grid::unit(15, "mm"))),
      row_names_gp = grid::gpar(fontsize = 8),
      column_names_gp = grid::gpar(fontsize = 8),
      taxon_renamer = RenameTaxa
)
}


psq_many <- ReShapeCorrel()

psq_many_Apis <- ReShapeCorrel() %>% 
      ps_filter(genus == "Apis")

psq_many_Bombus <- ReShapeCorrel() %>%
      ps_filter(genus == "Bombus")

########### split by year and genus

psq_many_2021 <- ReShapeCorrel() %>% 
      ps_filter(collection_year == "2021")

correlhm_genus2021 <- DrawCorHeatmap(psq_many_2021, c("Apis", "Bombus"))

psq_many_2022 <- ReShapeCorrel() %>% 
      ps_filter(collection_year == "2022")

correlhm_genus2022 <- DrawCorHeatmap(psq_many_2022, c("Apis", "Bombus"))

psq_many_2023 <- ReShapeCorrel() %>% 
      ps_filter(collection_year == "2023")

correlhm_genus2023 <- DrawCorHeatmap(psq_many_2023, c("Apis", "Bombus"))

########### split by month and genus (only months with both genera sampled)

psq_many_month <- ReShapeCorrel() %>% 
      ps_filter(collection_month == "July")

correlhm_genusJuly <- DrawCorHeatmap(psq_many_month, c("Apis", "Bombus"))

psq_many_month <- ReShapeCorrel() %>% 
      ps_filter(collection_month == "August")

correlhm_genusAugust <- DrawCorHeatmap(psq_many_month, c("Apis", "Bombus"))

psq_many_month <- ReShapeCorrel() %>% 
      ps_filter(collection_month == "September")

correlhm_genusSeptember <- DrawCorHeatmap(psq_many_month, c("Apis", "Bombus"))


correlhm_genus <- DrawCorHeatmap(psq_many, c("Apis", 
      "Bombus"
      ))

png(paste0("FigR8", ".png"), width=18, height=16, res=300, unit="cm")
correlhm_genus
dev.off()
pdf(paste0("FigR8", ".pdf"), width=7, height=6)
correlhm_genus
dev.off()

correlhmA_year <- DrawCorHeatmap(psq_many_Apis, c("y_2021", "y_2022", "y_2023"))
correlhmA_month <- DrawCorHeatmap(psq_many_Apis, c("May", "June", "July", "August", "September", "October", "November"))

correlhmB_year <- DrawCorHeatmap(psq_many_Bombus, c("y_2021", "y_2022", "y_2023"))
correlhmB_month <- DrawCorHeatmap(psq_many_Bombus, c("July", "August", "September"))

correlhmA_apiary <- DrawCorHeatmap(psq_many_Apis, c("Crosby", "Golf", "Vet", "Colony"))
correlhmB_apiary <- DrawCorHeatmap(psq_many_Bombus, c("Crosby", "Golf", "Vet", "Colony"))

correlhmA_distance <- DrawCorHeatmap(psq_many_Apis, c("d_100", "d_500", "d_1500"))
correlhmB_distance <- DrawCorHeatmap(psq_many_Bombus, c("d_100", "d_500", "d_1500"))

correlhmA_flowers <- psq_many_Apis %>% 
                        ps_filter(collection_year != "2021") %>%
                        DrawCorHeatmap(., c(
      "Agastache",
      "Arctium",
      "Asclepias",
      "Calamintha",
      "Chamaecrista",
      "Cirsium",
      "Dalea",
      "Eutrochium",
      "Heliopsis",
      "Liatris",
      "Linaria",
      "Lotus",
      "Medicago",
      "Melilotus",
      "Monarda",
      "Nepeta",
      "Pycnanthemum",
      "Ratibida",
      "Securigera",
      "Solidago",
      "Sonchus",
      "Trifolium"
))

correlhmB_flowers <- psq_many_Bombus %>% 
                        ps_filter(collection_year != "2021") %>%
                        DrawCorHeatmap(., c(
"Agastache",
"Allium",
"Asclepias",
"Calamintha",
"Chamaecrista",
"Cirsium",
"Dalea",
"Eutrochium",
#"Helianthus", # 2021 only
"Heliopsis",
"Liatris",
"Lotus",
"Medicago",
"Melilotus",
"Monarda",
"Nepeta",
#"Persicaria", # 2021 only
"Pycnanthemum",
"Rudbeckia",
"Securigera",
"Silphium",
"Solidago",
"Sonchus",
#"Symphyotrichum", # not in plant viruses, only in 2021
"Trifolium"
      ))


# List objects starting with "correlhm"
list_correlhm <- ls() %>% grep("^correlhm", ., value = TRUE)

# Loop through each correlhm object
for (obj_name in list_correlhm) {
  # Perform the specified operations
  assign(obj_name, get(obj_name) %>% ComplexHeatmap::draw() %>% grid::grid.grabExpr(), envir = .GlobalEnv)
}


png(paste0("FigR7", ".png"), width=48, height=32, res=300, unit="cm")
cowplot::plot_grid(correlhm_genus2021, correlhm_genus2022, correlhm_genus2023,
correlhm_genusJuly, correlhm_genusAugust, correlhm_genusSeptember,
labels=c("A", "B", "C", "D", "E", "F"),
ncol = 3, nrow = 2)
dev.off()

pdf(paste0("FigR7", ".pdf"), width=18, height=12)
cowplot::plot_grid(correlhm_genus2021, correlhm_genus2022, correlhm_genus2023,
correlhm_genusJuly, correlhm_genusAugust, correlhm_genusSeptember,
labels=c("A", "B", "C", "D", "E", "F"),
ncol = 3, nrow = 2)
dev.off()



png(paste0("FigR9", ".png"), width=36, height=32, res=300, unit="cm")
cowplot::plot_grid(correlhmA_year, correlhmB_year, correlhmA_month, correlhmB_month, labels=c("A", "B", "C", "D"))
dev.off()
pdf(paste0("FigR9", ".pdf"), width=14, height=12)
cowplot::plot_grid(correlhmA_year, correlhmB_year, correlhmA_month, correlhmB_month, labels=c("A", "B", "C", "D"))
dev.off()

png(paste0("FigR10", ".png"), width=36, height=32, res=300, unit="cm")
cowplot::plot_grid(correlhmA_apiary, correlhmB_apiary, correlhmA_distance, correlhmB_distance, labels=c("A", "B", "C", "D"))
dev.off()
pdf(paste0("FigR10", ".pdf"), width=14, height=12)
cowplot::plot_grid(correlhmA_apiary, correlhmB_apiary, correlhmA_distance, correlhmB_distance, labels=c("A", "B", "C", "D"))
dev.off()

png(paste0("FigR11", ".png"), width=48, height=16, res=300, unit="cm")
cowplot::plot_grid(correlhmA_flowers, correlhmB_flowers, labels=c("A", "B"))
dev.off()
pdf(paste0("FigR11", ".pdf"), width=21, height=6)
cowplot::plot_grid(correlhmA_flowers, correlhmB_flowers, labels=c("A", "B"))
dev.off()




      



###################################################################################
###################################################################################




############# breadth of coverage


count_matrix <- SelfSamtoolsDepthCovMappingByWindow %>% mutate(ContigWindow = paste(Contig, start, end, sep = ";")) %>% select(-start, -end, -Contig)
colnames(count_matrix) <- c("sample", "count", "OTU")
count_matrix <- count_matrix %>% pivot_wider(names_from = sample, values_from = count)
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
obj = filterData(obj, present = 1)
obj_CSS = cumNorm(obj, p=cumNormStat(obj))
taxa_obj_CSS = data.frame(MRcounts(obj_CSS, norm=TRUE, log=FALSE))

CoverageDepths <- taxa_obj_CSS %>% rownames_to_column("samples") %>% pivot_longer(-c(samples), names_to = "taxa", values_to = "counts") %>% filter(counts != 0)


CoverageDepths <- CoverageDepths %>% 
      separate(samples, into = c("Contig","start","end"), sep = ";") %>%
      select(-taxa)


##### total number of bases that exceed coverage thresholds (CSS). Since we are using windows for coverage, it must be multiplied by the window size (200) to get the number of bases
count_exceed <- function(x) { sum(x > 5) * 200 }
CoverageDepths5 <- CoverageDepths %>% aggregate(counts ~ Contig, FUN = count_exceed)
colnames(CoverageDepths5) <- c("Contig", "total_bases_over_5_coverage")

count_exceed <- function(x) { sum(x > 10) * 200 }
CoverageDepths10 <- CoverageDepths %>% aggregate(counts ~ Contig, FUN = count_exceed)
colnames(CoverageDepths10) <- c("Contig", "total_bases_over_10_coverage")

count_exceed <- function(x) { sum(x > 100) * 200 }
CoverageDepths100 <- CoverageDepths %>% aggregate(counts ~ Contig, FUN = count_exceed)
colnames(CoverageDepths100) <- c("Contig", "total_bases_over_100_coverage")

CoverageDepths <- CoverageDepths5 %>%
      left_join(CoverageDepths10, by='Contig') %>%
      left_join(CoverageDepths100, by='Contig') %>%
      left_join(AllContigLengths, by='Contig') %>%
      mutate(ContigSource = gsub(".contigs.fasta.*", "", Contig)) %>%
      mutate(ReadSource = gsub(".contigs.fasta.*", "", Contig))

CoverageDepths <- dfContigs %>% 
      select(contig, WasContigRemovedByDerep) %>% 
      distinct() %>%
      right_join(CoverageDepths, by=c('contig'='Contig'))


## add average coverage depth per base for each contig (total coveage per base for contig divided by contig length), then times 100 to get the percent of the contig length in bases that passed the coverage threshold in each category - e.g. 50 % of contig had at least 5X coverage
CoverageDepths <- transform(CoverageDepths, percent_of_contig_len_with_over_5X_coverage = total_bases_over_5_coverage / Length * 100)
CoverageDepths <- transform(CoverageDepths, percent_of_contig_len_with_over_10X_coverage = total_bases_over_10_coverage / Length * 100)
CoverageDepths <- transform(CoverageDepths, percent_of_contig_len_with_over_100X_coverage = total_bases_over_100_coverage / Length * 100)

## get the contig classifications info
dfContigsUniqByConTaxa <- dfContigs %>% 
      select(contig, RepresentativeName, SuperBin) %>%
      distinct()

CoverageDepths <- CoverageDepths %>% 
      left_join(dfContigsUniqByConTaxa, by='contig') %>%
      left_join(dfmd, by=c('ReadSource'='specific_read_source'))

CoverageDepths <- CoverageDepths %>%
      HostFilter()

### Reshape data by gathering data to make coverage a category

CoverageDepths <- gather(CoverageDepths, key = "cov_category", value = "percent_of_contig_len_with_over_nX_coverage", percent_of_contig_len_with_over_5X_coverage, percent_of_contig_len_with_over_10X_coverage, percent_of_contig_len_with_over_100X_coverage)

CoverageDepths$cov_category <- gsub("percent_of_contig_len_with_over_5X_coverage", "> 5x coverage", CoverageDepths$cov_category)
CoverageDepths$cov_category <- gsub("percent_of_contig_len_with_over_10X_coverage", "> 10x coverage", CoverageDepths$cov_category)
CoverageDepths$cov_category <- gsub("percent_of_contig_len_with_over_100X_coverage", "> 100x coverage", CoverageDepths$cov_category)


# Arrange the dataframe by the 'group' column and 'Name' column
CoverageDepths <- CoverageDepths %>% arrange(SuperBin, RepresentativeName)
# Set the order of the 'Name' column based on the alphabetical order of 'group'
CoverageDepths$RepresentativeName <- factor(CoverageDepths$RepresentativeName, levels = unique(CoverageDepths$RepresentativeName))

########## as categoryised coverage thresholds SCATTER

CoverageDepths %>% 
      TaxaFilter() %>%
      filter(Length > 1000) %>%
      filter(grepl('NO', WasContigRemovedByDerep)) %>% 
      filter(grepl('> 10x coverage', cov_category)) %>% 
ggplot(aes(x=Length, y=RepresentativeName)) +
geom_jitter(aes(color = genus, size=percent_of_contig_len_with_over_nX_coverage), alpha = 0.5, width = 0.1, height = 0) + 
facet_grid(SuperBin ~ genus, scales = "free_y", space = "free_y") +
scale_size_continuous(range = c(1,4)) +
theme(strip.text.y.right = element_text(angle = 0, face="bold")) +
theme(strip.text.x.top = element_text(face="bold")) +
theme(axis.text.x=element_text(angle = 90, face="bold")) +
theme(axis.text.y=element_text(face="bold")) +
theme(axis.title.y = element_blank()) +
theme(panel.spacing.x=unit(0.2, "lines")) +
theme(panel.spacing.y=unit(0.2, "lines")) +
theme(panel.background = element_rect(fill = NA, color = "black")) +
theme(strip.placement = "outside") +
theme(ggh4x.facet.nestline = element_line(colour = "black")) +
scale_x_continuous(limits = c(1000, 15000), breaks = seq(0, 15000, 2500)) +
labs(size = "% of contig length with > 10 x coverage",
      color = "Host genus")

ggsave(plot=last_plot(), paste0("FigR12", ".pdf"), dpi=300, scale=2, units = "cm")
ggsave(plot=last_plot(), paste0("FigR12", ".png"), dpi=300, scale=2, units = "cm")

#### contig length only
CoverageDepths_filt <- CoverageDepths %>% 
      TaxaFilter() %>%
      filter(Length > 1000) %>%
      filter(grepl('NO', WasContigRemovedByDerep)) %>% 
      select(-percent_of_contig_len_with_over_nX_coverage, -cov_category, -total_bases_over_5_coverage, -total_bases_over_10_coverage, -total_bases_over_100_coverage) %>%
      distinct() %>%
      filter(!grepl('Chronic bee|Cripavirus|Aparavirus|partiti-like|permutotetra|rhabdovirus|Agassiz|Mayfield|Loch|hasma-related|interruptus|picorna-like|Lake Sinai virus 3|Lake Sinai virus 6|Lake Sinai virus 1|virus Reo1|Unclassified sinaivirus|mosaic virus|Peanut|potexvirus|seed borne|betaflexiviridae|Apple|Prunus|mottle virus|Nepovirus|Turnip|Comovirus|comovirus', RepresentativeName))

CoverageDepths_filt$Length <- as.numeric(CoverageDepths_filt$Length)

##### test if contig length range is different
unique_viruses <- unique(CoverageDepths_filt$RepresentativeName)
results <- list()

for (virus in unique_viruses) {
  # Perform Wilcoxon rank-sum test
  result <- wilcox.test(Length ~ genus, 
                        data = subset(CoverageDepths_filt, RepresentativeName == virus))
  
  # Store results in a list
  results[[virus]] <- result
}

# Write results to a file
output_file <- "wilcox_test_results.txt"
sink(output_file)

# Print results
for (virus in unique_viruses) {
  cat("Virus:", virus, "\n")
  cat("W =", results[[virus]]$statistic, "\n")
  cat("p-value =", results[[virus]]$p.value, "\n\n")
}

# Close the file connection
sink()
#################################

CoverageDepthsBeeSwarm <- CoverageDepths_filt %>% ggplot(aes(genus, Length, fill=genus, color = genus)) + 
geom_beeswarm(size = 0.5) +
facet_nested_wrap(. ~ RepresentativeName, scales = "free", ncol = 3) +
guides(color = guide_legend(override.aes = list(size = 6))) +
      theme(legend.text=element_text(size=14),
      axis.text.y=element_text(size = 12),
      strip.text.x.top = element_text(face="bold", size = 14)
      )

ggsave(plot=CoverageDepthsBeeSwarm, paste0("Figure_S4", ".pdf"), dpi=300, width = 36, height = 48, units = "cm")
ggsave(plot=CoverageDepthsBeeSwarm, paste0("Figure_S4", ".png"), dpi=300, width = 36, height = 48, units = "cm")

########################################################
########################################################

###################################################################################
###################################################################################

## supplementary info

##pdf("Supplementary_Figures_S1-S6.pdf", width = 16, height = 16)


ggsave(plot=RDS_apiary_RDS_distance_RDS_flower, paste0("Figure_S1", ".png"), dpi=300, width = 36, height = 36, units = "cm")
ggsave(plot=RDS_apiary_RDS_distance_RDS_flower, paste0("Figure_S1", ".pdf"), dpi=300, width = 36, height = 36, units = "cm")

#RDS_apiary_RDS_distance_RDS_flower +
      #labs(caption = "Figure S1. PCA with redundancy analyses constraints of bee collection\n apiary (A), distance of collection from apiary (B), and flower genus collected\n from (C). Taxa influencing point distribution are labelled, as well as variables.\nGO = Ganda orthophasmavirus, BVV4 = Bombus-associated virus Vir4, AHV1 =\n Allermuir Hill virus 1, AHNV = Andrena haemorrhoa nege-like virus, BTI1 =\nBactrocera tryoni iflavirus 1, HPVO231 = Hymenopteran phasma-related virus\n OKIAV231, MV1 = Mayfield virus 1, BVR1 = Bombus-associated virus Reo1, HPV27 =\n Hubei picorna-like virus 27, VVAPV1 = Vespa velutina associated permutotetra-like\n virus 1, ARV = Agassiz Rock virus, CM = Cripavirus mortiferum, CR = Cripavirus\n ropadi, AA = Aparavirus apisacutum, AR1 = Apis rhabdovirus 1, HPV34 = Hubei\n partiti-like virus 34, VVPV2 = Vespa velutina partiti-like virus 2, LSV2 = Lake\n Sinai virus 2, HPVO233 = Hymenopteran phasma-related virus OKIAV233, US =\n Unclassified sinaivirus, LSV = Lake Sinai virus, LSV1 = Lake Sinai virus 1") +
      #theme(plot.caption = element_text(hjust = 0.5, size = 16))

#################################

png("Figure_S2.png", width = 24, height = 36, res = 300, units = "cm")

cowplot::plot_grid(correlhm_genus2021, correlhm_genus2022, correlhm_genus2023,
correlhm_genusJuly, correlhm_genusAugust, correlhm_genusSeptember,
labels=c("A", "B", "C", "D", "E", "F"),
ncol = 2, nrow = 3)

dev.off()

pdf("Figure_S2.pdf", width = 9.4, height = 14.2)

cowplot::plot_grid(correlhm_genus2021, correlhm_genus2022, correlhm_genus2023,
correlhm_genusJuly, correlhm_genusAugust, correlhm_genusSeptember,
labels=c("A", "B", "C", "D", "E", "F"),
ncol = 2, nrow = 3)

dev.off()

##################################

png("Figure_S3.png", width = 30, height = 36, res = 300, units = "cm")

      cowplot::plot_grid(
                  Unconstrained_PCA_2021,
                  Unconstrained_PCA_2022,
                  Unconstrained_PCA_2023,
                  Unconstrained_PCA_July,
                  Unconstrained_PCA_August,
                  Unconstrained_PCA_September,
                  labels = c('A','B','C','D','E','F'),
                  ncol = 2)
dev.off()

pdf("Figure_S3.pdf", width = 11.8, height = 14.2)

      cowplot::plot_grid(
                  Unconstrained_PCA_2021,
                  Unconstrained_PCA_2022,
                  Unconstrained_PCA_2023,
                  Unconstrained_PCA_July,
                  Unconstrained_PCA_August,
                  Unconstrained_PCA_September,
                  labels = c('A','B','C','D','E','F'),
                  ncol = 2)
dev.off()


########################################################################################
########################################################################################





####################### AGGREGATE BARCODES

df <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1BMMCJ6y7Slz9MeGH-5rkZFFZlx-c3Mof4L2fDmdkPAY/edit#gid=0", sheet="All_data")

df <- df %>% mutate(across(everything(), as.character)) ## ensure nothing is a stupid list

df <- df %>% filter(grepl('Spillover_FINAL', `Notes (Seq.)`)) ## filter only data for final spillover
df <- as.data.frame(df)

# Aggregate values
aggregated_df <- aggregate(`Path to data (Seq.)` ~ `Sample metadata code`, data = df, FUN = function(x) paste(unique(x), collapse = ";"))

aggregated_df$seqrun_count <- str_count(aggregated_df$`Path to data (Seq.)`, ";") + 1 ## count how many times each sampel seqeunced


# Write aggregated data to a CSV file
write.csv(aggregated_df, "aggregated_data.csv", row.names = FALSE)


list <- unique(df$`Path to data (Seq.)`) %>% data.frame()
write.csv(list, "binning_input_list.csv", row.names = FALSE)
