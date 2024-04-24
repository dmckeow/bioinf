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

#library(ggridges)
#library(ggalluvial)
#library(forcats)
#library(ggrepel)
#library(reshape) #
#library(plotrix) #
#library(tidyr) #
#library(viridisLite) #


theme_set(theme_bw())
set.seed(1234)


################################################
############## IMPORT DATASHEETS ###############
################################################
### save.image(file = here("SpilloverFinal.RData"))
### load(here("SpilloverFinal.RData"))

###### once big data part is done, you can skip the following section and jsut import the big data frames
## write.csv(dfMapping, here("dfMapping.csv"), row.names = FALSE)
## write.csv(dfContigs, here("dfContigs.csv"), row.names = FALSE)

## dfContigs <- read.csv(here("dfContigs.csv"), header=TRUE)
## dfMapping <- read.csv(here("dfMapping.csv"), header=TRUE)


###### Info for mapping reads vs contigs
ALLSamplesIdxstats <- read.csv(here("ALLSamples.idxstats"), header=TRUE, sep="\t")
ReadTotalPerSample <- read.csv(here("read_total_per_sample"), header=TRUE, sep="\t")
AllSamtoolsDepthCovMappingByWindow <- read.csv(here("AllSamtoolsDepthCovMappingByWindow.tsv"), header=TRUE, sep="\t")
##ReadsMappedByContigTotalBasesOver5Cov <- read.csv(here("tmp.ALLstats.list.count1"), header=TRUE, sep="\t")
##ReadsMappedByContigTotalBasesOver10Cov <- read.csv(here("tmp.ALLstats.list.count2"), header=TRUE, sep="\t")
##ReadsMappedByContigTotalBasesOver100Cov <- read.csv(here("tmp.ALLstats.list.count3"), header=TRUE, sep="\t")

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
      mutate(SuperBin=gsub("; $", "", SuperBin)) %>%
      mutate(SuperBin=gsub(".*; ", "", SuperBin)) %>%
      mutate(SuperBin=gsub(";$", "", SuperBin))
}


taxa_obj_CSS_HM_p <- taxa_obj_CSS_HM %>%
      TaxaFilter() %>%
ggplot(aes(x=Sample_metadata_code, y=RepresentativeName, fill=CSSCountReadsMapped)) + 
  geom_tile() +
  scale_fill_viridis(option="inferno", limits = c(1, 2000), oob = scales::oob_squish) +
  facet_nested(SuperBin ~ genus, space = "free", scales = "free") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.text.y=element_text(face="bold")) +
  theme(strip.text.y.right = element_text(angle = 0, hjust = 1, face="bold")) +
  theme(strip.text.x.top = element_text(face="bold")) +
theme(axis.title.y = element_blank()) +
theme(strip.placement = "outside") +
theme(ggh4x.facet.nestline = element_line(colour = "black")) +
labs(fill="# CSS norm. reads") +
xlab("389 Apis samples, 137 Bombus samples")


ggsave(plot=taxa_obj_CSS_HM_p, paste0("HeatmapReadsMappedVsClassifiedContigs", ".pdf"), dpi=300, scale=2, units = "cm")
ggsave(plot=taxa_obj_CSS_HM_p, paste0("HeatmapReadsMappedVsClassifiedContigs", ".png"), dpi=300, scale=2, units = "cm")

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



### calc for PCoA
PCoA_physeq <- physeq %>%
  tax_transform("identity", rank = "RepresentativeName") %>% # don't transform!
  dist_calc("bray") %>%
  ord_calc("PCoA")

### plot PCoA iris
PCoA_physeq %>% 
      ord_plot_iris(tax_level = "RepresentativeName", ord_plot = "above", anno_colour = "genus", n_taxa = 20)

ggsave(plot=last_plot(), paste0("OrdinationIrisGenusVirusTaxaReadsMappedVsClassifiedContigs", ".pdf"), dpi=300, scale=2, units = "cm")
ggsave(plot=last_plot(), paste0("OrdinationIrisGenusVirusTaxaReadsMappedVsClassifiedContigs", ".png"), dpi=300, scale=2, units = "cm")

### plot various PCoA plots
PlotOrdUnconsrained <- function() {
      PCoA_physeq %>%
  ord_plot(color = "genus",
      alpha = 0.8,
      size = 3,
      ) +
  ggside::theme_ggside_void() +
  theme(aspect.ratio = 1) +
  guides(fill = guide_legend(reverse = FALSE))
}

PlotOrdUnconsrainedOthers <- function(VAR) {
      PCoA_physeq %>%
  ord_plot(shape = "genus", color = VAR, alpha = 1.0, size = 3) +
  ggside::theme_ggside_void() +
  theme(aspect.ratio = 1) +
  guides(fill = guide_legend(reverse = FALSE)) +
  scale_shape_manual(values = c("Apis" = 17, "Bombus" = 16))
}

Ord_main <- PlotOrdUnconsrained()

ggsave(plot=Ord_main, paste0("Ordination_Year_ReadsMappedVsClassifiedContigs", ".pdf"), dpi=300, scale=2, units = "cm")
ggsave(plot=Ord_main, paste0("Ordination_Year_ReadsMappedVsClassifiedContigs", ".png"), dpi=300, scale=2, units = "cm")


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
  "Colony3237" = "#a6cee3",   # light blue
  "Colony34" = "#b2df8a",  # light green
  "Crosby" = "#33a02c",  # dark green
  "CoonRapids" = "#fb9a99",  # light red
  "CoonRapids1" = "#fdbf6f",  # light orange
  "Golf" = "#ff7f00",  # dark orange
  "Eagan" = "#cab2d6",  # light purple
  "Vet" = "#6a3d9a", # dark purple
  "Woodbury" = "#ffff99" # light yellow
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

cowplot::plot_grid(Ord_main, Ord_main_yearlab, Ord_main_monthlab, Ord_main_dist, Ord_main_apiary, labels = c('A','B','C','D','E'))

ggsave(plot=last_plot(), paste0("Ordination_OtherFactors_ReadsMappedVsClassifiedContigs", ".pdf"), dpi=300, scale=2, units = "cm")
ggsave(plot=last_plot(), paste0("Ordination_OtherFactors_ReadsMappedVsClassifiedContigs", ".png"), dpi=300, scale=2, units = "cm")


 PCA_main_taxa <- physeq %>%
      tax_transform("clr", rank = "RepresentativeName") %>%
      ord_calc(method="PCA"
            ) %>%
      ord_plot(
            plot_taxa = TRUE,
            color = "genus",
            size = 3, alpha = 0.8,
            constraint_lab_style = constraint_lab_style(colour = "black", type = "text", fontface = "bold", max_angle = 90, size = 3),
            constraint_vec_style  = vec_constraint(colour = "black")) +
      theme(aspect.ratio = 1)

#### Redundancy plot to test what variables contributed to communtiy similarity

RDS_time <- physeq %>%
      ps_mutate(
            y_2021 = as.numeric(collection_year == "2021"),
            y_2022 = as.numeric(collection_year == "2022"),
            y_2023 = as.numeric(collection_year == "2023"),
            May = as.numeric(collection_month == "May"),
            Jun = as.numeric(collection_month == "June"),
            Jul = as.numeric(collection_month == "July"),
            Aug = as.numeric(collection_month == "August"),
            Sep = as.numeric(collection_month == "September"),
            Oct = as.numeric(collection_month == "October"),
            Nov = as.numeric(collection_month == "November"),
            flower_or_colony = ifelse(grepl("Colony", distance), "from colony", "from flower")
            ) %>%
      tax_transform("clr", rank = "RepresentativeName") %>%
      ord_calc(method="RDA", constraints = c(
            "y_2021",
            "y_2022",
            "y_2023",
            "May",
            "Jun",
            "Jul",
            "Aug",
            "Sep",
            "Oct",
            "Nov"
            )) %>%
      ord_plot(plot_taxa = TRUE, color = "genus", size = 3, alpha = 0.8, shape = "flower_or_colony", constraint_lab_style = constraint_lab_style(colour = "black", type = "text", fontface = "bold", max_angle = 90, size = 3), constraint_vec_style  = vec_constraint(colour = "black")) +
      theme(aspect.ratio = 1) +
      scale_shape_manual(values=c(15,16))

ggsave(plot=RDS_time, paste0("OrdinationRDS_Time_ReadsMappedVsClassifiedContigs", ".pdf"), dpi=300, scale=2, units = "cm")
ggsave(plot=RDS_time, paste0("OrdinationRDS_Time_ReadsMappedVsClassifiedContigs", ".png"), dpi=300, scale=2, units = "cm")

RDS_site <- physeq %>%
      ps_mutate(
            s_Crosby = as.numeric(apiary == "Crosby"),
            s_Golf = as.numeric(apiary == "Golf"),
            s_Vet = as.numeric(apiary == "Vet"),
            d_100 = as.numeric(distance == "100"),
            d_500 = as.numeric(distance == "500"),
            d_1500 = as.numeric(distance == "1500"),
            flower_or_colony = ifelse(grepl("Colony", distance), "from colony", "from flower")
            ) %>%
      tax_transform("clr", rank = "RepresentativeName") %>%
      ord_calc(method="RDA", constraints = c(
            "s_Crosby",
            "s_Golf",
            "s_Vet",
            "d_100",
            "d_500",
            "d_1500"
            )) %>%
      ord_plot(plot_taxa = TRUE, color = "genus", size = 3, alpha = 0.8, shape = "flower_or_colony", constraint_lab_style = constraint_lab_style(colour = "black", type = "text", fontface = "bold", max_angle = 90, size = 3), constraint_vec_style  = vec_constraint(colour = "black")) +
      theme(aspect.ratio = 1) +
      scale_shape_manual(values=c(15,16))

RDS_flower <- physeq %>%
      ps_filter(distance != "Colony") %>%
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
            Chamaecrista = ifelse(is.na(flower_genus), 0, as.numeric(flower_genus == "Chamaecrista"))
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
      ord_plot(plot_taxa = 1:10, color = "genus", size = 3, alpha = 0.8,
            tax_lab_style = constraint_lab_style(colour = "red", type = "text", fontface = "bold", max_angle = 90, size = 3),
            tax_vec_style_all  = vec_constraint(colour = "red"),
            constraint_lab_style = constraint_lab_style(colour = "black", type = "text", fontface = "bold", max_angle = 90, size = 3),
            constraint_vec_style  = vec_constraint(colour = "black")) +
      theme(aspect.ratio = 1) +
      scale_shape_manual(values=c(15,16))

cowplot::plot_grid(PCA_main_taxa, RDS_time, RDS_site, RDS_flower, labels = c('A','B','C','D'))

ggsave(plot=last_plot(), paste0("Ordinations_Multiplot_ReadsMappedVsClassifiedContigs", ".pdf"), dpi=300, scale=2, units = "cm")
ggsave(plot=last_plot(), paste0("Ordinations_Multiplot_ReadsMappedVsClassifiedContigs", ".png"), dpi=300, scale=2, units = "cm")

##############################################################################
### Permanova via adonsi2

dfmd_PAN <- dfmd[dfmd$Sample_metadata_code %in% rownames(t(taxa_obj_CSS_PCA)), ] %>%
      select(-seqrun, -real_analyses_name, -specific_read_source) %>%
      distinct()
permanova <- adonis2(t(taxa_obj_CSS_PCA) ~ genus * collection_year * collection_month * apiary * distance, data = dfmd_PAN, method="bray")

permanova2 <- adonis2(t(taxa_obj_CSS_PCA) ~ genus * apiary * distance * flower_genus, data = dfmd_PAN, method="bray", na.action = "na.omit")

write.table(permanova, "permanova.txt", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(permanova2, "permanova2.txt", quote = FALSE, row.names = TRUE, col.names = TRUE)

#### composition plot
CompoPlot <- function(filter_var) {
  physeq %>%
  ps_filter(collection_year == filter_var) %>%
  comp_barplot(
    tax_level = "SuperBin", n_taxa = 19, other_name = "Other",
    palette = distinct_palette(n = 19, add = "grey90"),
    merge_other = FALSE,
    bar_outline_colour = NA,
    facet_by = "genus",
    tax_order = "prev",
    order_with_all_taxa = TRUE,
  ) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
 theme(aspect.ratio = 2)

}

#########or plto as faceted single figure

myPal <- tax_palette(
  data = physeq, rank = "SuperBin", n = 12, pal = "brewerPlus"
)

topPal <- tax_palette(physeq, pal = "brewerPlus", rank = "RepresentativeName", n = 12, add = NA)
names(topPal) <- sort(names(topPal))
topPal["Other"] <- "black"

CompoPlot <- function(FILTER_VAR, TAXLEV, MERGE, PALETTE) {
      physeq %>%
      ps_filter(genus == FILTER_VAR) %>%
      ps_mutate(
            season = recode(collection_month, May = "Spring", June = "Summer", July = "Summer", August = "Summer", September = "Autumn", October = "Autumn", November = "Autumn"),
            CompGroup = paste0(collection_year, season, apiary),
            flower_or_colony = ifelse(grepl("Colony", distance), "from colony", "from flower"),
            CompGroupForC = paste0(collection_year, flower_or_colony)
            ) %>%
      phyloseq::merge_samples(MERGE) %>%
  comp_barplot(
    tax_level = TAXLEV, n_taxa = 12,
    bar_width = 0.7,
    merge_other = FALSE, bar_outline_colour = NA,
    tax_transform_for_plot = "compositional",
    group_by = "collection_year",
    palette = PALETTE,
    sample_order = "bray",
    tax_order = names(PALETTE),
    other_name = "Other"
  )
}

CompoApis <- CompoPlot("Apis", "SuperBin", "CompGroup", myPal)
CompoBombus <- CompoPlot("Bombus", "SuperBin", "CompGroup", myPal)

CompoApisFC <- CompoPlot("Apis", "SuperBin", "CompGroupForC", myPal)
CompoBombusFC <- CompoPlot("Bombus", "SuperBin", "CompGroupForC", myPal)

CompoPlots <- c(CompoApis, CompoBombus, CompoApisFC, CompoBombusFC)

patch <- patchwork::wrap_plots(CompoPlots, nrow = 4, guides = 'collect') &
      plot_annotation(tag_levels = 'A')

patch & coord_flip()

ggsave(plot=last_plot(), paste0("BarplotCompositionReads_HigherTax_MappedVsClassifiedContigs", ".pdf"), dpi=300, scale=2, units = "cm")
ggsave(plot=last_plot(), paste0("BarplotCompositionReads_HigherTax_MappedVsClassifiedContigs", ".png"), dpi=300, scale=2, units = "cm")


CompoApis <- CompoPlot("Apis", "RepresentativeName", "CompGroup", topPal)
CompoBombus <- CompoPlot("Bombus", "RepresentativeName", "CompGroup", topPal)

CompoApisFC <- CompoPlot("Apis", "RepresentativeName", "CompGroupForC", topPal)
CompoBombusFC <- CompoPlot("Bombus", "RepresentativeName", "CompGroupForC", topPal)

CompoPlots <- c(CompoApis, CompoBombus, CompoApisFC, CompoBombusFC)


patch <- patchwork::wrap_plots(CompoPlots, nrow = 4, guides = 'collect') &
      plot_annotation(tag_levels = 'A')

patch & coord_flip()

ggsave(plot=last_plot(), paste0("BarplotCompositionReads_LowerTax_MappedVsClassifiedContigs", ".pdf"), dpi=300, scale=2, units = "cm")
ggsave(plot=last_plot(), paste0("BarplotCompositionReads_LowerTax_MappedVsClassifiedContigs", ".png"), dpi=300, scale=2, units = "cm")



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

png(paste0("HeatmapCompositionReadsMappedVsClassifiedContigs", ".png"), width=18, height=12, res=300, unit="cm")
DrawCompHeatmap()
dev.off()
pdf(paste0("HeatmapCompositionReadsMappedVsClassifiedContigs", ".pdf"), width=8, height=8)
DrawCompHeatmap()
dev.off()

######### correlations heatmap

# set up the data with numerical variables and filter to top taxa
ReShapeCorrel <- function() {
      physeq %>%
  ps_mutate(
      Apis = if_else(genus == "Apis", true = 1, false = 0),
      Bombus = if_else(genus == "Bombus", true = 1, false = 0),
      bimaculatus = if_else(species == "bimaculatus", true = 1, false = 0),
      citrinus = if_else(species == "citrinus", true = 1, false = 0),
      griseocollis = if_else(species == "griseocollis", true = 1, false = 0),
      impatiens = if_else(species == "impatiens", true = 1, false = 0),
      rufocinctus = if_else(species == "rufocinctus", true = 1, false = 0),
      vagans = if_else(species == "vagans", true = 1, false = 0),
    Colony = if_else(distance == "Colony", true = 1, false = 0),
    d_100 = if_else(distance == "100", true = 1, false = 0),
    d_500 = if_else(distance == "500", true = 1, false = 0),
    d_1500 = if_else(distance == "1500", true = 1, false = 0),
    y_2021 = if_else(collection_year == "2021", true = 1, false = 0),
    y_2022 = if_else(collection_year == "2022", true = 1, false = 0),
    y_2023 = if_else(collection_year == "2023", true = 1, false = 0),
      Spring = ifelse(grepl("May", collection_month), 1, 0),
      Summer = ifelse(grepl("June|July|August", collection_month), 1, 0),
      Autumn = ifelse(grepl("September|October|November", collection_month), 1, 0),
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

psq_many <- ReShapeCorrel()

psq_many_Apis <- ReShapeCorrel() %>% 
      ps_filter(genus == "Apis")

psq_many_Bombus <- ReShapeCorrel() %>%
      ps_filter(genus == "Bombus")



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
      column_names_gp = grid::gpar(fontsize = 8)
)
}

correlhm_genusspecies <- DrawCorHeatmap(psq_many, c("Apis", "Bombus", "bimaculatus", "citrinus", "griseocollis", "impatiens", "rufocinctus", "vagans"))

correlhmA_time <- DrawCorHeatmap(psq_many_Apis, c("y_2021", "y_2022", "y_2023", "Spring", "Summer", "Autumn"))
correlhmB_time <- DrawCorHeatmap(psq_many_Bombus, c("y_2021", "y_2022", "y_2023", "Summer", "Autumn"))

correlhmA_site <- DrawCorHeatmap(psq_many_Apis, c("Crosby", "Golf", "Vet", "Colony", "d_100", "d_500", "d_1500"))
correlhmB_site <- DrawCorHeatmap(psq_many_Bombus, c("Crosby", "Golf", "Vet", "Colony", "d_100", "d_500", "d_1500"))

correlhmA_flowers <- DrawCorHeatmap(psq_many_Apis, c(
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

correlhmB_flowers <- DrawCorHeatmap(psq_many_Bombus, c(
"Agastache",
"Allium",
"Asclepias",
"Calamintha",
"Chamaecrista",
"Cirsium",
"Dalea",
"Eutrochium",
"Helianthus",
"Heliopsis",
"Liatris",
"Lotus",
"Medicago",
"Melilotus",
"Monarda",
"Nepeta",
"Persicaria",
"Pycnanthemum",
"Rudbeckia",
"Securigera",
"Silphium",
"Solidago",
"Sonchus",
"Symphyotrichum",
"Trifolium"
      ))

correlhmA_time <- correlhmA_time %>% ComplexHeatmap::draw() %>% grid::grid.grabExpr()
correlhmB_time <- correlhmB_time %>% ComplexHeatmap::draw() %>% grid::grid.grabExpr()

correlhmA_site <- correlhmA_site %>% ComplexHeatmap::draw() %>% grid::grid.grabExpr()
correlhmB_site <- correlhmB_site %>% ComplexHeatmap::draw() %>% grid::grid.grabExpr()

correlhmA_flowers <- correlhmA_flowers %>% ComplexHeatmap::draw() %>% grid::grid.grabExpr()
correlhmB_flowers <- correlhmB_flowers %>% ComplexHeatmap::draw() %>% grid::grid.grabExpr()

png(paste0("correlhm_genusspecies", ".png"), width=18, height=16, res=300, unit="cm")
correlhm_genusspecies
dev.off()
pdf(paste0("correlhm_genusspecies", ".pdf"), width=7, height=6)
correlhm_genusspecies
dev.off()

png(paste0("correlhm_time", ".png"), width=36, height=16, res=300, unit="cm")
cowplot::plot_grid(correlhmA_time, correlhmB_time, labels=c("A", "B"))
dev.off()
pdf(paste0("correlhm_time", ".pdf"), width=14, height=6)
cowplot::plot_grid(correlhmA_time, correlhmB_time, labels=c("A", "B"))
dev.off()

png(paste0("correlhm_site", ".png"), width=36, height=16, res=300, unit="cm")
cowplot::plot_grid(correlhmA_site, correlhmB_site, labels=c("A", "B"))
dev.off()
pdf(paste0("correlhm_site", ".pdf"), width=14, height=6)
cowplot::plot_grid(correlhmA_site, correlhmB_site, labels=c("A", "B"))
dev.off()

png(paste0("correlhm_flowers", ".png"), width=48, height=16, res=300, unit="cm")
cowplot::plot_grid(correlhmA_flowers, correlhmB_flowers, labels=c("A", "B"))
dev.off()
pdf(paste0("correlhm_flowers", ".pdf"), width=21, height=6)
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

ggsave(plot=last_plot(), paste0("ScatterCoverageDepthAcrossContigsSelfReadsOnly", ".pdf"), dpi=300, scale=2, units = "cm")
ggsave(plot=last_plot(), paste0("ScatterCoverageDepthAcrossContigsSelfReadsOnly", ".png"), dpi=300, scale=2, units = "cm")

#### contig length only
CoverageDepths_filt <- CoverageDepths %>% 
      TaxaFilter() %>%
      filter(Length > 1000) %>%
      filter(grepl('NO', WasContigRemovedByDerep)) %>% 
      select(-percent_of_contig_len_with_over_nX_coverage, -cov_category, -total_bases_over_5_coverage, -total_bases_over_10_coverage, -total_bases_over_100_coverage) %>%
      distinct() %>%
      filter(!grepl('Chronic bee|Cripavirus|Aparavirus|partiti-like|permutotetra|rhabdovirus|Agassiz|Mayfield|Loch|hasma-related|interruptus|picorna-like|Lake Sinai virus 3|Lake Sinai virus 6|Lake Sinai virus 1|virus Reo1|Unclassified sinaivirus', RepresentativeName))

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

CoverageDepths_filt %>% ggplot(aes(genus, Length, fill=genus, color = genus)) + 
geom_beeswarm(size = 0.5) +
facet_nested_wrap(. ~ SuperBin + RepresentativeName, scales = "free", ncol = 3) 

ggsave(plot=last_plot(), paste0("ContigSizeBeeswarm", ".pdf"), dpi=300, scale=2, units = "cm")
ggsave(plot=last_plot(), paste0("ContigSizeBeeswarm", ".png"), dpi=300, scale=2, units = "cm")

########################################################
########################################################

###################################################################################
###################################################################################


################# CONTIG coverage across length visualisation vs REF genomes

#map <- merge(map_win, map_ref, by.x=c("specific_read_source", "reference"), by.y=c("specific_read_source", "reference"))

### add metadata
#map$specific_read_source <- gsub("^", "X", map$specific_read_source)
#map$specific_read_source <- gsub(".bctrimmedreads.fastq.gz.*", "", map$specific_read_source)
#map <- merge(map, dfmd, by.x="specific_read_source", by.y=0, all.x=TRUE)

#### calc RPKM
#map <- map %>% mutate(RPKM=paste0(avg_coverage_by_window/((ref_length/1e3)*(NumReadsInSpecificReadSource/1e6))))
#map$RPKM <- as.numeric(map$RPKM)
#map$avg_coverage_by_window <- as.numeric(map$avg_coverage_by_window)
#map <- map %>% mutate(RPKM=(round(RPKM, 1)))
#map <- map %>% mutate(avg_coverage_by_window=(round(avg_coverage_by_window, 1)))

#map <- map %>% mutate(avg_coverage_by_window = log10(avg_coverage_by_window + 1))
#map <- map %>% mutate(RPKM_log = log10(RPKM + 1))



#####################################################

#### average by genus, species, collection_month, collection_year, apiary, distance

#map <- transform(map, avg_coverage_by_window_AVG = ave(avg_coverage_by_window, ref_start, reference, genus, species, collection_month, collection_year, apiary, distance, FUN = ave))

#map <- transform(map, RPKM_AVG = ave(RPKM, ref_start, reference, genus, species, collection_month, collection_year, apiary, distance, FUN = ave))

#map <- map %>% group_by(ref_start, reference, genus, species, collection_month, collection_year, apiary, distance) %>% mutate(AVG_grouping = paste(specific_read_source, collapse = ";")) %>% ungroup()
#map <- data.frame(map)

#map <- map %>% mutate(RPKM_AVG=(round(RPKM_AVG, 1)))
#map <- map %>% mutate(avg_coverage_by_window_AVG=(round(avg_coverage_by_window_AVG, 1)))

#map <- map %>% mutate(avg_coverage_by_window_AVG_log = log10(avg_coverage_by_window_AVG + 1))
#map <- map %>% mutate(RPKM_AVG_log = log10(RPKM_AVG + 1))

#map <- map %>% mutate(RPKM_AVG_log=(round(RPKM_AVG_log, 1)))
#map <- map %>% mutate(avg_coverage_by_window_AVG_log=(round(avg_coverage_by_window_AVG_log, 1)))

#map$seqrun <- gsub("__barcode.*", "", map$specific_read_source)


################################################
#### simple heatmap for mapping vs reference (no interest in windows, just sum total per reference genome)
#map_ref <- map_ref %>% select(specific_read_source, NumReadsInSpecificReadSource) %>% unique()

#map_ref_tot <- merge(map_ref_tot, map_ref, by.x="specific_read_source", by.y="specific_read_source")


#map_ref_tot$mapping_reference <- gsub(".*.bctrimmedreads.fastq.gz.", "", map_ref_tot$specific_read_source)
#map_ref_tot$specific_read_source <- gsub(".bctrimmedreads.fastq.gz.*", "", map_ref_tot$specific_read_source)

#map_ref_tot_md <- merge(map_ref_tot, dfmd, by.x="specific_read_source", by.y="specific_read_source", all.x=TRUE)

#map_ref_tot_md <- map_ref_tot_md %>% mutate(reads_mapped_relative_librarysize=paste0(total_mapped_reads/NumReadsInSpecificReadSource))

#map_ref_tot_md$reads_mapped_relative_librarysize <- as.numeric(map_ref_tot_md$reads_mapped_relative_librarysize)

#map_ref_tot_md %>%
#ggplot(aes(x=specific_read_source, y=mapping_reference, fill=reads_mapped_relative_librarysize)) + 
 # geom_tile() +
  #scale_fill_viridis(option="inferno", limits = c(0, 1), oob = scales::oob_squish) +
  #facet_nested(. ~  genus, switch = "both", space = "free", scales = "free") +
  #theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  #theme(strip.text.y.left = element_text(angle = 0, hjust = 1)) +
#theme(strip.background.y = element_blank()) +
#theme(panel.border = element_blank()) +
#theme(axis.title = element_blank()) +
#theme(panel.spacing.x=unit(0.2, "lines")) +
#theme(panel.spacing.y=unit(0.1, "lines")) +
#theme(strip.placement = "outside") +
#theme(ggh4x.facet.nestline = element_line(colour = "black"))


#####################################################
#####################################################

########## as categoryised coverage thresholds
#map %>%
#filter(grepl('NC_002066.1___Sacbrood_virus', reference)) %>%
#distinct(ref_start, reference, genus, species, collection_month, collection_year, apiary, distance, .keep_all = TRUE) %>% 
#ggplot() +
#geom_line(aes(x=ref_start, y=RPKM_AVG_log, group=AVG_grouping, color=genus))
#geom_area(aes(x=ref_start, y=RPKM_log, group=genus, fill=genus), position = 'stack')


###for jsut as a heatmap
#map %>% 
#filter(grepl('Sacbrood|Deformed|Varroa|Black|LSV_1', reference)) %>%
#distinct(ref_start, reference, genus, species, collection_month, collection_year, apiary, distance, .keep_all = TRUE) %>% 
#ggplot(aes(x=ref_start, y=specific_read_source, fill=RPKM_log)) + 
 # geom_tile() +
  #scale_fill_viridis(option="inferno") +
  #facet_nested(genus ~ reference, switch = "both", space = "free", scales = "free") +
  #theme(axis.text.y=element_blank()) +
  #theme(strip.text.y.left = element_text(angle = 0, hjust = 1)) +
#theme(strip.background.y = element_blank()) +
#theme(panel.border = element_blank()) +
#theme(axis.title = element_blank()) +
#theme(panel.spacing.x=unit(0.2, "lines")) +
#theme(panel.spacing.y=unit(0.1, "lines")) +
#theme(strip.placement = "outside") +
#theme(ggh4x.facet.nestline = element_line(colour = "black"))




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
