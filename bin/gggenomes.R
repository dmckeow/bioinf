#!/usr/bin/env Rscript

library(tidyverse)
library(gggenomes)
#setwd("C:/Users/Dean Mckeown/Downloads")

setwd("C:/Users/Dean Mckeown/Downloads/emales_test2")

# parse sequence length and some metadata from fasta file
# note: ex() is just a helper to get stable paths to gggenomes example data
###emales_seqs <- read_seqs(ex("emales/emales.fna")) ## does NOT work on windows

emales_seqs <- read_fai("emales.seqkit.fai") %>%
extract(seq_desc, into = c("emale_type", "is_typespecies"), "=(\\S+) \\S+=(\\S+)", remove=F, convert=T) %>%
arrange(emale_type, length)

emales_seqs <- emales_seqs[1:6,] ## filter to 6 genomes for speed

######## get gene info from gff
emale_genes <- read_gff3("emales.gff") %>%
  mutate(gc_cont=as.numeric(gc_cont))             # per gene GC-content

##emale_genes %<>% mutate(seq_id_start_end = paste0(seq_id, ":", start-1, "-", end)) ## minus one because the cog file has prodigal coords that start at 1, not 0


# prefilter hits by minimum length and maximum divergence
emale_tirs_paf <- read_paf("emales-tirs.paf") %>%
  filter(seq_id == seq_id2 & map_length > 99 & de < 0.1)

emale_tirs <- bind_rows(
  select(emale_tirs_paf, seq_id=seq_id2, start=start2, end=end2, de))
  ##select(emale_tirs_paf, seq_id=seq_id1, start=start1, end=end1, de)) ## something is wrong with filtering TIRs this way - no start1

emale_links <- read_paf("emales.paf") ## worked as intended

emale_gc <- read_bed("emales-gc.tsv")

emale_cogs <- read_tsv("emales-cogs.tsv", col_names = c("feat_id", "cluster_id", "cluster_n"))

emale_cogs %<>% mutate(
  cluster_label = paste0(cluster_id, " (", cluster_n, ")"),
  cluster_label = fct_lump_min(cluster_label, 5, other_level = "rare"),
  cluster_label = fct_lump_min(cluster_label, 15, other_level = "medium"),
  cluster_label = fct_relevel(cluster_label, "rare", after=Inf))


###########################################################
### alternate blast imports

######## from diamond blastx
emale_blastx <- read_feats("emales.dmnd.blastx", format="blast")

############## swap coordinates if reverse
swap_indices <- emale_blastx$start > emale_blastx$end
# Use temporary variable for swapping
temp <- emale_blastx$start[swap_indices]
emale_blastx$start[swap_indices] <- emale_blastx$end[swap_indices]
emale_blastx$end[swap_indices] <- temp
emale_blastx$strand <- ifelse(swap_indices, "-", "+") ## add strand if the coords were swapped (-) or not (+)

emale_blastx %<>% 
  dplyr::rename(blast_desc = seq_id2) %>%
  select(seq_id, start, end, blast_desc) %>%
  distinct()


##emale_blast <- emale_genes %>%
  ##select(seq_id_start_end, feat_id) %>%
  ##merge(emale_blast, by.x="seq_id_start_end", by.y="feat_id", all.y=TRUE) %>%
  ##select(-seq_id_start_end) %>%
  ##as_tibble()
  

####################### searches based on the features - i.e. blast of proteins of predicted genes
emale_blast <- read_subfeats("emales_mavirus-blastp.tsv", format="blast")
emale_blast %<>%
  filter(evalue < 1e-3) %>%
  select(feat_id, start, end, feat_id2) %>%
  left_join(read_tsv("mavirus.faa.tsv", col_names = c("feat_id2", "blast_hit", "blast_desc")))

########### from interproscan of protein
emale_blast <- read_subfeats("emales.iprscan.pseudoblastp", format="blast")
emale_blast %<>%
  select(feat_id, start, end, feat_id2) %>%
  rename(blast_desc = feat_id2) %>%
  filter(grepl('Pfam::', blast_desc))


###########################################################
###########################################################


####### repair the names to match gff, if gene name format is contig:start-end
##emale_blast <- emale_genes %>%
  ##select(seq_id_start_end, feat_id) %>%
  ##merge(emale_blast, by.x="seq_id_start_end", by.y="feat_id", all.y=TRUE) %>%
  ##select(-seq_id_start_end) %>%
  ##as_tibble()

gggenomes(
  genes = emale_genes, seqs = emales_seqs, links=emale_links,
  feats = list(emale_tirs, gc=emale_gc, emale_blastx)) %>%
  add_clusters(emale_cogs) %>%
  add_subfeats(emale_blast) %>%
  sync() +
  geom_feat(position="identity", linewidth=6) +
  geom_seq() +
  geom_link(offset = c(0.3, 0.2), color="white", alpha=.3) +
  geom_bin_label() +
  geom_feat(aes(color="terminal inverted repeat"), data=feats(emale_tirs), linewidth=5, position="pile") +
  geom_feat(aes(color=blast_desc), data=feats(emale_blastx), linewidth=5, position="identity") +
  geom_feat(aes(color=blast_desc), data=feats(emale_blast), linewidth=5, position="pile") +
  geom_gene(aes(fill=cluster_label)) +
  geom_wiggle(aes(z=score, linetype="GC-content"), feats(gc), fill="lavenderblush4", position=position_nudge(y=-.2), height = .2) +
  scale_fill_brewer("Genes", palette="Set3", na.value="cornsilk3") +
  scale_color_brewer("Blast hits & Features", palette="Paired", na.value="cornsilk3") +
  scale_linetype("Graphs") +
  ggtitle(expression(paste("Endogenous mavirus-like elements of ", italic("C. burkhardae"))))


  ##geom_feat_note(aes(label=blast_desc), data=feats(emale_blast), nudge_y=.1, vjust=0)

