library(ggh4x)
library(ggplot2)
library(tidyverse)

df_reads <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/12Z5JHjL_PFjUF3MLHlOUOjBD8viFRN7FUVJhgLQQCyQ/edit?gid=711310775#gid=711310775", range = "NumReadsMappedSelfViral")

df_contigs <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1NhxQtQRc7nGmX0t0FwOZJAziMn3r0FshtpehCa-Lo5U/edit?gid=476945537#gid=476945537", range = "dfContigs_VirusOnly")

df <- df_contigs |>
  group_by(
    SeqSample
  ) |>
  summarise(
    n_dist_vtax = n_distinct(kaiju_taxonomy.y),
    n_dist_vgenes = n_distinct(sseqid.y),
    sum_vgenes = sum(viral_genes)
  ) |>
  right_join(df_reads, by = join_by(SeqSample == specific_read_source))

df |>
  ggplot(aes(
    num_reads_in_library,
    NumReadsMapped_to_ViralContig_SelfMapping,
    color = genus)) +
  geom_point() +
  geom_smooth()