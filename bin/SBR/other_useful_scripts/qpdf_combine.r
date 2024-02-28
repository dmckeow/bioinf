#!/usr/bin/env Rscript

##### set working directory
setwd("/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/IVEX_WGA")

##### load ggplot
library(qpdf)

#### concatenate PDF files
qpdf::pdf_combine(input = c(
  "IVEX_nucmerrep_Porterinema-fluviatile_contig7.pdf",
  "IVEX_nucmerrep_Porterinema-fluviatile_contig12.pdf",
  "IVEX_nucmerrep_Porterinema-fluviatile_contig21.pdf",
  "IVEX_nucmerrep_Porterinema-fluviatile_contig22.pdf",
  "IVEX_nucmerrep_Porterinema-fluviatile_contig40.pdf",
  "IVEX_nucmerrep_Porterinema-fluviatile_contig94.pdf"),
  output = "IVEX_nucmerrep_Porterinema-fluviatile_contigs7_12_21_22_40_94.pdf"
)
