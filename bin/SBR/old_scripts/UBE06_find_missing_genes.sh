#!/bin/bash

## $1=UBE05_D001*.gff
## $2="contig"
esvlist="/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/V3/UBE06_EsV1/EsV1_gene_list_Ex"
esvgff="/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/V3/UBE06_EsV1/UBE06_B006_Ectocarpus_siliculosus_virus_1_genome.gff"

grep -wF $2 $1 | cut -d ";" -f 47 | sort -Vu | awk '$0 ~ "=E[0-9]*"' | sed 's/vaa_ssalltitles=//g' | grep -vwF -f - $esvlist | sed 's/E/product=/g' | grep -wF -f - $esvgff | cut -d ";" -f 8,10,11,12 | awk -F ";" '$0 ~ "ncvog_prot_name=NA"' | cut -d ";" -f 1 > "$2"_missing_genes;
grep -wF $2 $1 | cut -d ";" -f 47 | sort -Vu | awk '$0 ~ "=E[0-9]*"' | sed 's/vaa_ssalltitles=//g' | grep -vwF -f - $esvlist | sed 's/E/product=/g' | grep -wF -f - $esvgff | cut -d ";" -f 8,10,11,12 | awk -F ";" '$0 !~ "ncvog_prot_name=NA"' | sed -E 's/(product=[0-9]*);.*uncharacterized.*/\1/g' >> "$2"_missing_genes;
sort -t ";" -V -k 1,1 "$2"_missing_genes -o "$2"_missing_genes;


