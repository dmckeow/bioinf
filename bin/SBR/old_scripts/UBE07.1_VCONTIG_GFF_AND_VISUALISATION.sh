#!/bin/bash
##### this script prepares a final gff file as standard and annotations from UBE05,
##### run script UBE07.1 then UBE07.2 as qsub
##### important outputs of UBE07: UBE07_*.gff (for usual .gff stuff), UBE07_*UGENE.gff & UBE07_*UGENE.fa (merge in UGENE to gb for annotation visualisation) and UBE07_E006_key (annotation abbreviation key)
##### important inputs from other scripts:
UBE05="/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/UBE05/UBE05_002"
ORIGINAL_GFFS="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/02__FINAL_ANNOTATIONS/*.gff"

#awk -F'\t' '$5 ~ /^ID=mRNA/ {print $5 >> "UBE07_A001"}' $UBE05; ## extract names of all UBE genes from all genomes
#find /projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/ -name "UBE02_FINAL_ALL*" > UBE07_A002; ## get all UBE gff info for all genomes
#while read line; do grep -wF -f UBE07_A001 $line >> UBE07_A003; done < UBE07_A002; ## make gff (but missing ; field delimiter) for UBE genes for all genomes

##### the following steps label all UBE with a category 'phexvi=' based on their BLAST result summaries from UBE05 (A_A_A, B_B_B, etc). These will be added to gff as ;features but for the UGENE.gff they will replace $3 - this will determine annotation colour in UGENE
#cut -f 5,43-47 $UBE05 | awk -F'\t' 'BEGIN {OFS="\t"} ; {if($0 ~ "A_A_A" && $0 !~ "B_B_B" && $2 !~ "NCVOG") print $1,"phexvi=viral"}' > UBE07_A004; ## viral; UBE likely viral; AAA, no BBB, not NCVOG
#cut -f 5,43-47 $UBE05 | sed 's/NCVOG[0-9]*/&\_/g' | awk -F'\t' 'BEGIN {OFS="\t"} ; {if($0 ~ "A_A_A" && $0 !~ "B_B_B" && $2 ~ "NCVOG" && $2 !~ /NCVOG0022|NCVOG0038|NCVOG0262|NCVOG0249|NCVOG0023/) print $1,"phexvi=ncvog"}' >> UBE07_A004; ## ncvog; UBE with conserved viral gene hit (NCVOG); AAA, no BBB, with NCVOG, but not 1 of 5 core genes
#cut -f 5,43-47 $UBE05 | sed 's/NCVOG[0-9]*/&\_/g' | awk -F'\t' 'BEGIN {OFS="\t"} ; {if($0 ~ "A_A_A" && $0 !~ "B_B_B" && $2 ~ "NCVOG" && $2 ~ /NCVOG0022|NCVOG0038|NCVOG0262|NCVOG0249|NCVOG0023/) print $1,"phexvi=ncvogcore"}' >> UBE07_A004; ## ncvogcore; UBE with conserved viral gene hit (NCVOG) - AAA, no BBB, with NCVOG, is 1 of 5 core genes
#cut -f 5,44-47 $UBE05 | awk -F'\t' 'BEGIN {OFS="\t"} ; {if($0 ~ "A_A_A" && $0 ~ "B_B_B") print $1,"phexvi=hostviral"}' >> UBE07_A004; ## hostviral; UBE host or virus; AAA and BBB
#cut -f 5,44-47 $UBE05 | awk -F'\t' 'BEGIN {OFS="\t"} ; {if($0 !~ "A_A_A" && $0 ~ "B_B_B") print $1,"phexvi=host"}' >> UBE07_A004; ## host; UBE host; BBB, no AAA
#cut -f 5,44-47 $UBE05 | awk -F'\t' 'BEGIN {OFS="\t"} ; {if($0 ~ "C_C_C" && $0 !~ "A_A_A" && $0 !~ "B_B_B") print $1,"phexvi=orfan"}' >> UBE07_A004; ## orfan; UBE with no hits; CCC, no AAA and no BBB

##### the following steps labels UBE genes with NCVOG codes and functions (these will be UGENE names for ncvog/ncvogcore) and merges them with UBE gene names and phexvi category
#cut -f 43 $UBE05 | sed 's/NCVOG/ncvog=/g' | sed 's/ncvog=[0-9]*/&\tncvog_function=/g' | sed -E 's/;|:|\[|\]|\|| |\.|\/|\(|\)|,|-/_/g' | sed 's/___/__/g' | sed 's/__/_/g' | sed 's/^$/ncvog=none\tncvog_function=none/g' > UBE07_B001; ## reformat NCVOG info

#cut -f 36 $UBE05 | sed 's/ /_/g' | sed 's/\]/\t/g' | cut -f 1 | sed 's/_|\[/\t/g' | sed 's/\.[0-9]*_/&\t/1' | sed 's/\[/\t/g' | sed 's/_\t/\t/g' > UBE07_B001.1; ## get AA BLAST accession,product name, organism and reformat
#awk -F'\t' 'BEGIN{OFS=";"} ; {print "UBE_accession="$1,"UBE_product="$2,"UBE_organism="$3}' UBE07_B001.1 | sed 's/=;/=-;/g' | sed 's/=$/=-/g' > UBE07_B001.2; ## add qualifiers to AA BLAST info (sed step inserts - for any genes with blank AA results - almost orfans)

#cut -f 5 $UBE05 | paste - UBE07_B001 UBE07_B001.2 | sort -V -k 1,1 > UBE07_B002; ## get UBE gene names and merge with UBE07 qualifiers (ncvog_function,AA BLAST info, etc), sorting is important in these steps
#sort -V -k 1,1 UBE07_A004 | paste - UBE07_B002 | cut -f 2,4,5,6 | sed 's/\t/;/g' > UBE07_B003; ## extract and reformat the ncvog and phexvi labels
#sort -V -k 9,9 UBE07_A003 | cut -f 1-8 > UBE07_B004; ## get the \t delimited fields of UBE genes
#sort -V -k 9,9 UBE07_A003 | cut -f 9-15 | sed 's/\t/;/g' | paste UBE07_B004 - | paste -d ";" - UBE07_B003 > UBE07_B005; ## get the ; delimited fields of UBE genes & merge with \t delimited fields; this is the final formatted gff file for the UBE genes

##### the following steps get all other features including non-UBE genes from contigs with UBE clusters
#awk -F'\t' '{print $1}' UBE07_A003 | uniq > UBE07_C001 ## list all UBE contigs
#grep -wF -f UBE07_C001 $ORIGINAL_GFFS | sed 's/\/projet.*://g' > UBE07_C002; ## get all features of UBE contigs
#awk -F'\t|;' '{print $9}' UBE07_A003 > UBE07_C003; ## list all UBE genes
grep -vwF -f UBE07_C003 UBE07_C002 | awk -F'\t|;' 'BEGIN {OFS=";"} ; {if($3 ~ "mRNA") print $0,"phexvi=non_UBE","ncvog=NA","ncvog_function=NA","UBE_accession=-","UBE_product=-","UBE_organism=-"}' > UBE07_C004; ## get all mRNA features from UBE contigs which are non-UBE and add blank UBE qualifiers
awk -F'\t|;' '{if($3 !~ "mRNA") print $0}' UBE07_C002 > UBE07_D001; ## get all the non-mRNA features of UBE contigs

##### the following step makes the final standard gff files which include ALL GENES AND FEATURES (including UBE07 annotations for UBE genes) of CONTIGS WITH UBE GENES ONLY
cat UBE07_B005 UBE07_C004 UBE07_D001 | sort -V -k 1,1 -k 4,5 -k 3,3 | awk -F'ontig' '{print > $1"_UBE07.gff"}';
for file in *_c_UBE07.gff; do mv "$file" "${file//_c_UBE07.gff/_UBE07.gff}"; done;  for file in *_C_UBE07.gff; do mv "$file" "${file//_C_UBE07.gff/_UBE07.gff}"; done; ## fix file names

##### the following steps make the same files as previous step, but replace field $3 with phexvi categories (and rename to stay gff standard), so UGENE will use them to assign annotation colours, AND these steps abbreviate all NCVOG_function names and shorten the gene names, so that they will be visible on annotations
##### BUT FIRST we must make abbreviations for NCVOG functions to look nice in UGENE
awk -F'\t|;' 'BEGIN {OFS="\t"} ; {if($17 !~ "=NA" && $17 !~ "=none") print $18}' *_UBE07.gff | sort -V | uniq | sed '1d' | sed 's/ncvog_function=//g' > UBE07_E001; ## get non-redundant NCVOG function name list
sed 's/^../& \t/' UBE07_E001 | awk -F'\t' '{print tolower($1)}' > UBE07_E002 ## make draft abbreviations from 1st 2 letters
paste UBE07_E001 UBE07_E002 | sort -V -k 2,2 > UBE07_E003; ## merge to edit abbreviations in place manually
awk -F'\t' '{print $2}' UBE07_E003 | uniq -d > UBE07_E004; ## make list of all redundant or too broad abbreviations

