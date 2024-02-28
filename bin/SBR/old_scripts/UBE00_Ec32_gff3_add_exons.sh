#!/bin/bash
## Command line variables: $1 chromosome files(s)
## converts (Ec32 gff3):
## chr_07	ORCAE	mRNA	513	9851	.	-	.	ID=Ec-07_000010.1;Parent=Ec-07_000010;gene_id=Ec-07_000010;Name=Ec-07_000010.1;length=2130
## chr_07	ORCAE	exon	513	690	.	-	.	ID=Ec-07_000010.1.12;Parent=Ec-07_000010.1;gene_id=Ec-07_000010.1;Name=Ec-07_000010.1;Ontology_term=SO:0000202
## to add exons count at end (fields needed for analyses in []):
## [chr_07]	ORCAE	mRNA	[513]	[9851]	.	-	.	[ID=Ec-07_000010.1];Parent=Ec-07_000010;gene_id=Ec-07_000010;Name=Ec-07_000010.1;length=2130;[exons=1]

######extract 1 file per chromosome of all full gene mRNAs ($3 "\t") sorted by Name= ($4 ";")

awk -F '[\t;]' '$3 ~ /mRNA/ {print > "mRNA_"$1}' ${1};
sort -V -t ";" -k 4,4 mRNA_* | awk '{print > "mRNA2_"$1}';


awk -F '[\t;]' '$3 ~ /exon/ {print > $9}' ${1};	## exract all exons per gene to 1 file per gene

grep -c exon Name* | cat -n > exon_counts;
sed -i 's/\:/\;exons\=/' exon_counts;
awk -F '\t' '{print $2 > "exon_counts2"}' exon_counts;

cut -d ";" -f4 mRNA2_* > mRNA_names; ## makes list of mRNA Name=s
grep -f mRNA_names exon_counts2 | sort -V -t ";" -k 1,1 > exon_counts3;   ### extracts only exons numbers with matching mRNAs (by Name=)
paste -d ";" mRNA2_* exon_counts3 > mRNA_exon_counts; ## merge exon counts to each mRNA
cut -d ";" -f 1-5,7 mRNA_exon_counts | awk '{print > $1}';
sed -i 's/\t/\;/g' chr*;

