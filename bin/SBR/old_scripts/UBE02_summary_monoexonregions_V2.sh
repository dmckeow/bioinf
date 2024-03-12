#!/bin/bash

##### CL VARIABLES for CONTIGS:
## $1 species name e.g. Ecto-Species4-M
## $2 genome fasta
## $3 proteome fasta
##### CL VARIABLES for CHROMOSOMES:
## $1 species name
## $2 directory/ of all input .gff files
## $3 input genome fasta file $4 proteome fasta
##### Input checklist (in working directory):
## UBE01_FINAL_[contig]
## [species][contig] (generated by UBE01)
##### FUNCTION: summarises counts of UBE01 output: uni/bi exonic gene clusters, genes, other genes, cluster sizes
#### CONTIG/CHROMOSOME instructions: @1 and @2 find and replace contig/chromosome option as required
## $2 genome.fasta must be in a file you have write permissions in, as bedtools make an index file there ##
#### Important outputs:
## UBE02_FINAL_ALL_[species]
## xx*
## UBE02_CN010_${1} & UBE02_CL006_${1}
## [species]_UBE_genes_NT.fa & [species]_UBE_genes_AA.fa

############### summarise contig/chromosome info ################
grep -o -c START UBE01_FINAL_* > UBE02_CN001; grep -o -c ID UBE01_FINAL_* > UBE02_CN002; ## COUNT1_${1}  COUNT2_${1}
for file in UBE01_FINAL_*; do echo "$file" >> UBE02_CN003; done; ## FINAL_* ; >> LIST_FINAL_${1}
sed -i 's/UBE01_FINAL_//g' UBE02_CN003; ## LIST_FINAL_${1}
for file in ${1}*ontig*; do echo "$file" >> UBE02_CN004; done; #########@1 CONTIG ${1}_contig* OR CHROMOSOME ${2}chr* ############# LIST_CON1_${1}
grep -f UBE02_CN004 UBE02_CN003 > UBE02_CN005; ## LIST_CON1_${1} LIST_FINAL_${1} > LIST_CON2_${1};
 
for file in $(cat UBE02_CN005); do echo "$file" >> UBE02_CN006; grep -o -c ID "$file" >> UBE02_CN007; done; ## cat LIST_CON2_${1} >> COUNT3.1_${1} >> COUNT3.2_${1}
paste UBE02_CN006 UBE02_CN007 > UBE02_CN008; ## COUNT3.1_${1} COUNT3.2_${1} > COUNT3_${1}
paste UBE02_CN001 UBE02_CN002 UBE02_CN008 > UBE02_CN009; ## COUNT1_${1} COUNT2_${1} COUNT3_${1} > COUNT4_${1};
sed -i 's/\:/\t/g' UBE02_CN009; cut -f 1,2,4,6 UBE02_CN009 > UBE02_CN010_${1}; ## COUNT4_${1} COUNT5_${1}
## field headings are $1 contig/chromosome $2 number of UBE gene clusters $3 number of UBE genes $4 total number of all genes ##

############### summarise cluster info ###############

cat UBE01_FINAL_* > UBE02_FINAL_ALL_${1}; csplit -s -z UBE02_FINAL_ALL_${1} /START/ '{*}'; ## ALL_FINAL_${1}

for file in xx*; do echo "$file" >> UBE02_CL001; grep -c ID "$file" >> UBE02_CL002; head -1 "$file" >> UBE02_CL003; tail -1 "$file" >> UBE02_CL004; done; ## COUNTCL1_${1} COUNTCL2_${1} COUNTCL3_${1} COUNTCL4_${1}

paste UBE02_CL001 UBE02_CL002 UBE02_CL003 UBE02_CL004 > UBE02_CL005; sed -i 's/\t/;/g' UBE02_CL005; cut -d';' -f 1,2,3,6,11,23,27 UBE02_CL005 > UBE02_CL006_${1}; #########@2 CONTIG: 1,2,3,6,11,23,27 OR CHROMOSOME: 1,2,3,6,11,22,26 ## paste COUNTCL1_${1} COUNTCL2_${1} COUNTCL3_${1} COUNTCL4_${1} > COUNTCL5_${1}; COUNTCL5_${1} > COUNTCL6_${1}

#### field headings are $1 xx cluster #, $2 # genes, $3 chromosome/contig, $4 start loci, $5 start gene, $6 end loci, $7 end gene ##

######### extract fastas of all UBE genes  ##############
### first tidy up ALL_FINAL for bedtools (tab delimited, all lines same # fields)
sed -i 's/;START\|;END//g' UBE02_FINAL_ALL_${1}; sed -i 's/\;/\t/g' UBE02_FINAL_ALL_${1};

### PREP BED FILE FOR AMINO ACID & BEDTOOLS ### AAn > UBE02_AAn
cut -f 9 UBE02_FINAL_ALL_${1} > UBE02_AA1; sed -i 's/ID=//g' UBE02_AA1;
cut -f 3 UBE02_FINAL_ALL_${1} > UBE02_AA2; sed -i 's/mRNA/1/g' UBE02_AA2;
cut -f 13 UBE02_FINAL_ALL_${1} > UBE02_AA3; sed -i 's/cds_size=//g' UBE02_AA3; sed -i 's/length=//g' UBE02_AA3; ## AA2.5
awk '{t=($1/3); print t}' UBE02_AA3 > UBE02_AA4; paste UBE02_AA1 UBE02_AA2 UBE02_AA4 > UBE02_AA5; ## AA3 > UBE02_AA4; AA4 > UBE02_AA5

## !NT!: CONTIG: -fi ${2} OR CHROMOSOME: -fi ${3} 
bedtools getfasta -fi $2 -bed UBE02_FINAL_ALL_${1} -fo /projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/UBE_fastas/UBE02_${1}_UBE_genes_NT.fa;
## !AA!: CONTIG: -fi ${3} OR CHROMOSOME: -fi ${4}
bedtools getfasta -fi $3 -bed UBE02_AA5 -fo /projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/UBE_fastas/UBE02_${1}_UBE_genes_AA.fa;

