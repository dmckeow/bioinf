#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 24
#SBATCH --mem 48GB
#
#
#SBATCH --nodes=1

###### MANUAL STEPS AND PREREQUISITES ######
## 1. Must have run all three IVEX000 scripts, IVEX001.5 (if using), and IVEX001, because this script will run on all genomes together

### run from project folder (which contains archives input script tmp finalresult)

##### command line arguments required #####

################### FINAL OUTPUTS ##########################
## IVEX002_final_001_ncvog_count_contigs_heatmap - for ggplot R heatmap
## IVEX002_final_001_ncvog_count_contigs - for manual viewing/editing in excel etc
## IVEX002_final_002_ncvog_count_genomes_heatmap - for ggplot R heatmap
## IVEX002_final_002_ncvog_count_genomes - for manual viewing/editing in excel etc
## IVEX002_final_003_virus_taxa_percent_count_genomes_heatmap - for ggplot R heatmap
## IVEX002_final_004_rbitscore_scatterplot_[N] - for ggplot R scattergraph
## IVEX002_final_005_contig_summary - for manual viewing/editing in excel etc
## *.gff.gz and IVEX002_final_006.fa.gz - for manual viewing in Geneious, etc
## IVEX002_final_006_clusters - for finding viral gene clusters for manual viewing (not really useful)

################### INPUTS ##########################
###### load software ######
module load seqkit/0.14.0;

##### input files for analysis #####
pwd=$(pwd) ## get path to project folder
FGF=$pwd"/finalresult/IVEX001/IVEX001_final_*.gff"
fa=$pwd"/input/IVEX000_genomes/"
VS=$pwd"/finalresult/IVEX001.5/" ## optional - only needed if using virsorter

##### NCVOG defintion files #####
NCVKEY=$pwd"/input/IVEX000_NCVOG" ## made by script IVEX000_NCVOG.sh

############################### SCRIPT ###############################################

############### STEP 001 ############################
###### NCVOG count per contig (absent NCVOGs (except group1-5) and contigs with no NCVOGs not included)
#### for viewing in excel, etc
cd $pwd"/finalresult/IVEX002"

## list all NCVOGs present in genomes (with viral rbitscore > cellular rbitscore and 0.2)
sed -e 's/vrbitscore=//g' -e 's/rbitscore=//g' $FGF | awk -F "\t|;" 'BEGIN{OFS=";"};{if($3 =="mRNA" && $16 !~ "NA" && $28 !~ "NA" && $18 ~ /ncvog_id=[0-9]+/ && $16 > $28 && $16 >= 0.2) print $18,$21,$22,$23}' | sort -Vu > IVEX002__001_000
## find all contigs with NCVOG hits (with viral rbitscore > cellular rbitscore), split into tmp files by contig
for f in $FGF; do sed -e 's/vrbitscore=//g' -e 's/rbitscore=//g' $f | awk -F "\t|;" -v f=$f 'BEGIN{OFS=";"};{if($3 =="mRNA" && $16 !~ "NA" && $28 !~ "NA" && $18 ~ /ncvog_id=[0-9]+/ && $16 > $28 && $16 >= 0.2) print f,$1,$18,$21,$22,$23 > "tmp_"$1}'; done
### list NCVOGs from Groups1-5 and/or EsV-1/FsV-158
awk -F"\t" '{if($0 ~ "ncvog_group=Group_" || $0 ~ "Ectocarpus_siliculosus_virus_1" || $0 ~ "Feldmannia_species_virus") print $5";"$7";"$8";"$9}' $NCVKEY | sort -Vu > IVEX002__001_001
## reduce to only NCVOGs from Groups1-5 and/or EsV-1/FsV-158 and present in your genomes
cat IVEX002__001_000 IVEX002__001_001 | sort -Vu > tmp && mv tmp IVEX002__001_001
### merge genome and contig name, in case names are ambiguous
## fix merge of genome and contig name
for f in tmp_*; do sed -i -E 's/^.+IVEX001_final_(.+)\.gff;(.+)(;ncvog_id=.+$)/\1__\2\3/g' $f; done
## remove name column and rename fiel per contig to include genome name
for f in tmp_*; do awk -F ";" 'BEGIN{OFS=";"};{print $2,$3,$4,$5 > "2tmp_"$1}' $f; done
#### prep counts per NCVOG, each file contig
for f in 2tmp_*; do cat IVEX002__001_001 $f | sort -t";" -k 2,2V -k 1,1V | uniq -c | sed -E 's/ +([0-9]+) /\1\t/g' | awk -F"\t" -v f=$f '{print f"\t"$1-1"\t"$2}' > tmp && mv tmp $f; done
rm -f tmp_*; for f in 2tmp_*; do mv "$f" "${f//2tmp_/tmp_}"; done

#### prep counts for EsV-1 and FsV-158 references
awk -F "\t" 'BEGIN{OFS=";"};{if($1 ~"Ectocarpus_siliculosus_virus_1") print $5,$7,$8,$9}' $NCVKEY | cat - IVEX002__001_001 | sort -t";" -k 2,2V -k 1,1V | uniq -c | sed -Ee 's/ +([0-9]+) /\1\t/g' | awk -F"\t" '{print "EsV-1\t"$1-1"\t"$2}' > tmp_EsV-1
awk -F "\t" 'BEGIN{OFS=";"};{if($1 ~"Feldmannia_species_virus") print $5,$7,$8,$9}' $NCVKEY | cat - IVEX002__001_001 | sort -t";" -k 2,2V -k 1,1V | uniq -c | sed -Ee 's/ +([0-9]+) /\1\t/g' | awk -F"\t" '{print "FsV-158\t"$1-1"\t"$2}' > tmp_FsV-158
## get first column for table (ncvog label)
sort -t";" -k 2,2V -k 1,1V IVEX002__001_001 | sed -e 's/ncvog_id=/_/g' -Ee 's/ncvog_group=Group_|ncvog_legend=|ncvog_abbrev=//g' -e 's/ncvog_group=Not_assigned/na/g' > IVEX002__001_002

##### get header for table (contig label)
cat tmp_* | cut -f 1 | awk '!a[$0]++' | sed -e 's/$/;/g' -Ee 's/^2tmp_|^tmp_//g' | tr -d '\n' | sed 's/;$/\n/g' | sed 's/^/ncvog_id;ncvog_group;ncvog_abbrev;ncvog_legend;/g' > IVEX002__001_003

#### prepare heatmap for R (for manually picking out specific contigs later)
cat tmp_* | sed 's/^2tmp_/tmp_/g' | sort -t ";" -k 2,2Vr | sed -e 's/;/\t/g' -e '/ncvog_group=Not_assigned/d' -Ee 's/ncvog_group=Group_|ncvog_id=|ncvog_legend=|ncvog_abbrev=//g' -e 's/\t0\t/\t\t/g' -e 's/^tmp_//g' -e 's/^EsV-1\t/EsV-1__EsV-1\t/g' -e 's/^tmp_//g' -e 's/^FsV-158\t/FsV-158__FsV-158\t/g' > tmp2_IVEX002_final_001_ncvog_count_contigs_heatmap

#### remove contigs with no NCVOGs and split by genome
awk -F"\t" '{if($2 !="") print $1}' tmp2_IVEX002_final_001_ncvog_count_contigs_heatmap | sort -Vu | grep -wF -f - tmp2_IVEX002_final_001_ncvog_count_contigs_heatmap | awk -F "__" '{print $0 > "IVEX002_final_001_"$1"_ncvog_count_contigs_heatmap"}'

### rename reference genome heatmaps
mv IVEX002_final_001_EsV-1_ncvog_count_contigs_heatmap tmp3_EsV-1
mv IVEX002_final_001_FsV-158_ncvog_count_contigs_heatmap tmp3_FsV-158

## remerge with reference and add header
for f in IVEX002_final_001_*_ncvog_count_contigs_heatmap; do cat tmp3_* $f | sed -z 's/^/contig\tcount\tncvog\tgroup\tabbrev\tlegend\n/1' | sed -E 's/^.+__//1' > tmp && mv tmp $f; done
for f in IVEX002_final_001_*_ncvog_count_contigs_heatmap; do awk -F "\t" -v f=$f '{print $0"\t"f}' $f | sed -Ee 's/\tIVEX002_final_001_(.+)_ncvog_count_contigs_heatmap/\t\1/g' -Ee 's/(\tlegend\t).+/\1genome/g' > tmp && mv tmp $f; done

#### merge into final file for excel viewing
for f in tmp_*; do cut -f 2 $f > tmp; cp tmp $f; done;
find . -maxdepth 1 -type f -name tmp_\* | sort | split -l 1000 -d - tmp_lists;
for f in tmp_lists*; do paste -d ";" $(cat $f) > tmp_merge_$(basename $f); done; paste -d ";" IVEX002__001_002 tmp_merge* > IVEX002_final_001_ncvog_count_contigs

cat IVEX002__001_003 IVEX002_final_001_ncvog_count_contigs > tmp && mv tmp IVEX002_final_001_ncvog_count_contigs

rm -f tmp*
### these steps swap the rows and columns (transpoase) IVEX002_final_001_ncvog_count_contigs and add explanation header
split -d -l 1 IVEX002_final_001_ncvog_count_contigs tmp_
sed -i 's/;/\n/g' tmp_*;
find . -maxdepth 1 -type f -name tmp_\* | sort | split -l 1000 -d - tmp_lists;
for f in tmp_lists*; do paste -d ";" $(cat $f) > tmp_merge_$(basename $f); done; paste -d ";" tmp_merge* > tmp; mv tmp IVEX002_final_001_ncvog_count_contigs;

sed -i -z 's/^/## no. of copies of core Nucleocytoviricota genes (NCVOGs) per contig (with viral rbitscore > cellular rbitscore and 0.2). Group refers to level of gene conservation, 1 being the most conserved core genes.\n/g' IVEX002_final_001_ncvog_count_contigs
rm -f tmp*

############### STEP 002 ############################
##### counts of virus contigs checklist
## get genome, contig rbitscores, and virus hit type key file for all contigs
for f in $FGF; do sed -Ee 's/vrbitscore=|ncvog_group=|vlineage=//g' -e 's/rbitscore=//g' $f | awk -F "\t|;" -v f=$f 'BEGIN{OFS="\t"};{if($3 =="mRNA" && $16 !~ "NA" && $28 !~ "NA" && $10 ~ "vsalltitles=" && $10 !~ "vsalltitles=none" && $16 >= 0.2) print f,$1,$16,$28,$21,$17}'; done | sed -Ee 's/Group_1|Group_2|Group_3/Group_1-3/g' -Ee 's/Group_4|Group_5/Group_4-5/g' -Ee 's/.*IVEX001_final_(.+)\.gff\t/\1\t/g' > tmp_IVEX002__005_000
## list all genomes/contigs
for f in $FGF; do awk -F"\t" '{print FILENAME"\t"$1}' $f; done | sort -Vu | sed -E 's/.*IVEX001_final_(.+)\.gff\t/\1\t/g' > tmp_IVEX002__005_001
## get count of each contig per category below
awk -F"\t" '{if($3 > $4) print $1"\t"$2}' tmp_IVEX002__005_000 | cat - tmp_IVEX002__005_001 | sort -V | uniq -c | sed -E 's/ +([0-9]+) .+/\1/g' | awk '{print $1-1}' > tmp_IVEX002__005_002_H
awk -F"\t" '{if($3 == $4) print $1"\t"$2}' tmp_IVEX002__005_000 | cat - tmp_IVEX002__005_001 | sort -V | uniq -c | sed -E 's/ +([0-9]+) .+/\1/g' | awk '{print $1-1}' > tmp_IVEX002__005_002_G
awk -F"\t" '{if($3 < $4) print $1"\t"$2}' tmp_IVEX002__005_000 | cat - tmp_IVEX002__005_001 | sort -V | uniq -c | sed -E 's/ +([0-9]+) .+/\1/g' | awk '{print $1-1}' > tmp_IVEX002__005_002_F
awk -F "\t" '{if($5 =="Group_1-3" && $3 > $4) print $1"\t"$2}' tmp_IVEX002__005_000 | cat - tmp_IVEX002__005_001 | sort -V | uniq -c | sed -E 's/ +([0-9]+) .+/\1/g' | awk '{print $1-1}' > tmp_IVEX002__005_002_E
awk -F "\t" '{if($5 =="Group_4-5" && $3 > $4) print $1"\t"$2}' tmp_IVEX002__005_000 | cat - tmp_IVEX002__005_001 | sort -V | uniq -c | sed -E 's/ +([0-9]+) .+/\1/g' | awk '{print $1-1}' > tmp_IVEX002__005_002_D
awk -F "\t" '{if($5 =="Not_assigned" && $3 > $4) print $1"\t"$2}' tmp_IVEX002__005_000 | cat - tmp_IVEX002__005_001 | sort -V | uniq -c | sed -E 's/ +([0-9]+) .+/\1/g' | awk '{print $1-1}' > tmp_IVEX002__005_002_C
awk -F "\t" '{if($5 =="none") print $1"\t"$2}' tmp_IVEX002__005_000 | cat - tmp_IVEX002__005_001 | sort -V | uniq -c | sed -E 's/ +([0-9]+) .+/\1/g' | awk '{print $1-1}' > tmp_IVEX002__005_002_B
cut -f 1,2 tmp_IVEX002__005_000 | cat - tmp_IVEX002__005_001 | sort -V | uniq -c | sed -E 's/ +([0-9]+) .+/\1/g' > tmp_IVEX002__005_002_A

## merge counts
paste -d ";" tmp_IVEX002__005_001 tmp_IVEX002__005_002_* | sed 's/\t/;/g' | sort -t ";" -k 1,1V -k 10,10nr -k 7,7nr  > IVEX002_final_005_contig_summary

#######################

###### keep only contigs with >2 viral genes
awk -F ";" '{if($3 >= 2) print $0}' IVEX002_final_005_contig_summary > tmp && mv tmp IVEX002_final_005_contig_summary
