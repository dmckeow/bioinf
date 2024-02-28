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

############### STEP 006 ############################
#### reformat gffs for manual inspection

#### get viral rbitscore minus cell rbitscore
for f in $FGF; do awk -F "\t|;" '{print $16"\t"$28}' $f | sed -Ee 's/vrbitscore=|rbitscore=//g' -Ee 's/NA//g' | awk -F "\t|;" '{if($1 ~ /[0-9]/) print $1-$2; else print ""}' > tmp_"$(basename $f .gff)"; done

#### alter type field based on viral minus cell scores; -1.0 to -0.2 cell, -0.2 to 0.2 virus or cell (even if no cellular hit), 0.2 to 1.0 virus
## this determines colour in Geneious
### an important on proteins with neither viral nor cellular hits; they are only VIRAL ORFans, because proteins without an initial viral hit were NOT blasted against NR - these proteins are not true ORFans because they have hits to the NR (cellular) database

for f in $FGF; do paste tmp_"$(basename $f .gff)" $f | awk -F"\t" 'BEGIN{OFS="\t"};{if($1 !~ /[0-9]/) gsub(/mRNA/,"mRNA_ORFan",$4)}1' | awk -F"\t" 'BEGIN{OFS="\t"};{if($1 ~ /[0-9]/ && $1 >= -0.2 && $1 <= 0.2) gsub(/mRNA/,"mRNA_cellular_or_viral",$4)}1' | awk -F"\t" 'BEGIN{OFS="\t"};{if($1 ~ /[0-9]/ && $1 > 0.2) gsub(/mRNA/,"mRNA_viral",$4)}1' | awk -F"\t" 'BEGIN{OFS="\t"};{if($1 ~ /[0-9]/ && $1 < -0.2) gsub(/mRNA/,"mRNA_cellular",$4)}1' | cut --complement -f 1 | sed -Ee 's/=none$|=NA$/=-/g' -Ee 's/=none;|=NA;/=-;/g' -Ee 's/=none\t|=NA\t/=-\t/g' > IVEX002_final_006_$(basename $f); done

## fix file names
for f in IVEX002_final_006_IVEX001_final_*; do mv "$f" "${f//IVEX002_final_006_IVEX001_final_/}"; done

### reduce to only those contigs kept in contig summary
sed '1,1d' IVEX002_final_005_contig_summary | awk -F ";" '{print $2 > "tmp_"$1}'
for f in *.gff; do grep -wF -f tmp_$(basename $f .gff) $f > tmp && mv tmp $f; done

#### add NCVOG group to 1-5 in brackets
for f in *.gff; do awk -F";" 'BEGIN{OFS=";"};{if($1 ~ /mRNA_cellular_or_viral|mRNA_cellular|mRNA_viral/) gsub(/$/,$13")",$14)}1' $f | awk -F";" 'BEGIN{OFS=";"};{if($1 ~ /mRNA_cellular_or_viral|mRNA_cellular|mRNA_viral/) gsub(/Name=-ncvog_group=Not_assigned)/,"Name=u",$14)}1' | awk -F";" 'BEGIN{OFS=";"};{if($1 ~ /mRNA_cellular_or_viral|mRNA_cellular|mRNA_viral/) gsub(/ncvog_group=Not_assigned)|ncvog_group=-)/,"",$14)}1' | awk -F";" 'BEGIN{OFS=";"};{if($1 ~ /mRNA_cellular_or_viral|mRNA_cellular|mRNA_viral/) gsub(/ncvog_group=Group_/,"(",$14)}1' > tmp && mv tmp $f; done

rm -f *.gff.gz; gzip -f *.gff

######### get reduced fasta for only contigs in gffs
for f in *.gff.gz; do zcat $f | cut -f 1 | sort -Vu | seqkit grep -n -f - "$fa"$(basename $f .gff.gz).fa.gz > $(basename $f .gff.gz).fa; done

rm -f *.fa.gz; gzip -f *.fa

## file to quickly find viral clusters
for f in $FGF; do sed -Ee 's/vrbitscore=|ncvog_group=|vlineage=|ncvog_legend=//g' -e 's/rbitscore=//g' $f | awk -F "\t|;" -v f=$f 'BEGIN{OFS=";"};{if($3 =="mRNA" && $16 !~ "NA" && $28 !~ "NA" && $10 ~ "vsalltitles=" && $10 !~ "vsalltitles=none" && $16 >= 0.2) print f,$1,$9,$4,$5,$16,$28,$10,$23}'; done | sed -Ee 's/^.*\/*IVEX001_final_(.+)\.gff;/\1;/g' -e 's/ID=|vsalltitles=//g' | awk -F";" 'BEGIN{OFS=";"};{if($1 == b && $2 == c && $4-a < 30000) print $1,$2,$3,$4,$5,$6-$7,$8,$9; else print "#####"$1,$2,$3,$4,$5,$6-$7,$8,$9} {a=$5} {b=$1} {c=$2}' > IVEX002_final_006_clusters

rm -f tmp*
