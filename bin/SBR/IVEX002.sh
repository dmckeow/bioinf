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
cd $pwd"/finalresult/IVEX002_SEP22"

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
###### NCVOG count per genome (absent NCVOGs (except group1-5) and genomes with no NCVOGs not included)
#### for excel, etc and heatmap in R

## same IVEX002__001_000 used
## find all ncvog hits (with viral rbitscore > cellular rbitscore), split into tmp files by genome
for f in $FGF; do sed -e 's/vrbitscore=//g' -e 's/rbitscore=//g' $f | awk -F "\t|;" 'BEGIN{OFS=";"};{if($3 =="mRNA" && $16 !~ "NA" && $28 !~ "NA" && $18 ~ /ncvog_id=[0-9]+/ && $16 > $28 && $16 >= 0.2) print $18,$21,$22,$23}' > tmp_$(basename $f .gff); done
## same IVEX002__001_001 used
#### prep counts per NCVOG, each file genome
for f in tmp_*; do cat IVEX002__001_001 $f | sort -t";" -k 2,2V -k 1,1V | uniq -c | sed -E 's/ +([0-9]+) /\1\t/g' | awk -F"\t" -v f=$f '{print f"\t"$1-1"\t"$2}' > tmp && mv tmp $f; done

#### prep counts for EsV-1 and FsV-158 references (same as in STEP 001)
awk -F "\t" 'BEGIN{OFS=";"};{if($1 ~"Ectocarpus_siliculosus_virus_1") print $5,$7,$8,$9}' $NCVKEY | cat - IVEX002__001_001 | sort -t";" -k 2,2V -k 1,1V | uniq -c | sed -Ee 's/ +([0-9]+) /\1\t/g' | awk -F"\t" '{print "EsV-1\t"$1-1"\t"$2}' > tmp_EsV-1
awk -F "\t" 'BEGIN{OFS=";"};{if($1 ~"Feldmannia_species_virus") print $5,$7,$8,$9}' $NCVKEY | cat - IVEX002__001_001 | sort -t";" -k 2,2V -k 1,1V | uniq -c | sed -Ee 's/ +([0-9]+) /\1\t/g' | awk -F"\t" '{print "FsV-158\t"$1-1"\t"$2}' > tmp_FsV-158

## get first column for table (ncvog label)
sort -t";" -k 2,2V -k 1,1V IVEX002__001_001 | sed -e 's/ncvog_id=/_/g' -Ee 's/ncvog_group=Group_|ncvog_legend=|ncvog_abbrev=//g' -e 's/ncvog_group=Not_assigned/na/g' > IVEX002__002_002
## get header for table (genome label)
cat tmp_* | cut -f 1 | awk '!a[$0]++' | sed 's/$/;/g' | tr -d '\n' | sed 's/;$/\n/g' | sed 's/^/ncvog_id;ncvog_group;ncvog_abbrev;ncvog_legend;/g' | sed 's/tmp_IVEX001_final_//g' > IVEX002__002_003

#### prepare heatmap for R (for all genomes)
cat tmp_* | sort -t ";" -k 2,2Vr | sed -e 's/;/\t/g' -e '/ncvog_group=Not_assigned/d' -Ee 's/ncvog_group=Group_|ncvog_id=|ncvog_abbrev=|ncvog_legend=//g' -e 's/\t0\t/\t\t/g' -e 's/tmp_IVEX001_final_//g' | sed -z 's/^/genome\tcount\tncvog\tgroup\tabbrev\tlegend\n/1' > IVEX002_final_002_ncvog_count_genomes_heatmap

#### merge into final file for excel viewing
for f in tmp_*; do cut -f 2 $f > tmp && mv tmp $f; done
find . -maxdepth 1 -type f -name tmp_\* | sort | split -l 1000 -d - tmp_lists;
for f in tmp_lists*; do paste -d ";" $(cat $f) > tmp_merge_$(basename $f); done; paste -d ";" IVEX002__002_002 tmp_merge* > IVEX002_final_002_ncvog_count_genomes
cat IVEX002__002_003 IVEX002_final_002_ncvog_count_genomes > tmp && mv tmp IVEX002_final_002_ncvog_count_genomes
rm -f tmp*
### these steps swap the rows and columns (transpose) IVEX002_final_001_ncvog_count_genomes and add explanation header
split -d -l 1 IVEX002_final_002_ncvog_count_genomes tmp_
sed -i 's/;/\n/g' tmp_*;
find . -maxdepth 1 -type f -name tmp_\* | sort | split -l 1000 -d - tmp_lists;
for f in tmp_lists*; do paste -d ";" $(cat $f) > tmp_merge_$(basename $f); done; paste -d ";" tmp_merge* > tmp; mv tmp IVEX002_final_002_ncvog_count_genomes;
sed -i -z 's/^/## no. of copies of core Nucleocytoviricota genes (NCVOGs) per genome (with viral rbitscore > cellular rbitscore and 0.2). Group refers to level of gene conservation, 1 being the most conserved core genes.\n/g' IVEX002_final_002_ncvog_count_genomes
rm -f tmp*

############### STEP 003 ############################
#### make heatmap of viral taxa counts per genome

#### per genome, get viral taxonomy lineage for every viral gene (with viral rbitscore > 0.2 cellular rbitscore 0)
for f in $FGF; do sed -Ee 's/vrbitscore=|vlineage=//g' -e 's/rbitscore=//g' $f | awk -F "\t|;" 'BEGIN{OFS=";"};{if($3 =="mRNA" && $16 !~ "NA" && $28 !~ "NA" && $10 ~ "vsalltitles=" && $10 !~ "vsalltitles=none" && $16 > $28 && $16 >= 0.2 && $17 != "none") print $17}' > tmp_$(basename $f .gff); done

###### reformat to family,genus or equivalent category if possible
for f in tmp_*; do sed -Ee 's/(.*-)([A-Za-z]+virus)(.*-.*)/\2-\1\2\3/g' -Ee 's/(.*-)([A-Za-z]+viridae)(.*-.*)/\2-\1\2\3/g' -Ee 's/(.*-)([A-Za-z]+viricota)(.*-.*)/\2-\1\2\3/g' $f | cut -d "-" -f 1-3 | sed -E 's/-Viruses|Viruses-//g' | awk -F "-" '{if($2 =="") print "unclassified-"$0; else print $0}'| awk -F "-" '{if($3 =="") print "unclassified-"$0; else print $0}' | awk -F "-" '{if($0 ~ "viricota" && $1 =="unclassified") print $2"-"$3"-"$1; else print $0}' > tmp && mv tmp $f; done

## non-redundant list of all viral taxa present in genomes
cat tmp_* | sort -Vu > IVEX002__003_000

#### prepare v taxa counts per genome
for f in tmp_*; do cat IVEX002__003_000 $f | sort -V | uniq -c | sed -E 's/ +([0-9]+) /\1\t/g' | awk -F"\t" -v f=$f '{print f"\t"$1-1"\t"$2}' > tmp && mv tmp $f; done

#### prepare heatmap for R (for all genomes), hits for virus taxa count (convert to %)
for f in tmp_*; do total=$(awk -F "\t" '{sum+=$2;} END{print sum;}' $f); awk -v total=$total -v f=$f '{if(total > 0) print f"\t"$2"\t"($2/total)*100"\t"$3; else print f"\t0\t0\t"$3}' $f | sed 's/tmp_IVEX001_final_//g' | sed 's/\t0\t0\t/\t\t\t/g' ; done | sort -t $'\t' -k 3,3Vr | sed -z 's/^/genome\tcount\tpercent\ttaxa\n/g' > IVEX002_final_003_virus_taxa_percent_count_genomes_heatmap


#### separate the 2 taxa labels into different columns
awk -F"\t" 'BEGIN{OFS="\t"}; {gsub("-","\t",$4)}1' IVEX002_final_003_virus_taxa_percent_count_genomes_heatmap > tmp && mv tmp IVEX002_final_003_virus_taxa_percent_count_genomes_heatmap
sed -i 's/\ttaxa$/\ttaxa1\ttaxa2\ttaxa3/g' IVEX002_final_003_virus_taxa_percent_count_genomes_heatmap
rm -f tmp*

############### STEP 004 ############################
##### plot scattergraphs of rbitscore virus vs cellular per genome

for f in $FGF; do sed -Ee 's/vrbitscore=|ncvog_group=|vlineage=//g' -e 's/rbitscore=//g' $f | awk -F "\t|;" -v f=$f 'BEGIN{OFS="\t"};{if($3 =="mRNA" && $16 !~ "NA" && $28 !~ "NA" && $10 ~ "vsalltitles=" && $10 !~ "vsalltitles=none") print f,$16,$28,$21}' ; done | sed -E 's/.*IVEX001_final_(.+).gff/\1/g' | sed -Ee 's/Group_1|Group_2|Group_3/NCVOG_groups_1-3/g' -Ee 's/Group_4|Group_5/Not_assigned/g' | awk -F "\t" 'BEGIN{OFS="\t"};{if($4 =="none") print $1,$2,$3,"non-NCV"; else print $0}' | awk -F "\t" 'BEGIN{OFS="\t"};{if($4 =="Not_assigned") print $1,$2,$3,"NCV"; else print $0}' | sort -t $'\t' -k 1,4V | sed -Ee 's/(\tnon-NCV)$/\1\t1/g' -Ee 's/(\tNCV)$/\1\t2/g' -Ee 's/(\tNCVOG_groups_1-3)$/\1\t3/g' | sed -z 's/^/genome\tvrbitscore\trbitscore\tncvog\tsort\n/' > IVEX002_final_004_rbitscore_scatterplot
sed '1,1d' IVEX002_final_004_rbitscore_scatterplot | awk -F"\t" '{print $0 > "IVEX002_final_004_rbitscore_scatterplot_"$1}'
sed -i -z 's/^/genome\tvrbitscore\trbitscore\tncvog\tsort\n/' IVEX002_final_004_rbitscore_scatterplot_*

rm -f tmp*

############### STEP 005 ############################
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

################# THIS SECTION IS OBSOLETE - USE IT IF YOU WANT TO INCLUDE VIRSORTER
#### list contigs which virsorter called viral
##for f in "$VS"IVEX001.5_final_*.gff; do
    ##awk -F "\t" -v f=$f '{print f";"$1}' $f | sort -Vu | sed -E 's/.*IVEX001\.5_final_(.+)\.gff(;.+)$/s\/\1\2_____\/virsorter\/g/g'
##done > IVEX002_virsorter_list

#### flag contigs with virsorter in contig summary
##cut -d ";" -f 1,2 IVEX002_final_005_contig_summary | sed 's/$/_____/g' > tmp
##split -d --number=l/20 tmp tmp_vs_
## parallel jobs
##for f in tmp_vs_*; do
  ##sed -f IVEX002_virsorter_list $f > tmp3_"$f" &
##done
##wait

##cat tmp3_tmp_vs_* | sed 's/.*_____$//g' | paste -d ";" IVEX002_final_005_contig_summary - > tmp && mv tmp IVEX002_final_005_contig_summary
#######################

###### keep only contigs with >2 viral genes (OR with virsorter flag if previous obsolete section used
awk -F ";" '{if($3 >= 2 || $0 ~ "virsorter") print $0}' IVEX002_final_005_contig_summary > tmp && mv tmp IVEX002_final_005_contig_summary

##### get contig sizes
##### summary file is counts of genes with virus hits only and rbitscore => 0.3\n
cut -d ";" -f 1,2 IVEX002_final_005_contig_summary | awk -F ";" '{print $2 > $1".tmp_xx1"}';
for f in *.tmp_xx1; do seqkit grep -n -f $f "$fa"$(basename $f .tmp_xx1).fa.gz > $(basename $f .tmp_xx1).tmp_xx2; done
for f in *.tmp_xx2; do csplit -z -s -n 12 -f "$f".tmp_xx3_ $f /\>/ {*}; done
cat *.tmp_xx2 | csplit -s -n 12 -f tmp_xx3_ - /\>/ {*}
for f in *.tmp_xx2.tmp_xx3_*; do ls $f; done > tmp_xx4; sed -i -E 's/\.tmp_xx2\.tmp_xx3_[0-9]+//g' tmp_xx4
for f in *.tmp_xx2.tmp_xx3_*; do head -1 $f; done > tmp_xx5; sed -i 's/>//g' tmp_xx5
for f in *.tmp_xx2.tmp_xx3_*; do seqkit stat -T $f; done > tmp_xx6; cut -f 8 tmp_xx6 | sed '/max_len/d' > 2tmp_xx6 && mv 2tmp_xx6 tmp_xx6
paste -d ";" tmp_xx4 tmp_xx5 tmp_xx6 | sort -t ";" -k 1,2V > tmp_xx7
sort -t ";" -k 1,2V IVEX002_final_005_contig_summary | cut -d ";" --complement -f 1,2 > tmp_xx8
sort -t ";" -k 1,2V IVEX002_final_005_contig_summary | cut -d ";" -f 1,2 | paste -d ";" - tmp_xx7 tmp_xx8 > tmp_xx9
awk -F ";" '{if($2 != $4) print $0}' tmp_xx9 > check_IVEX002_final_005_contig_summary
cut -d ";" --complement -f 1,2 tmp_xx9 | sort -t ";" -k 1,1V -k 11,11nr -k 8,8nr | sed -z 's/^/genome;contig;contig_size;total_all_virus_genes;non-ncvog(any_virus);ncvog_no_group;ncvog_group_4-5;ncvog_group_1-3;viral_rbitscore_lesser;rbitscores_equal;viral_rbitscore_greater;note;context;mummer_alignment\n/1' > IVEX002_final_005_contig_summary

for f in *tmp_xx*; do rm -f $f; done; rm -f tmp*

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
