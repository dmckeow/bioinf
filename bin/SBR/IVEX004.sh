#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition long
#SBATCH --cpus-per-task 8
#SBATCH --mem 64GB

###### MANUAL STEPS AND PREREQUISITES ######
## The purpose of this script is to BLASTP all proteins with NO viral hits against NR to identify host _proteins
## It is intended to run on a small set of contigs with known viral regions (take a long time if whole genome BLASTed)
## The main output is a modified version of IVEX002.sh gff file, with cell-only annotations added
## After running this, use the output gffs to generate figures from IVEX_cgview.sh and to view in Geneious
## 1. Ran scripts IVEX000-002 and have a small set of contigs of interest
## 2. list of contigs of interest that match names in gff

### run from project folder (which contains archives input script tmp finalresult)

####################### INPUTS #####################################
###### load software ######
module load seqkit/0.14.0; module load diamond/2.0.9

##### command line arguments required #####
## $1 = IVEX002 gff.gz
## $2 = contig list
## $3 = 001 (change 001 to 002 once all genomes run through STEP 001)
## $4 = EVE info file containing:
## genome;contig;EVE;EVE_start;EVE_end;context;contig_size;eve_size
## Chordaria-linearis;C-linearis_contig13;13.1;2214000;2562000;HVH;3542593;348000

##### input files for analysis #####
pwd=$(pwd) ## get path to project folder
FAA=$pwd"/input/IVEX000/*$(basename "$1" .gff.gz)*.faa" ## the protein fastas prepared by IVEX000.sh
FA=$pwd"/input/IVEX000_genomes/$(basename "$1" .gff.gz).fa.gz"
FGF=$pwd"/finalresult/IVEX001/IVEX001_final_*.gff"
evelist=$pwd"/$4"

###### input file for databases ######
DB="/shared/projects/phaeoexplorer_virus/db/" ## you may need to change this if large database files are too large
NR=$DB"nr.dmnd" ## see STEP 000 of script IVEX001.sh to generate this database file
NRFA="/db/nr/current/fasta/nr.fsa" ## this is located on server

##### NCBI taxonomy files #####
A2T="/db/accession2taxid/current/flat/prot.accession2taxid" ## this is located on server
NODE="/db/taxonomy/current/flat/nodes.dmp" ## this is located on server
NAME="/db/taxonomy/current/flat/names.dmp" ## this is located on server

##### NCVOG defintion files #####
NCVKEY=$pwd"/input/IVEX000_NCVOG" ## made by script IVEX000_NCVOG.sh

####################### SCRIPT #####################################
cd $pwd"/finalresult/IVEX004"

################ STEP 001 ##############
if [[ "$3" =~ "001" ]]; then

### get gff for all viral ORFan genes from contigs of interest
zcat $1 | awk '/\tmRNA_ORFan\t/' | grep -wF -f $2 - | sort -t $'\t' -k 9,9V > tmp1_IVEX004_$(basename "$1" .gff.gz).gff
### get protein fasta for all viral ORFan genes from contigs of interest
awk -F "\t|;" '{print $9}' tmp1_IVEX004_$(basename "$1" .gff.gz).gff | sed 's/ID=//g' | seqkit grep -n -f - $FAA > tmp1_IVEX004_$(basename "$1" .gff.gz).faa

##### diamond all viral ORFan genes from contigs of interest vs NR, excluding viruses
diamond blastp -k 10 -e 1e-5 --unal 1 --taxon-exclude 10239 -d $NR -q tmp1_IVEX004_$(basename "$1" .gff.gz).faa -o IVEX004_$(basename "$1" .gff.gz).NR.blast -f 6 qseqid salltitles evalue bitscore staxids

### reformat blast results for gff
sort -t $'\t' -k 1,1V -k 4,4nr IVEX004_$(basename "$1" .gff.gz).NR.blast | awk -F"\t" '!a[$1]++' | awk -F"\t" 'BEGIN{OFS="\t"};{if($2 !="\*") gsub(/[^a-zA-Z0-9]/,"_",$2)}1' | sed 's/\t/;/g' | sed -Ee 's/_+/_/g' -Ee 's/^_|_$//g' -Ee 's/;_|_;/;/g' | awk -F";" 'BEGIN{OFS=";"};{if($2 =="\*") print $1,"-","-","-","-"; else print $0}' | awk -F";" 'BEGIN{OFS=";"};{print "ID="$1,"salltitles="$2,"evalue="$3,"bitscore="$4,"staxids="$5,"rbitscore=-"}' | sed -Ee 's/=;/=-;/g' | sort -t $'\t' -k 1,1V > tmp2_$(basename "$1" .gff.gz)

### merge new BLAST results with corresponding gff line (inserted after field 23) and change to type mRNA_cellular if hit present for protein
sed 's/;/\t/g' tmp1_IVEX004_$(basename "$1" .gff.gz).gff | cut --complement -f 24-28 | cut -f 1-23 | paste - tmp2_$(basename "$1" .gff.gz) > tmp3_$(basename "$1" .gff.gz)
sed 's/;/\t/g' tmp1_IVEX004_$(basename "$1" .gff.gz).gff | cut --complement -f 24-28 | cut --complement -f 1-23 | paste tmp3_$(basename "$1" .gff.gz) - | sed 's/\t/;/9g' | cut --complement -d ";" -f 16 | awk -F"\t" 'BEGIN{OFS="\t"};{if($0 ~ /;bitscore=[0-9]/) gsub(/.+/,"mRNA_cellular",$3)}1' | sed -Ee 's/;$//g' > tmp_IVEX004_$(basename "$1" .gff.gz).gff

##### remove original gff lines of BLASTed proteins from input gff
### get gff lines for contigs of interest, except the "ORFan" ones, then add in "ORFan" lines with new BLAST results, making the final gff
zcat $1 | grep -wF -f $2 - | sed '/\tmRNA_ORFan\t/d' | cat - tmp_IVEX004_$(basename "$1" .gff.gz).gff | sort -t $'\t' -k 9,9V > IVEX004_$(basename "$1" .gff.gz).gff

### get reduced contig assembly fro convenience
seqkit grep -n -f $2 $FA > IVEX004_$(basename "$1" .gff.gz).fa

rm -f tmp*
fi

################ STEP 002 ##############
#### make EVE heatmap
if [[ "$3" =~ "002" ]]; then

#for f in $FGF; do awk -F"\t" -v a=$(basename "$f" .gff) '{if($3 =="mRNA") print a"\t"$0}' $f; done | sed 's/IVEX001_final_//g' > test_tmp2_1

#split -d -e -l 1 $evelist tmp_evelist_

#for f in tmp_evelist_*; do
  #awk -F "\t" -v a=$(cut -d ";" -f 1 $f) -v b=$(cut -d ";" -f 2 $f) '{if($0 ~ a"\t"b"\t") print $0}' test_tmp2_1 | awk -F "\t" -v c=$(cut -d ";" -f 4 $f) -v d=$(cut -d ";" -f 5 $f) -v e=$(cut -d ";" -f 3 $f) 'BEGIN{OFS="\t"};{if(c <= $5 && d >= $6) print $1,e,$3,$4,$5,$6,$7,$8,$9,$10}' > 2"$f" &
#done
#wait
#cat 2tmp_evelist_* > test_tmp2_2
rm -f *tmp_evelist_*
rm -f test_tmp2_1

###### NCVOG count per EVE (absent NCVOGs (except group1-5) and EVEs with no NCVOGs not included)
#### for viewing in excel, etc

## list all NCVOGs present in genomes (with viral rbitscore > cellular rbitscore and 0.2)
sed -e 's/vrbitscore=//g' -e 's/rbitscore=//g' $FGF | awk -F "\t|;" 'BEGIN{OFS=";"};{if($3 =="mRNA" && $16 !~ "NA" && $28 !~ "NA" && $18 ~ /ncvog_id=[0-9]+/ && $16 > $28 && $16 >= 0.2) print $18,$21,$22,$23}' | sort -Vu > IVEX004__001_000
## find all EVEs with NCVOG hits (with viral rbitscore > cellular rbitscore), split into tmp files by EVE
sed -e 's/vrbitscore=//g' -e 's/rbitscore=//g' test_tmp2_2 | awk -F "\t|;" 'BEGIN{OFS=";"};{if($4 =="mRNA" && $17 !~ "NA" && $29 !~ "NA" && $19 ~ /ncvog_id=[0-9]+/ && $17 > $29 && $17 >= 0.2) print $1,$2,$19,$22,$23,$24 > "tmp_"$1}'
### list NCVOGs from Groups1-5 and/or EsV-1/FsV-158
awk -F"\t" '{if($0 ~ "ncvog_group=Group_" || $0 ~ "Ectocarpus_siliculosus_virus_1" || $0 ~ "Feldmannia_species_virus") print $5";"$7";"$8";"$9}' $NCVKEY | sort -Vu > IVEX004__001_001
## reduce to only NCVOGs from Groups1-5 and/or EsV-1/FsV-158 and present in your genomes
cat IVEX004__001_000 IVEX004__001_001 | sort -Vu > tmp && mv tmp IVEX004__001_001
### merge genome and EVE name, in case names are ambiguous
## fix merge of genome and EVE name
for f in tmp_*; do sed -i -E 's/^(.+);(.+)(;ncvog_id=.+$)/\1__\2\3/g' $f; done
## remove name column and rename fiel per EVE to include genome name
for f in tmp_*; do awk -F ";" 'BEGIN{OFS=";"};{print $2,$3,$4,$5 > "2tmp_"$1}' $f; done
#### prep counts per NCVOG, each file EVE
for f in 2tmp_*; do cat IVEX004__001_001 $f | sort -t";" -k 2,2V -k 1,1V | uniq -c | sed -E 's/ +([0-9]+) /\1\t/g' | awk -F"\t" -v f=$f '{print f"\t"$1-1"\t"$2}' > tmp && mv tmp $f; done
rm -f tmp_*; for f in 2tmp_*; do mv "$f" "${f//2tmp_/tmp_}"; done

#### prep counts for EsV-1 and FsV-158 references
awk -F "\t" 'BEGIN{OFS=";"};{if($1 ~"Ectocarpus_siliculosus_virus_1") print $5,$7,$8,$9}' $NCVKEY | cat - IVEX004__001_001 | sort -t";" -k 2,2V -k 1,1V | uniq -c | sed -Ee 's/ +([0-9]+) /\1\t/g' | awk -F"\t" '{print "EsV-1\t"$1-1"\t"$2}' > tmp_EsV-1
awk -F "\t" 'BEGIN{OFS=";"};{if($1 ~"Feldmannia_species_virus") print $5,$7,$8,$9}' $NCVKEY | cat - IVEX004__001_001 | sort -t";" -k 2,2V -k 1,1V | uniq -c | sed -Ee 's/ +([0-9]+) /\1\t/g' | awk -F"\t" '{print "FsV-158\t"$1-1"\t"$2}' > tmp_FsV-158
## get first column for table (ncvog label)
sort -t";" -k 2,2V -k 1,1V IVEX004__001_001 | sed -e 's/ncvog_id=/_/g' -Ee 's/ncvog_group=Group_|ncvog_legend=|ncvog_abbrev=//g' -e 's/ncvog_group=Not_assigned/na/g' > IVEX004__001_002

##### get header for table (EVE label)
cat tmp_* | cut -f 1 | awk '!a[$0]++' | sed -e 's/$/;/g' -Ee 's/^2tmp_|^tmp_//g' | tr -d '\n' | sed 's/;$/\n/g' | sed 's/^/ncvog_id;ncvog_group;ncvog_abbrev;ncvog_legend;/g' > IVEX004__001_003

#### prepare heatmap for R (for manually picking out specific EVEs later)
cat tmp_* | sed 's/^2tmp_/tmp_/g' | sort -t ";" -k 2,2Vr | sed -e 's/;/\t/g' -e '/ncvog_group=Not_assigned/d' -Ee 's/ncvog_group=Group_|ncvog_id=|ncvog_legend=|ncvog_abbrev=//g' -e 's/\t0\t/\t\t/g' -e 's/^tmp_//g' -e 's/^EsV-1\t/EsV-1__EsV-1\t/g' -e 's/^tmp_//g' -e 's/^FsV-158\t/FsV-158__FsV-158\t/g' > tmp2_IVEX004_final_001_ncvog_count_EVEs_heatmap

#### remove EVEs with no NCVOGs and split by genome
awk -F"\t" '{if($2 !="") print $1}' tmp2_IVEX004_final_001_ncvog_count_EVEs_heatmap | sort -Vu | grep -wF -f - tmp2_IVEX004_final_001_ncvog_count_EVEs_heatmap | awk -F "__" '{print $0 > "IVEX004_final_001_"$1"_ncvog_count_EVEs_heatmap"}'

### rename reference genome heatmaps
mv IVEX004_final_001_EsV-1_ncvog_count_EVEs_heatmap tmp3_EsV-1
mv IVEX004_final_001_FsV-158_ncvog_count_EVEs_heatmap tmp3_FsV-158

## remerge with reference and add header
for f in IVEX004_final_001_*_ncvog_count_EVEs_heatmap; do cat tmp3_* $f | sed -z 's/^/EVE\tcount\tncvog\tgroup\tabbrev\tlegend\n/1' | sed -E 's/^.+__//1' > tmp && mv tmp $f; done
for f in IVEX004_final_001_*_ncvog_count_EVEs_heatmap; do awk -F "\t" -v f=$f '{print $0"\t"f}' $f | sed -Ee 's/\tIVEX004_final_001_(.+)_ncvog_count_EVEs_heatmap/\t\1/g' -Ee 's/(\tlegend\t).+/\1genome/g' > tmp && mv tmp $f; done

## add in extra EVE info from eve_list
awk -F";" '{print "s/"$1";"$3"___/"$6";"$7";"$8"/g"}' $evelist > tmp_fix;
for f in IVEX004_final_001_*_ncvog_count_EVEs_heatmap; do sed '1,1d' $f | awk -F"\t" '{print $7";"$1"___"}' | sed -f tmp_fix - | sed 's/.*___/;;/g' | sed -z 's/^/context;contig_size;EVE_size\n/1' | paste -d ";" $f - | sed 's/;/\t/g' > tmp && mv tmp $f; done

#### merge into final file for excel viewing
for f in tmp_*; do cut -f 2 $f > tmp; cp tmp $f; done;
find . -maxdepth 1 -type f -name tmp_\* | sort | split -l 1000 -d - tmp_lists;
for f in tmp_lists*; do paste -d ";" $(cat $f) > tmp_merge_$(basename $f); done; paste -d ";" IVEX004__001_002 tmp_merge* > IVEX004_final_001_ncvog_count_EVEs

cat IVEX004__001_003 IVEX004_final_001_ncvog_count_EVEs > tmp && mv tmp IVEX004_final_001_ncvog_count_EVEs

rm -f tmp*
### these steps swap the rows and columns (transpoase) IVEX004_final_001_ncvog_count_EVEs and add explanation header
split -d -l 1 IVEX004_final_001_ncvog_count_EVEs tmp_
sed -i 's/;/\n/g' tmp_*;
find . -maxdepth 1 -type f -name tmp_\* | sort | split -l 1000 -d - tmp_lists;
for f in tmp_lists*; do paste -d ";" $(cat $f) > tmp_merge_$(basename $f); done; paste -d ";" tmp_merge* > tmp; mv tmp IVEX004_final_001_ncvog_count_EVEs;

sed -i -z 's/^/## no. of copies of core Nucleocytoviricota genes (NCVOGs) per EVE (with viral rbitscore > cellular rbitscore and 0.2). Group refers to level of gene conservation, 1 being the most conserved core genes.\n/g' IVEX004_final_001_ncvog_count_EVEs
rm -f tmp*

fi
############### STEP 003 ############################
if [[ "$3" =~ "003" ]]; then
##### counts of virus EVEs checklist
## get genome, EVE rbitscores, and virus hit type key file for all EVEs
sed -Ee 's/vrbitscore=|ncvog_group=|vlineage=//g' -e 's/rbitscore=//g' test_tmp2_2 | awk -F "\t|;" 'BEGIN{OFS="\t"};{if($4 =="mRNA" && $17 !~ "NA" && $29 !~ "NA" && $11 ~ "vsalltitles=" && $11 !~ "vsalltitles=none" && $17 >= 0.2) print $1,$2,$17,$29,$22,$18}' | sed -Ee 's/Group_1|Group_2|Group_3/Group_1-3/g' -Ee 's/Group_4|Group_5/Group_4-5/g' > tmp_IVEX004__005_000
## list all genomes/EVEs
cut -f 1,2 test_tmp2_2 | sort -Vu > tmp_IVEX004__005_001
## get count of each EVE per category below
awk -F"\t" '{if($3 > $4) print $1"\t"$2}' tmp_IVEX004__005_000 | cat - tmp_IVEX004__005_001 | sort -V | uniq -c | sed -E 's/ +([0-9]+) .+/\1/g' | awk '{print $1-1}' > tmp_IVEX004__005_002_H
awk -F"\t" '{if($3 == $4) print $1"\t"$2}' tmp_IVEX004__005_000 | cat - tmp_IVEX004__005_001 | sort -V | uniq -c | sed -E 's/ +([0-9]+) .+/\1/g' | awk '{print $1-1}' > tmp_IVEX004__005_002_G
awk -F"\t" '{if($3 < $4) print $1"\t"$2}' tmp_IVEX004__005_000 | cat - tmp_IVEX004__005_001 | sort -V | uniq -c | sed -E 's/ +([0-9]+) .+/\1/g' | awk '{print $1-1}' > tmp_IVEX004__005_002_F
awk -F "\t" '{if($5 =="Group_1-3" && $3 > $4) print $1"\t"$2}' tmp_IVEX004__005_000 | cat - tmp_IVEX004__005_001 | sort -V | uniq -c | sed -E 's/ +([0-9]+) .+/\1/g' | awk '{print $1-1}' > tmp_IVEX004__005_002_E
awk -F "\t" '{if($5 =="Group_4-5" && $3 > $4) print $1"\t"$2}' tmp_IVEX004__005_000 | cat - tmp_IVEX004__005_001 | sort -V | uniq -c | sed -E 's/ +([0-9]+) .+/\1/g' | awk '{print $1-1}' > tmp_IVEX004__005_002_D
awk -F "\t" '{if($5 =="Not_assigned" && $3 > $4) print $1"\t"$2}' tmp_IVEX004__005_000 | cat - tmp_IVEX004__005_001 | sort -V | uniq -c | sed -E 's/ +([0-9]+) .+/\1/g' | awk '{print $1-1}' > tmp_IVEX004__005_002_C
awk -F "\t" '{if($5 =="none") print $1"\t"$2}' tmp_IVEX004__005_000 | cat - tmp_IVEX004__005_001 | sort -V | uniq -c | sed -E 's/ +([0-9]+) .+/\1/g' | awk '{print $1-1}' > tmp_IVEX004__005_002_B
cut -f 1,2 tmp_IVEX004__005_000 | cat - tmp_IVEX004__005_001 | sort -V | uniq -c | sed -E 's/ +([0-9]+) .+/\1/g' > tmp_IVEX004__005_002_A

## merge counts
paste -d ";" tmp_IVEX004__005_001 tmp_IVEX004__005_002_* | sed 's/\t/;/g' | sort -t ";" -k 1,1V -k 10,10nr -k 7,7nr  > IVEX004_final_005_EVEs_summary

#######################

###### keep only EVEs with >2 viral genes
awk -F ";" '{if($3 >= 2) print $0}' IVEX004_final_005_EVEs_summary | sort -t ";" -k 1,1V -k 11,11nr -k 8,8nr | sed -z 's/^/genome;EVE;total_all_virus_genes;non-ncvog(any_virus);ncvog_no_group;ncvog_group_4-5;ncvog_group_1-3;viral_rbitscore_lesser;rbitscores_equal;viral_rbitscore_greater\n/1' > tmp && mv tmp IVEX004_final_005_EVEs_summary

rm -f tmp*
fi
