#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition long
#SBATCH --cpus-per-task 16
#SBATCH --mem 64GB

###### MANUAL STEPS AND PREREQUISITES ######
## 1. Must have run all three IVEX000 scripts
## 2. If you already have a diamond db (.dmnd) of RVDB and NCBI's nr, built WITH $NODE, $NAME, and $A2T, deactivate STEP 000 with #
## 3. generate inputs $VT using Taxonkit (not on cluster currently): list all taxids from viruses (10239) and brown algae (2870), then get full lineage for each taxid [taxid]TAB[kingdom;phylum;etc] - see Taxonkit guide - very easy) use taxonkit lineage [file listing taxids of all virus taxa on NCBI] | sed -e 's/[^A-Za-z0-9\t;]/_/g' -Ee 's/_+/_/g' -Ee 's/^_|_$//g'
## 4. like IVEX000.sh, this script runs on genomes in groups (which were set in script IVEX000.sh); run this script for each rungroup - it can be run simultaneously for all groups

### run from project folder (which contains archives input script tmp finalresult)

####################### INPUTS #####################################
###### load software ######
module load seqkit/0.14.0; module load diamond/2.0.9

##### command line arguments required #####
## $1 = rungroup number e.g. 1

##### input files for analysis #####
pwd=$(pwd) ## get path to project folder
DB="/shared/projects/phaeoexplorer_virus/db/" ## you may need to change this if large database files are too large
FA=$pwd"/input/IVEX000/gms"$1"_*.faa" ## the protein fastas prepared by IVEX000.sh
GF=$pwd"/input/IVEX000/gms"$1"_*.gff" ## the annotations prepared by IVEX000.sh

###### input file for databases ######
RVFA=$DB"U-RVDBv21.0-prot.fasta" ## download latest fasta from https://rvdb-prot.pasteur.fr/ and unzip
RVDB=$DB"U-RVDBv21.0-prot.dmnd" ## see STEP 000 to generate this database file

NCVFA=$DB"NCVOG.fa" ## see STEP 000 to download this file
NCVOG=$DB"NCVOG.dmnd" ## see STEP 000 to generate this database file

NRFA="/db/nr/current/fasta/nr.fsa" ## this is located on server
NR=$DB"nr.dmnd" ## see STEP 000 to generate this database file

##### NCVOG defintion files #####
NCVKEY=$pwd"/input/IVEX000_NCVOG" ## made by script IVEX000_NCVOG.sh

##### NCBI taxonomy files #####
A2T="/db/accession2taxid/current/flat/prot.accession2taxid" ## this is located on server
NODE="/db/taxonomy/current/flat/nodes.dmp" ## this is located on server
NAME="/db/taxonomy/current/flat/names.dmp" ## this is located on server
VT=$DB"taxonkit_viruses_full_lineage" ## see MANUAL STEPS AND PREREQUISITES point #3


####################### SCRIPT #####################################

################ STEP 000 ##############
##### to create and setup NCVOG database:
##cd $DB; rm -f NCVOG.fa.tar.gz* NCVOG.fa; wget https://ftp.ncbi.nih.gov/pub/wolf/COGs/NCVOG/NCVOG.fa.tar.gz
##tar -xf NCVOG.fa.tar.gz; rm -f NCVOG.fa.tar.gz;
##diamond makedb --in NCVOG.fa --db NCVOG.dmnd;
##cd $pwd

##### create diamond RVDB protein database WITH taxonomy info (names must be returned to NCBI-format first)
### takes 30 min to make - only make if required - ENSURE all fasta and taxid inputs are some same update/version
##cd $DB; sed -i -E 's/>.+\|.+\|(.+)\|.+\|.+\|(.+)/>\1 \2/g' $RVFA; ## reformat RVDB fasta deflines to NCBI format
##diamond makedb --in $RVFA --taxonmap $A2T --taxonnodes $NODE --taxonnames $NAME --db $RVDB;
##cd $pwd

##### create diamond nr protein database WITH taxonomy info
### takes 1-2 hours and 150-200 gb - only make if required - ENSURE all fasta and taxid inputs are same update/version
##cd $DB; diamond makedb --in $NRFA --taxonmap $A2T --taxonnodes $NODE --taxonnames $NAME --db $NR;
##cd $pwd

################ STEP 001 ##############

### create and move into rungroup subdirectory
mkdir $pwd"/finalresult/IVEX001/"IVEX001_gms"$1"; cd $pwd"/finalresult/IVEX001/"IVEX001_gms"$1"

##### prepare key for ncvog_gi to ncvog_id ncvog_name, ncvog_function, ncvog_group, ncvog_abbrev, ncvog_legend replacement
awk -F"\t" 'BEGIN{OFS=";"};{print $2"\t"$5,$4,$6,$7,$8,$9}' $NCVKEY | sed -E 's/ncvog_gi=([0-9]+)\t(.+)/s\/_____\1_____\/\2\/g/g' > IVEX001__001_000;

##### diamond ALL proteins vs virus-only databases
### vs NCVOG
for f in $FA; do
	diamond blastp -k 10 -e 1e-3 --unal 1 -d $NCVOG -q $f --more-sensitive -o IVEX001_$(basename $f .faa).NCVOG.blast -f 6 qseqid salltitles bitscore;
done

### vs RVDB
for f in $FA; do
	diamond blastp -k 10 -e 1e-3 --unal 1 -d $RVDB -q $f --more-sensitive -o IVEX001_$(basename $f .faa).RVDB.blast -f 6 qseqid salltitles pident length evalue bitscore staxids;
done

########## keep only the best hit per gene, prepare NCVOG, name and function for each gene, remove bad characters from hit names
##### prepare RVDB blast
sort -t $'\t' -k 1,1V -k 6,6nr IVEX001_*.RVDB.blast | awk -F"\t" '!a[$1]++' | awk -F"\t" 'BEGIN{OFS="\t"};{if($2 =="\*") print $1,"none","none","none","none","none","none"; else print $0}' > IVEX001__001_002; ## reduce to best hit only and change no-hit fields to none
awk -F"\t" 'BEGIN{OFS="\t"};{if($2 !="none") gsub(/[^a-zA-Z0-9]/,"_",$2)}1' IVEX001__001_002 | sed 's/\t/;/g' | sed -Ee 's/_+/_/g' -Ee 's/^_|_$//g' -Ee 's/;_|_;/;/g' | awk -F";" 'BEGIN{OFS=";"};{print "qseqid="$1,"salltitles="$2,"pident="$3,"length="$4,"evalue="$5,"bitscore="$6,"staxids="$7}' | sed 's/=$/=none/g' > IVEX001__001_003;

###### prepare NCVOG blast
### labels genes with core NCLDV gene homology
sort -t $'\t' -k 1,1V -k 3,3nr IVEX001_*.NCVOG.blast | awk -F"\t" '!a[$1]++' | cut -f 2 | sed -E 's/^gi\|([0-9]+).*/_____\1_____/g' > IVEX001__001_004; ## get NCVOG gi
sed -f IVEX001__001_000 IVEX001__001_004 > tmp && mv tmp IVEX001__001_004; ## replace NCVOG with NCVOG name, function, etc
sed -i -E 's/^_____[0-9]+_____$|^\*$/ncvog_id=none;ncvog_name=none;ncvog_function=none;ncvog_group=none;ncvog_abbrev=none;ncvog_legend=none/g' IVEX001__001_004 ## assign empty none fields to non-NCVOG genes

################ STEP 002 ##############
#### NR blast of only genes with RVDB hits

##### extract only proteins with RVDB hits. To be blasted vs nr
awk -F";" '{if($2 !~"salltitles=none") print $1}' IVEX001__001_003 | sed 's/qseqid=//g' > IVEX001__002_000;
for f in $FA; do
	seqkit grep -n -f IVEX001__002_000 $f > IVEX001_$(basename $f .faa)_002_000.faa;
done;

##### diamond ALL proteins with RVDB hits vs nr, excluding viruses
for f in IVEX001_*_002_000.faa; do
	diamond blastp -k 15 -e 1e-5 --unal 1 --taxon-exclude 10239 -d $NR -q $f -o $(basename $f .faa).NR.blast -f 6 qseqid salltitles evalue bitscore staxids;
done

######## prepare NR blast
#### get rid of virus hits that escaped NR exclusion, by removing hits in which subject has hit vs RVDB
## Removing NR hits which have subjects with RVDB hits reduces false cellular hits which are actually integrated viral genes mislabelled as cellular in public databases
## However, the trade off is that genes which some HGT relationship between cells and virus may be given a cell score of 0, leading to a line of values at graph bottom; even for core genes, some cellular homologs can have higher scores than a viral homolog (e.g. NCLDV MCP)
## get fasta of NR subject hit proteins
cut -f 2 IVEX001_*.NR.blast | sed '/^\*$/d' | sed 's/ /\t/1' | cut -f 1 | sort -Vu > tmp
seqkit grep -f tmp $NRFA > tmp.faa
## blast NR subjects vs RVDB
diamond blastp -k 5 -e 1e-40 -d $RVDB -q tmp.faa -o IVEX001__002.hiddenvirus.blast -f 6 qseqid salltitles;

## get all hits in NR blast which have subjects with RVDB hits, mark with VIRHOM flag in salltitles
cut -f 1 IVEX001__002.hiddenvirus.blast | sort -Vu | grep -hwF -f - IVEX001_*.NR.blast | awk -F "\t" '{print $1"\t"$2"_VIRHOM\t"$3"\t"$4"\t"$5}' > tmp1
## get all non-duplicated (and therefore have NO RVDB hits) NR hits, keep unchanged
awk -F "\t" '{print $1"\t"$2"_VIRHOM\t"$3"\t"$4"\t"$5}' IVEX001_*.NR.blast | cat tmp1 - | sort -V | uniq -u | sed 's/_VIRHOM//g' > tmp2

## merge all VIRHOM and other hits. IMPORTANT - here we divide the VIRHOM hits by -1, to reduce (make negative for sorting step) their cellular score as it is very likely that these are hidden viral genes - it is arbitrary - ensure it is noted in methods/figures
awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4*-1"\t"$5}' tmp1 | cat - tmp2 > tmp_blast

## reduce to best hit only (by bitscore) and change no-hit fields to none (for proteins with no NR hits)
## if a negative viral bitscore is still the highest and therefore retained, it is treated as no hit, as it cannot be considered a reliable cellular hit
sort -t $'\t' -k 1,1V -k 4,4nr tmp_blast | awk -F"\t" '!a[$1]++' | awk -F"\t" 'BEGIN{OFS="\t"};{if($2 =="\*" || $2 ~ "_VIRHOM") print $1,"none","none","none","none"; else print $0}' > IVEX001__002_002;
awk -F"\t" 'BEGIN{OFS="\t"};{if($2 !="none") gsub(/[^a-zA-Z0-9]/,"_",$2)}1' IVEX001__002_002 | sed 's/\t/;/g' | sed -Ee 's/_+/_/g' -Ee 's/^_|_$//g' -Ee 's/;_|_;/;/g' | awk -F";" 'BEGIN{OFS=";"};{print "qseqid="$1,"salltitles="$2,"evalue="$3,"bitscore="$4,"staxids="$5}' | sed 's/=$/=none/g' > IVEX001__002_003;

## add back genome proteins with no RVDB hits
grep ">" $FA | sed -E 's/.+>(.+)/\1/g' > IVEX001__002_004;
grep ">" IVEX001_*_002_000.faa | sed -E 's/.+>(.+)/\1/g' > IVEX001__002_005;
sort -V IVEX001__002_004 IVEX001__002_005 | uniq -u | awk -F";" 'BEGIN{OFS=";"};{print "qseqid="$1,"salltitles=none","evalue=none","bitscore=none","staxids=none"}' >> IVEX001__002_003;
sort -t";" -k 1,1V IVEX001__002_003 -o IVEX001__002_003;

################ STEP 003 ##############
####### add in full taxonomy lineages #####
#### for viruses in RVDB hits
cut -d";" -f 7 IVEX001__001_003 | sed -Ee 's/^staxids=([0-9]+)/_____\1_____/g' -e 's/staxids=none//g' > IVEX001__003_001; ## prepare taxid field

#### reformat taxid key file: remove all characters except letters, digits, underscores, and tab. Reduces mutliple underscores in series to 1, removes underscores trailing fields and lines
sed -e 's/;/\t/g' -e 's/[^Aa-Zz0-9\t]/_/g' -Ee 's/_+/_/g' -Ee 's/^_|_$//g' -Ee 's/\t_|_\t/\t/g' $VT | sed -e 's/\t/;/g' -e 's/;/\t/1' -Ee 's/^([0-9]+)\t/_____\1_____\t/g' > IVEX001__003_002;
sort -Vu IVEX001__003_001 | grep -wF -f - IVEX001__003_002 > tmp && mv tmp IVEX001__003_002; ## reduce key to relevant lines
sed -i -E 's/(_____[0-9]+_____)\t(.+)/s\/\1\/\2\/g/g' IVEX001__003_002 ## reformat for sed -f
sed -f IVEX001__003_002 IVEX001__003_001 > IVEX001__003_003 ## replace taxids with full virus lineage
sed -i -e 's/^/vlineage=/g' -Ee 's/^vlineage=_____[0-9]+_____$|^vlineage=$/vlineage=none/g' -e 's/;/-/g' IVEX001__003_003 ## fix empty lines

################ STEP 004 ##############
##### self hit blast - query and db is proteins with RVDB hits
## make database for each genome
rm -f IVEX001_*_004_self.dmnd
for f in IVEX001_*_002_000.faa; do
	diamond makedb --in $f --db $(basename $f _002_000.faa)_004_self;
done

#### run self blasts
for f in IVEX001_*_002_000.faa; do
	diamond blastp -k 200 -e 20 --masking 0 --unal 1 -d $(basename $f _002_000.faa)_004_self.dmnd -q $f -o $(basename $f _002_000.faa).self.blast -f 6 qseqid sseqid bitscore;
done

####### prepare self hits
## get all the successful self-hits
awk -F"\t" '{if($1 == $2) print $0; else print $1"\tnone\tnone"}' IVEX001_*.self.blast | sort -t $'\t' -k 1,1V -k 3,3nr | awk -F"\t" '!a[$1]++' > IVEX001__004_003
## list proteins that failed to self hit
awk -F "\t" '{if($2 =="none" && $3 =="none") print $1}' IVEX001__004_003 > tmp;
## for proteins that failed to self hit, get highest scoring hit instead, add back to other hits
cut -f 1,3 IVEX001_*.self.blast | grep -wF -f tmp - | awk -F"\t" 'BEGIN{OFS="\t"};{print $1,$1,$2}' | sort -t $'\t' -k 1,1V -k 3,3nr | awk -F"\t" '!a[$1]++' >> IVEX001__004_003
## remove failed hit null line or empty line if all self hits succeeded
sed -i -e '/\tnone\tnone/d' -e '/^$/d' IVEX001__004_003
## add field labels
awk -F"\t" 'BEGIN{OFS=";"};{print "qseqid="$1,"sseqid="$2,"bitscore="$3}' IVEX001__004_003 > tmp && mv tmp IVEX001__004_003
## for proteins that failed to match any self-protein, change to none
sed -i 's/bitscore=-1/bitscore=none/g' IVEX001__004_003

## add back proteins with no RVDB hits
sort -V IVEX001__002_004 IVEX001__002_005 | uniq -u | awk -F";" 'BEGIN{OFS=";"};{print "qseqid="$1,"sseqid=none","bitscore=none"}' >> IVEX001__004_003;
sort -t";" -k 1,1V IVEX001__004_003 -o IVEX001__004_003;

####### calculate relative bitscores for RVDB and NR hits per gene (hit bitscore divided by self hit score)
## get bitscores RVDB, NR, self
paste -d ";" IVEX001__001_003 IVEX001__002_003 IVEX001__004_003 | cut -d";" -f 6,11,15 | sed -e 's/bitscore=//g' -e 's/none/0.0/g' | awk -F";" '{if($3 =="0.0") print "NA;NA;NA"; else print $0}' > IVEX001__004_004;
## calc score, any score above 1, set to 1
awk -F";" '{if($3 !~ "NA") print $1/$3; else print $0}' IVEX001__004_004 | awk -v n=1.0 '{if($1 > n) print "1.0"; else print $1}' | paste IVEX001__004_004 - | awk -F"\t" '{if($1 ~ "NA;") print "rbitscore=NA"; else print "rbitscore="$2}' > IVEX001__004_005; ## RVDB rbitscore
awk -F";" '{if($3 !~ "NA") print $2/$3; else print $0}' IVEX001__004_004 | awk -v n=1.0 '{if($1 > n) print "1.0"; else print $1}' | paste IVEX001__004_004 - | awk -F"\t" '{if($1 ~ "NA;") print "rbitscore=NA"; else print "rbitscore="$2}' > IVEX001__004_006; ## NR rbitscore

################ STEP 005 ##############
#### FINAL GFF ASSEMBLY
#### merge all outputs columns: viral RVDB blast, viral rbitscore, viral lineage, NCVOG info, cellular NR blast, cellular rbitscore
paste -d";" IVEX001__001_003 IVEX001__004_005 IVEX001__003_003 IVEX001__001_004 IVEX001__002_003 IVEX001__004_006 | cut -d";" --complement -f 1,16 | sed -e 's/^/v/1' -e 's/;/;v/1' -e 's/;/;v/2' -e 's/;/;v/3' -e 's/;/;v/4' -e 's/;/;v/5' -e 's/;/;v/6' > IVEX001__005_001
### list all proteins in final file (blasted vs RVDB)
cut -d ";" -f 1 IVEX001__001_003 | sed -e 's/qseqid=//g' > IVEX001__005_002
## extract from gff only mRNA lines and protein name present in final file
## ; is changed to tab in of additional fields, then the 1-9 (10 to initially include filename) tab delimited field of gff are pasted to final results
awk -F"\t" '{if($3 =="mRNA") print FILENAME"\t"$0}' $GF | grep -wF -f IVEX001__005_002 - | sed 's/;/\t/1' | sort -t $'\t' -k 10,10V | cut -f 1-10 | paste -d ";" - IVEX001__005_001 > tmp && mv tmp IVEX001__005_001
## additional fields beyond 9 are added to final file
awk -F"\t" '{if($3 =="mRNA") print $0}' $GF | grep -wF -f IVEX001__005_002 - | sed 's/;/\t/1' | sort -t $'\t' -k 9,9V | cut --complement -f 1-9 | paste -d ";" IVEX001__005_001 - | sed 's/;$//g' > tmp && mv tmp IVEX001__005_001
## get mRNA lines from gff that do not have proteins in final file (not blasted vs RVDB)
awk -F"\t" '{if($3 =="mRNA") print FILENAME"\t"$0}' $GF | grep -vwF -f IVEX001__005_002 - >> IVEX001__005_001
## get all non-mRNA lines
awk -F"\t" '{if($3 !="mRNA") print FILENAME"\t"$0}' $GF >> IVEX001__005_001
## fix genome name field, sort gff, and split by genome, then finally remove genome field 1
sed -E 's/^\/.+\/gms[0-9]+_(.+).gff\t/\1\t/g' IVEX001__005_001 | sort -t $'\t' -k 1,2V -k 5,6V > tmp && mv tmp IVEX001__005_001
awk -F"\t" '{print > "IVEX001_final_"$1".gff"}' IVEX001__005_001
for f in IVEX001_final_*.gff; do cut --complement -f 1 $f > tmp && mv tmp $f; done
rm -f IVEX001__005_001

###### make final check file which shows output line/gene counts and indicate where any errors may occur
### an important on proteins with neither viral nor cellular hits; they are only VIRAL ORFans, because proteins without an initial viral hit were NOT blasted against NR - these proteins are not true ORFans because they have hits to the NR (cellular) database
grep ">" $FA | wc -l | sed 's/^/## counts in each section must match or something has gone wrong\n## No. proteins in IVEX000 fasta, $FA\n/g' > IVEX001__005_final_check;
cut -f 1 IVEX001_*.RVDB.blast | sort -Vu | wc -l | sed 's/^/## No. proteins from $FA successfully blasted vs RVDB\n/g' >> IVEX001__005_final_check;
wc -l IVEX001__001_003 IVEX001__004_005 IVEX001__003_003 IVEX001__001_004 IVEX001__002_003 IVEX001__004_006 | sed -z 's/^/## No. lines in all files of final merge\n/g' >> IVEX001__005_final_check;
sed -i -z 's/$/\n\n/g' IVEX001__005_final_check
wc -l IVEX001__002_000 | sed 's/^/## (A) No. proteins with hits vs RVDB\n/g' >> IVEX001__005_final_check;
cut -f 1 IVEX001_*.NR.blast | sort -Vu | wc -l | sed 's/^/## (B) No. proteins with hits vs RVDB successfully blasted vs NR\n/g' >> IVEX001__005_final_check;
awk -F"\t" '{if($1 == $2) print $1}' IVEX001_*.self.blast | sort -Vu | wc -l | sed 's/^/## (C) No. proteins with hits vs RVDB successfully self blasted; can mismatch - self-hit failures were corrected if C+D+E = A or B\n/g' >> IVEX001__005_final_check;
awk -F";" '{if($2 !~ "sseqid=none" && $3 ~ "bitscore=none") print $0}' IVEX001__004_003 | wc -l | sed 's/^/## (D) No. proteins which failed to hit ANY other self-proteins\n/g' >> IVEX001__005_final_check;
awk -F";" '{if($2 !~ "sseqid=none" && $3 ~ "bitscore=none") print $0}' IVEX001__004_003 | wc -l > tmp;
awk -F"\t" '{if($1 == $2) print $0; else print $1"\tnone\tnone"}' IVEX001_*.self.blast | sort -t $'\t' -k 1,1V -k 3,3nr | awk -F"\t" '!a[$1]++' | awk '/\tnone\tnone/' | wc -l | paste - tmp | awk -F "\t" '{print $1-$2}' | sed 's/^/## (E) No. proteins that failed to self-hit, but highest self-protein hit was used instead\n/g' >> IVEX001__005_final_check;
awk '/rbitscore=[0-9]/' IVEX001__004_005 | wc -l | sed 's/^/## No. proteins with hits vs RVDB with rbitscore + (D) to correct self-hit failure\n/g' >> IVEX001__005_final_check;
awk '/rbitscore=[0-9]/' IVEX001__004_006 | wc -l | sed 's/^/## No. proteins with hits vs NR with rbitscore + (D) to correct self-hit failure\n/g' >> IVEX001__005_final_check;
cat IVEX001_final_*.gff | awk -F"\t|;" '{if($3 =="mRNA" && $10 ~ "vsalltitles=" && $10 !~ "vsalltitles=none" ) print $1}' - | wc -l | sed 's/^/## No. proteins with hits vs RVDB in final gff\n/g' >> IVEX001__005_final_check;
sed -i -z 's/$/\n\n/g' IVEX001__005_final_check
cat IVEX001_final_*.gff | wc -l | sed 's/^/## No. lines final gff\n/g' >> IVEX001__005_final_check;
cat $GF | wc -l | sed 's/^/## No. lines in IVEX000 gff, $GF\n/g' >> IVEX001__005_final_check;
rm -f tmp*

mv IVEX001_final_*.gff ..
