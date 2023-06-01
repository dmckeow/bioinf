#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition long
#SBATCH --cpus-per-task 12
#SBATCH --mem 128GB

####################### INPUTS #####################################
###### load software ######
module unload blast
module load seqkit/0.14.0; module load blast/2.9.0; module load diamond/2.0.9; module load mafft/7.407; module load fasttree/2.1.10; module load bedtools/2.30.0
GeneMarkS="/projet/fr2424/sib/dmckeown/genemark_suite_linux_64/gmsuite/gmsn.pl"

##### command line arguments required #####
## $1 = viral protein fasta to be query

##### input files for analysis #####
PHEX_FA="/shared/projects/phaeoexplorer_virus/phaeoex_screen/input/IVEX000_genomes/*.fa.gz"
STRA_FA="/scratch2/fr2424/sib/dmckeown/db/stramenopiles/fasta/*.fna"

#VFA="/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/IVEX_stramenopiles/IVEX_stramenopiles_vproteins.faa"
#STRA_KEY="/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/IVEX_reports/paper_all_genomes/NCBI_stramenopile_genomes"

###### input file for databases ######
STRADB="/scratch2/fr2424/sib/dmckeown/db/stramenopiles/fasta/stra"
NR="/db/nt/current/blast/nt"
##### NCVOG defintion files #####
NCVKEY="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/IVEX_NCVOG_006"


### ON 30.04.2022 did : mv IVEX_stramenopiles_rg* /scratch2/fr2424/sib/dmckeown/IVEX_stramenopiles

####################### SCRIPT #####################################
############################################################# V2 #####################################################################
########### RUN FIRST BEFORE RUNNING SCRIPT:

#### identify all genomes to identified, and get a list of all nucleotide fasta accessions for each genome - this will be used to restrict the BLAST search to only these subjects
#### fix your fasta deflines so they are consistent (keep record of the original deflines in case you need them later)

#### create a fasta file of all the viral proteins to be used as blast queries
### ensure all relevant info, especially protein type/family can be easily identified in the deflines e.g.:
## >0023_D5|gi|109287999|ref|YP_654693.1|hypothetical_protein_MIV121R|Invertebrate_iridescent_virus_3
#### ensure the fasta deflines have problem characters removed e.g.:
##for f in vlist_*.faa; do sed -i -e 's/[^A-Za-z0-9>|[\.]/_/g' -Ee 's/_+/_/g' -e 's/_$//g' -Ee 's/\[/|/g' -e 's/_\||\|_/|/g' $f; done


################ STEP 000 ##############
### set name for outputs
n1=$(echo "$1" | sed 's/.*\///g' | cut -d "." -f 1)

### make database of local files (phaeoexplorer)
#for f in $PHEX_FA; do g=$(basename "$f" .fa.gz); zcat $f | awk -v h=$g '{gsub(">",">"h"|")}1'; done > tmp.fa
#for f in $STRA_FA; do g=$(basename "$f" .fna); cat $f | awk -v h=$g '{gsub(">",">"h"|")}1'; done >> tmp.fa

#makeblastdb -in tmp.fa -title "Phaeoexplorer and NCBI non-phaeophyceae Stramenopiles" -dbtype nucl -out phexstra

### run blast (tblastn uses protein query vs translated nucleotide subject) vs local files database
#tblastn -db "/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/IVEX_stramenopiles/phexstra" -query $1 -out "$n1".PHEXSTRA.blast -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore salltitles" -max_target_seqs 2000 -evalue 10e-5

rm -f tmp*

################ STEP 001 ##############
### process tBLASTn results
### flip coordinates of reverse matches and mark all as forward or reverse
#awk -F"\t" 'BEGIN{OFS="\t"};{if($7 > $8) print $1,$2,$3,$4,$5,$6,$8,$7,$9,$10,$11,"-"; else print $0,"+"}' "$n1".PHEXSTRA.blast > tmp_IVEX_stramenopiles1

### remove redundant matches, where they align to exactly same region and keeping the one with the lowest evalue and highest bitscore
#sort -t $'\t' -k 2,2V -k 7,7n -k 8,8n -k 9,9n -k 10,10nr tmp_IVEX_stramenopiles1 | awk -F"\t" '!a[$2,$7,$8]++' > tmp_IVEX_stramenopiles2

#cp tmp_IVEX_stramenopiles2 tmp_IVEX_stramenopiles2.5

### identify subjects that overlap with previous one (in order of start and end, same contig/genome)
#be=0
#af=1
#until [[ "$be" -eq "$af" ]]
#  do
#      be=`wc -l tmp_IVEX_stramenopiles2 | sed 's/ .*//g'`
#      sort -t $'\t' -k 2,2V -k 7,7n -k 8,8n tmp_IVEX_stramenopiles2 | awk -F"\t" 'BEGIN{OFS="\t"};{if($2 == A && $7 < b) print $0,"OVERLAP"; else print $0} {A=$2} {b=$8}' | sed '/\tOVERLAP/d' > tmp && mv tmp tmp_IVEX_stramenopiles2
#      af=`wc -l tmp_IVEX_stramenopiles2 | sed 's/ .*//g'`
#    sleep 1
#done

### merge all blast results back together, in order of start and end, now with blocks of overlapping frames marked
#cat tmp_IVEX_stramenopiles2.5 tmp_IVEX_stramenopiles2 | sort -V | uniq -u | sed 's/$/\tOVERLAP/g' | cat - tmp_IVEX_stramenopiles2 | sort -t $'\t' -k 2,2V -k 7,7n -k 8,8n > tmp_IVEX_stramenopiles3

### tag the first frame of each overlapping block, each start marked by OVERLAP_START (same done for frames with no overlaps)
#awk -F"\t" '{if($13=="") print $0"\tOVERLAP_START"; else print $0}' tmp_IVEX_stramenopiles3 > tmp_IVEX_stramenopiles4

### split by overlapping blocks of frames subjects
#csplit -z -s -f tmp_IVEX_stramenopiles4_ tmp_IVEX_stramenopiles4 /OVERLAP_START/ {*}

### keep the frame with the lowest evalue and highest bitscore per overlap block
#for f in tmp_IVEX_stramenopiles4_*; do sort -t $'\t' -k 9,9n -k 10,10nr $f | cut --complement -f 13 | awk -F"\t" '!a[$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11]++' > 2_"$f"; done
#cat 2_tmp_IVEX_stramenopiles4_* > IVEX_stramenopiles.blast.reduced

#for f in tmp_*; do rm -f $f; done; for f in *_tmp_*; do rm -f $f; done;

### make quick summary of reduced blast results - initial look at data
#sed -E 's/[0-9]+\.[0-9]+e-/1.0e-/g' IVEX_stramenopiles.blast.reduced | awk -F"\t" '{if($9 <= 1.0e-40) print $0"\tlow"; else print $0}' | awk -F"\t" '{if($9 <= 1.0e-20 && $0 !~ /\tlow$/) print $0"\tmed"; else print $0}' | awk -F"\t" '{if($9 <= 1.0e-5 && $0 !~ /\tlow$|\tmed$/) print $0"\thigh"; else print $0}' | cut -f 1,2,13 | sed 's/|/\t/g' | cut -f 1,8,10 | sort -V | uniq -c | sed -Ee 's/ +/\t/g' -e 's/^\t//g' | sed -z 's/^/count\tcoregene\tgenome\tevalue\n/1' > IVEX_stramenopiles.blast.reduced.quicksummary
#sed -i -e 's/\thigh\t/\t<=e-05\t/g' -e 's/\tmed\t/\t<=e-20\t/g' -e 's/\tlow\t/\t<=e-40\t/g' IVEX_stramenopiles.blast.reduced.quicksummary

######### MANUAL STEP 001 ###############
##### add taxa fields here
## AND add categories with zero values
## count   coregene        genome  evalue  fulltaxa        taxa1   taxa2   taxa3
## 0       0022_MCP        Ascophyllum-nodosum_MALE        <=e-05  cellular_organisms;Eukaryota;Sar;Stramenopiles;Ochrophyta;PX_clade;Phaeophyceae;order;family;genus;species                                    Ochrophyta       PX_clade        Phaeophyceae

###### summarise as percentages of genomes per taxa1_taxa3 groups
#sed -Ee 's/^0\t/Z_Z_Z\t/g' -Ee 's/^[0-9]+\t/1\t/g' -Ee 's/Z_Z_Z\t/0\t/g' IVEX_stramenopiles.blast.reduced.quicksummary > IVEX_stramenopiles.blast.reduced.quicksummary2 ## make presence or absence instead of gene per geomne count
#cut -f 3,6,8 IVEX_stramenopiles.blast.reduced.quicksummary2 | sed '1,1d' | sort -Vu | cut -f 2,3 | sort -V | uniq -c | sed -Ee 's/ +/\t/g' -e 's/^\t//g' | awk -F"\t" '{print "s/"$2"-"$3"/"$1"/g"}' > fixer
#cut -f 6,8 IVEX_stramenopiles.blast.reduced.quicksummary2 | awk -F"\t" '{print $1"-"$2}' > tofix
#sed -f fixer tofix > fixed; sed -i 's/taxa1-taxa3/grouptotal_taxa1_taxa3/g' fixed
#paste IVEX_stramenopiles.blast.reduced.quicksummary2 fixed | cut --complement -f 3,5 > IVEX_stramenopiles.blast.reduced.quicksummary3
#awk -F"\t" '{if($1 !~ "count") print $0"\t"($1/$7)*100; else print $0"\tpercent_of_grouptotal_taxa1_taxa3"}' IVEX_stramenopiles.blast.reduced.quicksummary3 > tmp && mv tmp IVEX_stramenopiles.blast.reduced.quicksummary3

#### aggregate into percent per group
#awk -F"\t" 'BEGIN{OFS="_"};{print $0"\t"$2,$3,$4,$5,$6}' IVEX_stramenopiles.blast.reduced.quicksummary3 | sed 's/<=//2' | awk -F"\t" '{print $0 > "tmp_agg_"$9}'
#for f in tmp_agg_*; do awk -F "\t" '{sum+=$8;} END{print $0"\t"sum}' $f | cut --complement -f 1,8,9 > tmp && mv tmp $f ; done
#cat tmp_agg_*  | sed '/^coregene\t/d' | sed -z 's/^/coregene\tevalue\ttaxa1\ttaxa2\ttaxa3\tgrouptotal_taxa1_taxa3\tpercent_of_grouptotal_taxa1_taxa3\n/1' > IVEX_stramenopiles.blast.reduced.quicksummary4

############################################################# MANUAL STEP 002 #####################################################################
###### R figures can be made BEFORE doing this, but if this step is done, generate the figures again with the new files from this step

##### align each protein family with viral proteins, remove everything that doesnt align or is too short
##### before aligning, examine the .reduced file and decide which proteins are good matches and remove all proteins that are too short/bad scoring - as save as .reduced.curated
#### then run quick alignment and phylogeny and inspect tree/alignment to see which proteins to include (select in Dendroscope)
#### add public proteins!


### make local fasta of all input genomes (temporary)
for f in $PHEX_FA; do g=$(basename "$f" .fa.gz); zcat $f | awk -v h=$g '{gsub(">",">"h"|")}1'; done > tmp.fa
for f in $STRA_FA; do g=$(basename "$f" .fna); cat $f | awk -v h=$g '{gsub(">",">"h"|")}1'; done >> tmp.fa
sed -i -E 's/(>PUBLIC_Undaria-pinnatifida\|LG)([0-9])$/\10\2/g' tmp.fa
cut -f 2,7,8 IVEX_stramenopiles.blast.reduced | awk -F"\t" '{print $1"\t"$2-1"\t"$3-1}' > IVEX_stramenopiles_hits.bed

bedtools getfasta -fi tmp.fa -bed IVEX_stramenopiles_hits.bed -fo IVEX_stramenopiles_hits.fa

#for f in IVEX_stramenopiles_blast_*.faa; do
#ls $f | sed -E 's/IVEX_stramenopiles_blast_(.+)\.faa/\1/g' | grep -f - $VFA | sed 's/>//g' | seqkit grep -n -f - $VFA | cat - $f | mafft --leavegappyregion --reorder --auto - > $(basename $f .faa).aln
#done

#for f in IVEX_stramenopiles_blast_*.aln; do FastTree $f > $(basename $f .aln).fasttree; done
