#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition long
#SBATCH --cpus-per-task 16
#SBATCH --mem 32GB

####################### INPUTS #####################################
###### load software ######
module load seqkit/0.14.0; module load mafft/7.407; module load fasttree/2.1.10; module load raxml-ng/0.9.0; module load raxml/8.2.12;

######### IMPORTANT
### this script is mostly an example workflow to get proteins from public and private data and prepare them for phylogeny
### the exception is the commands in STEP002 - section 3, which run the alignment and tree generation

##### command line arguments required #####
## $1 = [GENE].uniprot.phex.fa
## $2 = number of bootstraps

##### input files for analysis #####
pwd=$(pwd) ## get path to project folder
FA=$pwd"/input/IVEX000/gms"$1"_*.faa" ## the protein fastas prepared by IVEX000.sh
GF=$pwd"/input/IVEX000/gms"$1"_*.gff" ## the annotations prepared by IVEX000.sh

###### input file for databases ######

##### NCVOG defintion files #####
NCVKEY=$pwd"/input/IVEX000_NCVOG" ## made by script IVEX000_NCVOG.sh

########################## STEP 001 #################################
############# get representative members of protein family from public databases

####################################################################################
########## 1. get proteins from public database
## 1a. go to https://www.ebi.ac.uk/interpro/search/sequence/ and search protein sequence from desired protein family

## 1b. select a protein domain that represents the protein family and click "View [pfam accession] in Pfam"

## 1c. select the "species" tab, then select desired taxa on "sunburst" or "tree"

## 1d. download accessions or fasta of selected taxa to subdirectory named after gene of interest
## to get accession list from fasta:
##### grep ">" selected_sequences.fa | sed -e 's/\//\t/1' -Ee 's/>(.*)\t.*/\1/g' > selected_sequences.list

####################################################################################
########## 2. get better version of the public fasta with full species names
## 2a. go to https://www.uniprot.org/uploadlists/ and upload your selected_sequences.list and submit

## 2b. select all entries and download FASTA(canonical)

## 2c. change fasta deflines to species name[TAB]accession[TAB]gi AND rename
##### sed -Ee 's/>.+\|(.+)\|.+(OS=.+) (OX=[0-9]+) .+/>\2\t\1\t\3/g' -Ee 's/OS=|OX=//g' uniprot-yourlist_[etc].fasta > selected_sequences.uniprot.fa

## 2d. reformat fasta deflines to species name[dash]accession[dash]gi and remove all problematic characters (anything except letters, digits, dash, & TAB)
##### awk '{gsub(/[^a-zA-Z0-9\t]/,"_")}1' selected_sequences.uniprot.fa | sed -Ee 's/_+/_/g' -e 's/^_/>/g' -Ee 's/_$//g' -Ee 's/\t_|_\t/\t/g' -e 's/\t/-/g' > tmp && mv tmp selected_sequences.uniprot.fa

####################################################################################
########## 3. get any missing taxa or important members
## 3a. prepare taxid key file $VT
##### sed -E 's/^([0-9]+)\t/_____\1_____\t/g' /projet/fr2424/sib/dmckeown/db/taxonkit/taxonkit_viruses_full_lineage > tmp

## 3b. prepare list of all species with full taxonomy lineage present in your fasta - use this to manually check who is missing
##### grep ">" selected_sequences.uniprot.fa | cut -d "-" -f 3 | sort -Vu | sed -E 's/([0-9]+)/_____\1_____/g' | grep -wF -f - tmp | cut -f 2 | sort -V > selected_sequences.uniprot.specieslist

## 3c. you may want to get all taxonomy lineages from your [TAXA] of interest to compare your list to
##### grep -wF "[TAXA]" tmp | cut -f 2 | sort -V | cut -d ";" --complement -f [HIGHER_TAXA_FIELDS]

## 3d. using list files or other information, fetch missing proteins and reformat to match other fastas (keep in separate file until merge reccommended)
## ideally find them through: https://www.uniprot.org/uploadlists/ and reformat as with steps 2c and 2d:
##### sed -i -Ee 's/>.+\|(.+)\|.+(OS=.+) (OX=[0-9]+) .+/>\2\t\1\t\3/g' -Ee 's/OS=|OX=//g' selected_sequences.uniprot.missing.fa
##### awk '{gsub(/[^a-zA-Z0-9\t]/,"_")}1' selected_sequences.uniprot.missing.fa | sed 's/^_/>/g' | sed -Ee 's/_+/_/g' -Ee 's/_$//g' -Ee 's/\t_|_\t/\t/g' -e 's/\t/-/g' > tmp && mv tmp selected_sequences.uniprot.missing.fa

## 3e. merge auto and missing sequences and sort
##### cat selected_sequences.uniprot.fa selected_sequences.uniprot.missing.fa | seqkit sort -N - > selected_sequences.uniprot.complete.fa

####################################################################################
########## 4. remove redundant proteins (FROM SAME GENOME) such as paralogs, etc, or badly aligned proteins
#### CHECK FOR IF GENES HAVE INTRONS AND USE CORRECT TRANSLATION
## 4a. to check for multiple proteins from same genome - if no result, then you can probably skip to 4d.
##### grep ">" selected_sequences.uniprot.complete.fa | sed 's/>//g' > to_keep
## list duplicated proteins within same genome
##### cut -d "-" -f 1 to_keep | sort -V | uniq -d

## 4b. do quick alignment and tree - manually examine tree , 4a output and remove redundant proteins from to_keep
##### mafft --leavegappyregion --reorder --auto selected_sequences.uniprot.complete.fa > selected_sequences.uniprot.complete.aln
##### FastTree selected_sequences.uniprot.complete.aln > selected_sequences.uniprot.complete.fasttree

## 4c. remove redundant seqs to make final fasta
##### seqkit grep -n -f to_keep selected_sequences.uniprot.complete.fa > selected_sequences.uniprot.final.fa

## 4d. save accessions and gis for supplementary info and remove from fasta
##### grep ">" selected_sequences.uniprot.final.fa | sed -e 's/>//g' -e 's/-/\t/g' | sed -z 's/^/species\taccession\tgi\n/1' > selected_sequences.uniprot.final.accessions
##### sed -i -E 's/(>.+)-.+-.+/\1/g' selected_sequences.uniprot.final.fa

## 4e. rename final outputs
## selected_sequences.uniprot.final.accessions to [gene_name].uniprot.final.accessions
## selected_sequences.uniprot.final.fa to [gene_name].uniprot.final.fa

####################################################################################
############# all trees together
########## 5. final name fix (do this once steps 1-5 done for each protein family, for consistency and concatenation)

## 5a. count occurrance of each genome across all protein fastas - those with same no. of copies as no. proteins are good.
## The rest are those genomes missing a certain protein due to biology or error
## check that it is ok to remove certain genomes entirely so that the same species are in all proteins (for consistency purposes)
## do not remove species of interest that biologically lack a certain gene - this will be done for the concatenation step later
## check for variations of names that should be the same
##### grep ">" *uniprot.final.fa | sed 's/:>/\t/g' | cut -f 2 | sort -V | uniq -c | sed -E 's/ +([0-9]+) /\1\t/g' | sort -k 2,2V -k 1,1n > counts
## get all deflines below number of proteins families (n)
##### grep ">" *uniprot.final.fa | sed 's/:>/\t/g' | cut -f 2 | sort -V | uniq -c | sed -E 's/ +([0-9]+) /\1\t/g' | sort -k 2,2V -k 1,1n | awk -F"\t" '{if($1 < n) print $0}'
## fix the names in fastas (keep backup copies) so they all match in fastas
## add in any missing sequences to .fa and accessions
## make list of seq deflines to keep, removing those that lack some proteins due to error, but don't matter to the tree

## 5b. manually check the names of selected_sequences.uniprot.final.fa and fix anything too long or mispelled, etc as these names will be on your tree

## 5c. add fixed names to accessions legend file
##### grep ">" selected_sequences.uniprot.final.fa | sed -e 's/>//g' -ze 's/^/tree_name\n/1' | paste - selected_sequences.uniprot.final.accessions > tmp && mv tmp selected_sequences.uniprot.final.accessions

######################
############# 6. prep for concatenation
##### grep all deflines and keep only those that occur N times (N = number of proteins concatenated)
##### ensure each fasta is sorted the same

########################## STEP 002 #################################
############# get putative members of protein family from your data
########### 1. get initial gff and proteins
## 1a. get all gff lines for your protein (N) and prep for excel view

#####for f in ../IVEX002/*.gff.gz; do n="$(basename "$f" .gff.gz)"; zcat $f | awk -F "\t|;" -v n=$n 'BEGIN{OFS="\t"};{if($3 ~"mRNA" && $17 ~ "Nucleocytoviricota") print n,$0}' | sed -e 's/;/\t/g' -Ee 's/[Aa-Zz]+_[Aa-Zz]+=|[Aa-Zz]+=//g' | cut -f 1,2,10-17,19,24,25,29; done | sort -t $'\t' -k 6,6nr | sed -z 's/^/genome\tcontig\tID\tvsalltitles\tvpident\tvlength\tvevalue\tvbitscore\tvstaxids\tvrbitscore\tncvog_id\tncvog_legend\tsalltitles\trbitscore\n/g' > IVEX003_final_gff

## 1b. check for proteins with too short blast alignment lengths, or low bitscores
## mark these proteins in excel copy

## genome name added to protein names in case genomes share identical protein names
## remove them from original (split is useful for this) and awk -F ";" '{print $1"__"$3}' > to_keep_1

## 1c. get fasta of to_keep_1 - check in case of identical proteins names between genomes
## recommended to add genome name to these protein names that may be confused with another genome
##### touch to_keep_1.faa; for f in /projet/fr2424/sib/dmckeown/phaeoex_screen/input/GeneMarkS/gms*_*.faa; do srun awk '{gsub(">",FILENAME"__")}1' $f | sed -E 's/.+\/gms[0-9]+_(.+)\.faa/>\1/g' | seqkit grep -n -f to_keep_1 - >> to_keep_1.faa; done

########### 2. test alignment and tree of initial gff and proteins

## 2a. list duplicated proteins within same genome - you may wish to remove if redundancy is high
##### grep ">" to_keep_1.faa | sort -V

## 2b. do quick alignment and tree - manually examine alignment and tree
## just remove any proteins that are short enough to mean loss of info; partial proteins, etc
## keeping paralogs from same genome is ok if you lack any experimental data showing which is the "real" protein
##### mafft --leavegappyregion --reorder --auto to_keep_1.faa > to_keep_1.aln
##### FastTree to_keep_1.aln > to_keep_1.fasttree
## save good proteins names as to_keep_2 to fetch with seqkit

## 2c. remove bad seqs to make fasta2
##### seqkit grep -n -f to_keep_2 to_keep_1.faa > [GENE].phex.final.faa

## 2d. fix any proteins in this file (for example if truncated with another protein, etc)

## 2e. fix names (for display in tree) and record changes in
##### grep ">" [GENE].phex.final.faa | sed 's/>//g' > [GENE].phex.final.names
##### cp [GENE].phex.final.names [GENE].phex.final.new_names
## manually put new name in empty sed replace in new_names and merge with old names into name fix
##### paste -d "/" [GENE].phex.final.names [GENE].phex.final.new_names | sed -E 's/(.+)/s\/\1\/g/g' > name_fix
## then do
##### sed -f name_fix [GENE].phex.final.faa > [GENE].phex.FINAL.faa

## 2f. merge original and edited names key
## record key of name changes made
##### sed -Ee 's/s\/|\/g//g' -e 's/\//\t/g' name_fix > DNApolB_name_change_key

########### 3. add public proteins, run phylogeny
## 3a. add public proteins
##### somewhere do mkdir trees and move there
##### cat [GENE].phex.FINAL.faa ../../public/[GENE]/[GENE].uniprot.FINAL.fa > [GENE].phex.public.FINAL.faa
#cd $pwd"/finalresult/IVEX003"

## 3b. align (RUN IN SBATCH)
#mafft --leavegappyregion --reorder --auto $1 > $(basename $1 .faa).aln

## 3c. run tree

## construct ML phylogenetic trees (run with 1 cpus per task)
## adjust bootstraps 50 for ok tree, 200 for good tree
#FastTree $(basename $1 .faa).aln > $(basename $1 .faa).fasttree
#raxml-ng --all --msa $(basename $1 .faa).aln --model LG+G8+F --bs-trees $2 --prefix $(basename $1 .faa)_"$2"bs.raxml-ng; ## normal tree

#### old RaxML if -ng too slow
#raxmlHPC-PTHREADS -T 16 -d -f ae -m PROTGAMMAGTR -p $RANDOM -x $RANDOM -N $2 -s $(basename $1 .faa).aln -n $(basename "$1" .fa)_"$2"bs.raxml;
n1=$(echo "$1" | cut -d "." -f 1)
raxmlHPC-PTHREADS -T 16 -d -f ae -m GAMMAGTR -p $RANDOM -x $RANDOM -N $2 -s ${n1}.aln -n ${n1}_${2}bs.raxml;
