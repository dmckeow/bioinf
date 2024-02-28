#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition long
#SBATCH --cpus-per-task 1
#SBATCH --mem 32GB

###### load software ######
module load diamond/0.9.36; module load seqkit/0.14.0; module load gblocks/0.91b; module load mafft/7.407; module load raxml-ng/0.9.0; module load fasttree/2.1.10; module load raxml/8.2.12;

###### input file variables ######
P1="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/00__PUBLIC_GENOMES_PROTEOMES/*.fa"
P2="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/04__FINAL_PROTEOMES/phex_1-20/*.fa"
P3="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/04__FINAL_PROTEOMES/phex_21-33/*.fa"
RVDB="/projet/fr2424/sib/dmckeown/db/virus/fastas/U-RVDBv20.0-prot.fasta"

###### A - get fa of all NCV per core gene from public databases ######
#srun grep ">" /projet/fr2424/sib/dmckeown/db/virus/fastas/NCVOG.fa > UBE07_A000; ## get NCVOG deflines
#srun cut -f 1,7 /projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/NCVOG_NCVOGdef | sed 's/NCVOG/ncvog=/g' > UBE07_A001; ## get list gi and NCVOG

### get gi list for 5 group I core genes: ###
#srun awk -F "\t" '$2 =="ncvog=0022" {print "gi|"$1"|"}' UBE07_A001 | grep -f - UBE07_A000 | sed 's/ /_/g' | sed 's/>//g' > UBE07_A002_0022_MCP;
#srun awk -F "\t" '$2 =="ncvog=0023" {print "gi|"$1"|"}' UBE07_A001 | grep -f - UBE07_A000 | sed 's/ /_/g' | sed 's/>//g' > UBE07_A002_0023_D5;
#srun awk -F "\t" '$2 =="ncvog=0038" {print "gi|"$1"|"}' UBE07_A001 | grep -f - UBE07_A000 | sed 's/ /_/g' | sed 's/>//g' > UBE07_A002_0038_polB;
#srun awk -F "\t" '$2 =="ncvog=0076" {print "gi|"$1"|"}' UBE07_A001 | grep -f - UBE07_A000 | sed 's/ /_/g' | sed 's/>//g' > UBE07_A002_0076_A18;
#srun awk -F "\t" '$2 =="ncvog=0249" {print "gi|"$1"|"}' UBE07_A001 | grep -f - UBE07_A000 | sed 's/ /_/g' | sed 's/>//g' > UBE07_A002_0249_A32;
#srun awk -F "\t" '$2 =="ncvog=0262" {print "gi|"$1"|"}' UBE07_A001 | grep -f - UBE07_A000 | sed 's/ /_/g' | sed 's/>//g' > UBE07_A002_0262_VLTF3;

#srun sed 's/ /_/g' /projet/fr2424/sib/dmckeown/db/virus/fastas/NCVOG.fa | seqkit grep -n -f UBE07_A002_0022_MCP - > UBE07_A003_0022_MCP.fa;
#srun sed 's/ /_/g' /projet/fr2424/sib/dmckeown/db/virus/fastas/NCVOG.fa | seqkit grep -n -f UBE07_A002_0023_D5 - > UBE07_A003_0023_D5.fa;
#srun sed 's/ /_/g' /projet/fr2424/sib/dmckeown/db/virus/fastas/NCVOG.fa | seqkit grep -n -f UBE07_A002_0038_polB - > UBE07_A003_0038_polB.fa;
#srun sed 's/ /_/g' /projet/fr2424/sib/dmckeown/db/virus/fastas/NCVOG.fa | seqkit grep -n -f UBE07_A002_0076_A18 - > UBE07_A003_0076_A18.fa;
#srun sed 's/ /_/g' /projet/fr2424/sib/dmckeown/db/virus/fastas/NCVOG.fa | seqkit grep -n -f UBE07_A002_0249_A32 - > UBE07_A003_0249_A32.fa;
#srun sed 's/ /_/g' /projet/fr2424/sib/dmckeown/db/virus/fastas/NCVOG.fa | seqkit grep -n -f UBE07_A002_0262_VLTF3 - > UBE07_A003_0262_VLTF3.fa;

#srun diamond blastp -d /projet/fr2424/sib/dmckeown/db/virus/U-RVDBv20.0-prot.dmnd -q UBE07_A003_0022_MCP.fa --more-sensitive -o UBE07_A003_0022_MCP.blast -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles;
#srun diamond blastp -d /projet/fr2424/sib/dmckeown/db/virus/U-RVDBv20.0-prot.dmnd -q UBE07_A003_0023_D5.fa --more-sensitive -o UBE07_A003_0023_D5.blast -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles;
#srun diamond blastp -d /projet/fr2424/sib/dmckeown/db/virus/U-RVDBv20.0-prot.dmnd -q UBE07_A003_0038_polB.fa --more-sensitive -o UBE07_A003_0038_polB.blast -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles;
#srun diamond blastp -d /projet/fr2424/sib/dmckeown/db/virus/U-RVDBv20.0-prot.dmnd -q UBE07_A003_0076_A18.fa --more-sensitive -o UBE07_A003_0076_A18.blast -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles;
#srun diamond blastp -d /projet/fr2424/sib/dmckeown/db/virus/U-RVDBv20.0-prot.dmnd -q UBE07_A003_0249_A32.fa --more-sensitive -o UBE07_A003_0249_A32.blast -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles;
#srun diamond blastp -d /projet/fr2424/sib/dmckeown/db/virus/U-RVDBv20.0-prot.dmnd -q UBE07_A003_0262_VLTF3.fa --more-sensitive -o UBE07_A003_0262_VLTF3.blast -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles;

#srun cut -f 13 UBE07_A003_0022_MCP.blast | sort -Vu | sed -E '/^.*viral metagenome.*$/d' | sed 's/ /_/g' > UBE07_A004_0022_MCP; srun sed 's/ /_/g' $RVDB | seqkit grep -n -f UBE07_A004_0022_MCP - > UBE07_A005_0022_MCP.fa;
#srun cut -f 13 UBE07_A003_0023_D5.blast | sort -Vu | sed -E '/^.*viral metagenome.*$/d' | sed 's/ /_/g' > UBE07_A004_0023_D5; srun sed 's/ /_/g' $RVDB | seqkit grep -n -f UBE07_A004_0023_D5 - > UBE07_A005_0023_D5.fa;
#srun cut -f 13 UBE07_A003_0038_polB.blast | sort -Vu | sed -E '/^.*viral metagenome.*$/d' | sed 's/ /_/g' > UBE07_A004_0038_polB; srun sed 's/ /_/g' $RVDB | seqkit grep -n -f UBE07_A004_0038_polB - > UBE07_A005_0038_polB.fa;
#srun cut -f 13 UBE07_A003_0076_A18.blast | sort -Vu | sed -E '/^.*viral metagenome.*$/d' | sed 's/ /_/g' > UBE07_A004_0076_A18; srun sed 's/ /_/g' $RVDB | seqkit grep -n -f UBE07_A004_0076_A18 - > UBE07_A005_0076_A18.fa;
#srun cut -f 13 UBE07_A003_0249_A32.blast | sort -Vu | sed -E '/^.*viral metagenome.*$/d' | sed 's/ /_/g' > UBE07_A004_0249_A32; srun sed 's/ /_/g' $RVDB | seqkit grep -n -f UBE07_A004_0249_A32 - > UBE07_A005_0249_A32.fa;
#srun cut -f 13 UBE07_A003_0262_VLTF3.blast | sort -Vu | sed -E '/^.*viral metagenome.*$/d' | sed 's/ /_/g' > UBE07_A004_0262_VLTF3; srun sed 's/ /_/g' $RVDB | seqkit grep -n -f UBE07_A004_0262_VLTF3 - > UBE07_A005_0262_VLTF3.fa;

### fix deflines to TAXON|ACCESSION and remove problematic characters:
#for file in UBE07_A005_*.fa; do srun sed -i -E 's/(>).+\|.+\|(.+[0-9]+\.[0-9]+).+\|.+\|.+\[(.+)\]/\1\3|\2/g' $file; done;
#for file in UBE07_A005_*.fa; do srun sed -i "s/ \|;\|:\|\,\|(\|)\|'\|\//_/g" $file; srun sed -i -E 's/_+/_/g' $file; done;

## remove duplicates by sequence
#for file in UBE07_A005_*.fa; do srun seqkit rmdup -s -o $(basename $file .fa).fa.reduced $file; done;

### MANUAL STEP - remove excessive numbers of highly similar sequences; especially those with multiple copies per genome
### MANUAL STEP - UBE07_A005_*.reduced.list;
### IN Dendroscope: inspect .reduced.fasttree and use 'List Selected Taxa' to select final list - save as UBE07_A005_*.reduced.list;
### use .reduced.list to get seqs from .fa.reduced
#for file in UBE07_A005_*.fa.reduced; do srun mafft --leavegappyregion --reorder --auto $file > $(basename $file .fa.reduced).reduced.aln; done;
#for file in UBE07_A005_*.reduced.aln; do srun FastTree $file > $(basename $file .reduced.aln).reduced.fasttree; done;

#for file in UBE07_A005_*.reduced.list; do seqkit grep -n -f $file $(basename $file .reduced.list).fa.reduced > UBE07_A006_$(basename $file .reduced.list).fa.reduced ; done;
#for file in UBE07_A006_UBE07_A005_*; do srun mv "$file" "${file//UBE07_A006_UBE07_A005_/UBE07_A006_}"; done;
### MANUAL STEP - add in any key taxa missing


###### B - get phex core gene homologs ######

### MANUAL STEP - get phex protein names of NCV core gene homologs; srun awk -F "\t" '$14 =="0022"' UBE05_B005_*_SUMMARY > UBE07_B000_0022_MCP;
### protein names to .fa:
#srun cut -f 4 UBE07_B000_0022_MCP | seqkit grep -f - $P1 $P2 $P3 > UBE07_B001_0022_MCP.fa;
#srun cut -f 4 UBE07_B000_0023_D5 | seqkit grep -f - $P1 $P2 $P3 > UBE07_B001_0023_D5.fa;
#srun cut -f 4 UBE07_B000_0038_polB | seqkit grep -f - $P1 $P2 $P3 > UBE07_B001_0038_polB.fa;
#srun cut -f 4 UBE07_B000_0076_A18 | seqkit grep -f - $P1 $P2 $P3 > UBE07_B001_0076_A18.fa;
#srun cut -f 4 UBE07_B000_0249_A32 | seqkit grep -f - $P1 $P2 $P3 > UBE07_B001_0249_A32.fa;
#srun cut -f 4 UBE07_B000_0262_VLTF3 | seqkit grep -f - $P1 $P2 $P3 > UBE07_B001_0262_VLTF3.fa;

### MANUAL STEP - delete seqs from UBE07_B001_*.fa that are too short
### check min seq lengths: srun seqkit stat UBE07_B001_*.fa;
### check test alignment UBE07_B001_*.fa > UBE07_B001_*.aln: for file in UBE07_B001_*.fa; do srun mafft --leavegappyregion --reorder --auto $file > $(basename $file .fa).aln; done;
### remove bad seqs UBE07_B001_*.fa and save as UBE07_B002_0022_MCP.fa;

### fix deflines to genome_contigN.gene.isoform and remove problematic characters:
#srun sed -i 's/mRNA\.\| assembled CDS//g' UBE07_B002_*.fa; ## fix phex deflines
#for file in UBE07_B002_*.fa; do srun sed -i "s/ \|;\|:\|\,\|(\|)\|'\|\//_/g" $file; srun sed -i -E 's/_+/_/g' $file; done;


###### C - alignment and phylogeny #######

## merge NCV and phex .fa
#srun cat UBE07_B002_0022_MCP.fa UBE07_A006_0022_MCP.fa.reduced > UBE07_C000_0022_MCP.fa;
#srun cat UBE07_B002_0023_D5.fa UBE07_A006_0023_D5.fa.reduced > UBE07_C000_0023_D5.fa;
#srun cat UBE07_B002_0038_polB.fa UBE07_A006_0038_polB.fa.reduced > UBE07_C000_0038_polB.fa;
#srun cat UBE07_B002_0076_A18.fa UBE07_A006_0076_A18.fa.reduced > UBE07_C000_0076_A18.fa;
#srun cat UBE07_B002_0249_A32.fa UBE07_A006_0249_A32.fa.reduced > UBE07_C000_0249_A32.fa;
#srun cat UBE07_B002_0262_VLTF3.fa UBE07_A006_0262_VLTF3.fa.reduced > UBE07_C000_0262_VLTF3.fa;

## ALIGNMENT. $1=input.fa
#srun mafft --leavegappyregion --reorder --auto $1 > $(basename $1 .fa).aln;

## Gblocks (may not work well for divergent sequences). $1=input.fa
#Gblocks $(basename $1 .fa).aln -b5=a -p=n;

## construct ML phylogenetic trees (run with 1 cpus per task)
srun raxml-ng --all --threads 1 --msa $(basename $1 .fa).aln --model LG+G8+F --bs-trees 200 --prefix $(basename $1 .fa).1; ## normal tree

#srun raxml-ng --threads 8 --all --msa $(basename $1 .fa).aln-gb --model LG+G8+F --bs-trees 200 --prefix $(basename $1 .fa).gb; ## Gblocks tree

#### old RaxML if -ng too slow
#raxmlHPC-PTHREADS -T 8 -d -f ae -m PROTGAMMAGTR -p $RANDOM -x $RANDOM -#200 -s $(basename $1 .fa).aln -n "$(basename $1 .fa)".raxml_old;








