#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --cpus-per-task 8
#SBATCH --mem 240GB

######### databases
NT="/db/nt/current/blast/nt"
NR="/db/nr/current/blast/nr"

####################### INPUTS #####################################
### commnad line arguments
## $1 input fasta query protein, NCLDVs
## $2 gene name (must match info in IVEX002 gffs)

## files
PHAA="/shared/projects/phaeoexplorer_virus/phaeoex_screen/input/IVEX000/*.faa"
SCR="/scratch2/fr2424/sib/dmckeown/"
PHGF="/shared/projects/phaeoexplorer_virus/phaeoex_screen/finalresult/IVEX002/*.gff.gz"

###### load software ######
module unload blast
module load blast; module load seqkit

n1=$(echo "$1" | sed 's/.*\///g' | cut -d "." -f 1)

### make database of local files (phaeoexplorer) (with filename included in deflines)
#for f in $PHAA; do g="$(basename "$f" .faa | sed -E 's/^gms[0-9]+_//g')"; cat $f | awk -v h=$g '{gsub(">",">"h"|")}1'; done > "$SCR"tmp.fa
#makeblastdb -in "$SCR"tmp.fa -title "Phaeoexplorer proteomes" -dbtype prot -out "$SCR"phaa
PHAADB=""$SCR"phaa"

##### get all homologs to gene of interest from phex dataset
#for f in $PHGF; do n="$(basename $f .gff.gz)"; zcat $f | awk -F "\t|;" -v a=$2 '{if($0 ~ a) print $9}' | sed 's/ID=//g' | seqkit grep -f - /shared/projects/phaeoexplorer_virus/phaeoex_screen/input/IVEX000/gms*"$n".faa | awk -v n=$n '{if($0 ~ ">") print ">"n"|"$0; else print $0}' | sed -e 's/>//2' -e '/^_.*\.faa$/d' ; done > tmp_"$n1".faa
#cat tmp_"$n1".faa $1 > tmp_"$n1"_query.faa ## now we have a query fasta comprised of the NCLDV homologs AND the cellular homologs in the brown algal phex dataset
#grep ">" tmp_"$n1"_query.faa > "$n1".query.list ### keep deflines for reference
## this should return any and all viral/cellular homologs from nr with psiblast

############### input VS ALL Phaeoexplorer proteomes
######## to maximise possible homologs retreived from phex
#### psiblast with iterations
#psiblast -query tmp_"$n1"_query.faa -db $PHAADB -out "$n1".phaa.psiblast -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles" -max_target_seqs 5 -num_threads 12 -num_iterations 3

#cut -f 2 "$n1".phaa.psiblast | sort -Vu > tmp_"$n1"_list ## get nr list of seq IDs of all matches
### get fasta for matches, merge with queries, remove duplicates by name
#seqkit grep -f tmp_"$n1"_list "$SCR"tmp.fa | cat $1 - | seqkit rmdup - -n -o "$n1".phaa.psiblast.faa

############### input and Phaeoexplorer proteome matches VS ALL NR
#### psiblast with iterations
psiblast -query "$n1".phaa.psiblast.faa -db $NR -out "$n1".nr.psiblast -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles" -max_target_seqs 5 -num_threads 8 -num_iterations 3

#cat $1 "$n1"_psiblast.fa > tmp ## concatenate query fasta and match fasta
#seqkit rmdup -s tmp -o "$n1".phaa.psiblast.nr.fa ## remove duplicated sequences


##blastdbcmd -db $PHAADB -entry_batch tmp_list -outfmt “%f” -target_only -out "$n1".phaa.psiblast.faa ## get fasta of all matches in phex

#rm -f tmp*
