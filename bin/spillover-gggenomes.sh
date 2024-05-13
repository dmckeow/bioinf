#!/bin/bash

conda activate bioinftools

#### input variables
#InFna="/panfs/jay/groups/27/dcschroe/dmckeow/gggenomes/data-raw/emales/emales.fna" ## input fasta nucleotide
InFna="/home/dcschroe/dmckeow/data/house_ref_genomes/multiple_DWV.fna" ## input fasta nucleotide
DmndDb="/panfs/jay/groups/27/dcschroe/shared/bioinfdb/DMND/vog.dmnd" ## database for diamond

OutName=$(basename ${InFna%.*}) ## auto set output names based on input fasta basename

##BlastSubject1="/panfs/jay/groups/27/dcschroe/dmckeow/gggenomes/data-raw/emales/mavirus.faa"
rm -fr tmp_${OutName}; mkdir tmp_${OutName}

# Annotate genes - output needed is a gff file
prodigal -n -t ${OutName}.train -i $InFna
prodigal -t ${OutName}.train -i $InFna -o ${OutName}.gff -f gff -a ${OutName}.faa
sed -i -e '/>/ s/;.*//g' -e 's/.* ID=/>/g' ${OutName}.faa

# six frame translation
esl-translate $InFna > ${OutName}-tran.faa
sed -i 's/|/,/g' ${OutName}-tran.faa
sed -i 's/ /|/g' ${OutName}-tran.faa

# split into one genome per file
rm -fr ${OutName}.split
seqkit split -O ${OutName}.split -i $InFna  

# self-align opposite strands    

for f in ${OutName}.split/*; do
  minimap2 -c -B5 -O6 -E3 --rev-only $f $f > tmp_${OutName}/$(basename ${f%.*}).paf
done

cat tmp_${OutName}/*.paf > ${OutName}-tirs.paf


# All-vs-all alignment for synteny regions
minimap2 -X -N 50 -p 0.1 -c $InFna $InFna > ${OutName}.paf


######## get GC content
samtools faidx $InFna
cp ${InFna}.fai ${OutName}.fai

seqkit faidx $InFna
cp ${InFna}.fai ${OutName}.seqkit.fai

cut -f 1,2 ${InFna}.fai > tmp_${OutName}/${OutName}.fai
WINDOW=50
while IFS=$'\t' read -r sequence length; do
    for ((start=1; start<=length; start+=WINDOW)); do
        end=$((start+WINDOW-1))
        if [ "$end" -gt "$length" ]; then
            end="$length"
        fi
        echo -e "$sequence\t$start\t$end"
    done
done < tmp_${OutName}/${OutName}.fai > tmp_${OutName}/${OutName}.bed


bedtools nuc -fi $InFna -bed tmp_${OutName}/${OutName}.bed | cut -f 1-3,5 | sed 's/\t/\tGC.w50\t/3' | sed '1,1d' > ${OutName}-gc.tsv

mmseqs easy-cluster ${OutName}.faa ${OutName}-mmseqs ./tmp_${OutName} -e 1e-5 -c 0.7
mmseqs easy-cluster ${OutName}-tran.faa ${OutName}-mmseqs-tran ./tmp_${OutName}-tran -e 1e-5 -c 0.7

### give seqs ortholog group names based on mmseqs preds
awk -F "\t" '{if($1 == $2) print $1}' ${OutName}-mmseqs_cluster.tsv > tmp_${OutName}/${OutName}-cogs.tsv
counter=1
while read -r line; do
	group=$(printf '%03d' $counter)
    echo "/^$line\\t/ s/$/\tcog$group/g"
    ((counter++))
done < tmp_${OutName}/${OutName}-cogs.tsv | sed -f - ${OutName}-mmseqs_cluster.tsv > ${OutName}-cogs.tsv
cut -f 3 ${OutName}-cogs.tsv | sort -V | uniq -c | sed -E 's/^ +(.+) (.+)/\/\\t\2\$\/ s\/$\/\\t\1\/g/g' | sed -f - ${OutName}-cogs.tsv | cut -f 2-4 > tmp_${OutName}/${OutName}-cogs.tsv && mv tmp_${OutName}/${OutName}-cogs.tsv ${OutName}-cogs.tsv

### give seqs ortholog group names based on mmseqs translations
awk -F "\t" '{if($1 == $2) print $1}' ${OutName}-mmseqs-tran_cluster.tsv > tmp_${OutName}/${OutName}-tran-cogs.tsv
counter=1
while read -r line; do
  group=$(printf '%03d' $counter)
    echo "/^$line\\t/ s/$/\tcog$group/g"
    ((counter++))
done < tmp_${OutName}/${OutName}-tran-cogs.tsv | sed -f - ${OutName}-mmseqs-tran_cluster.tsv > ${OutName}-tran-cogs.tsv
cut -f 3 ${OutName}-tran-cogs.tsv | sort -V | uniq -c | sed -E 's/^ +(.+) (.+)/\/\\t\2\$\/ s\/$\/\\t\1\/g/g' | sed -f - ${OutName}-tran-cogs.tsv | cut -f 2-4 > tmp_${OutName}/${OutName}-tran-cogs.tsv && mv tmp_${OutName}/${OutName}-tran-cogs.tsv ${OutName}-tran-cogs.tsv


############# blast against database of interest
# mavirus.faa - published
##blastp -query ${OutName}.faa -subject $BlastSubject1 -outfmt 7 > ${OutName}_mavirus-blastp.tsv
##perl -ne 'if(/>(\S+) gene=(\S+) product=(.+)/){print join("\t", $1, $2, $3), "\n"}' \${BlastSubject1} > $(basename $BlastSubject1).tsv

########## interproscan
### amino acid, predictions
Q="sbatch --time=48:00:00 --cpus-per-task=12 --mem=48GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='module load interproscan;"
R="'"
echo "$Q sed 's/\*/X/g' ${OutName}.faa | interproscan.sh -t p -f TSV -i - -o ${OutName}.iprscan$R" > run_iprscan

### amino acid, frame translations
Q="sbatch --time=48:00:00 --cpus-per-task=12 --mem=48GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='module load interproscan;"
R="'"
echo "$Q sed 's/\*/X/g' ${OutName}-tran.faa | interproscan.sh -t p -f TSV -i - -o ${OutName}-tran.iprscan$R" >> run_iprscan

########## diamond
### amino acid, predictions
Q="sbatch --time=24:00:00 --cpus-per-task=12 --mem=48GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval \"\$(conda shell.bash hook)\"; conda activate bioinftools;"
R="'"
echo "$Q diamond blastp --evalue 1e-20 --sensitive -d $DmndDb -q ${OutName}.faa -f 6 -o ${OutName}.dmnd.blastp$R" > run_diamond

### amino acid, frame translations
Q="sbatch --time=24:00:00 --cpus-per-task=12 --mem=48GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval \"\$(conda shell.bash hook)\"; conda activate bioinftools;"
R="'"
echo "$Q diamond blastp --evalue 1e-20 --sensitive -d $DmndDb -q ${OutName}-tran.faa -f 6 -o ${OutName}-tran.dmnd.blastp$R" >> run_diamond

######### submit search jobs
bash run_iprscan
bash run_diamond


############ ONCE the search jobs are done:

### turn interproscan into mock blast table
for f in ${OutName}*.iprscan; do
  sed -E 's/:+/:/g' $f | awk -F "\t" 'BEGIN{OFS="\t"}; {print $1,$4"::"$5"::"$6,"0",$3,"0","0",$7,$8,"0","0",$9,"0"}' > ${f}.pseudoblast
done

### filter blast results to remove redundant hits
for f in ${OutName}*.dmnd.blast*; do
  awk -F "\t" '!a[$1,$7,$8]++' $f > tmp_${OutName}/${OutName}.dmnd.blast && mv tmp_${OutName}/${OutName}.dmnd.blast $f ## get rid of hist to same coords, keeping best hit

  awk -F "\t" 'BEGIN{OFS="\t"};{ gsub(/\..+/, "", $2) } 1' $f > tmp_${OutName}/${OutName}.dmnd.blast && mv tmp_${OutName}/${OutName}.dmnd.blast $f ## VOG specific
done

