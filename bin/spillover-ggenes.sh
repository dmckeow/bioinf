#!/bin/bash

## /panfs/jay/groups/27/dcschroe/dmckeow/data/Spillover_FINAL/phylogeny

conda activate bioinftools

#### run prodigal on fasta concatenated by binning
for Concat in *.fa; do 
  esl-translate $Concat > $Concat.translated_proteins.faa
  sed -i -E 's/>.+ source=(.+) coords=([0-9]+)\.\.([0-9]+) .*/>\1__coords\2__\3/g' $Concat.translated_proteins.faa
done

### run interproscan
rm -f run_iprscan; touch run_iprscan
Q="sbatch --time=12:00:00 --cpus-per-task=12 --mem=64GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='module load interproscan;"
R="'"
for f in *.translated_proteins.faa; do 
  OutName=$(basename $f .faa); echo "$Q sed 's/\*/X/g' ${OutName}.faa | interproscan.sh --iprlookup -t p -f GFF3 -i - -o ${OutName}.iprscan; sed -i -e '/^##FASTA/q' -e '/^#/d' ${OutName}.iprscan; sed -i -e '/^##FASTA/d' ${OutName}.iprscan$R" >> run_iprscan; 
done

bash run_iprscan

awk '/\tPfam\t/' *.translated_proteins.iprscan > all.translated_proteins.iprscan.all
sed -i 's/;/\t/g' all.translated_proteins.iprscan.all
sed -i -E 's/ +/_/g' all.translated_proteins.iprscan.all
sed -i -E 's/__coords([0-9]+)__([0-9]+)\t/\t\1\t\2\t/g' all.translated_proteins.iprscan.all

awk -F "\t" 'BEGIN{FS=OFS="\t"} $2 > $3 {temp=$2; $2=$3; $3=temp; $9="-"} 1' all.translated_proteins.iprscan.all | awk 'BEGIN{OFS="\t"} $9 == "-" {$9="0"} $9 == "+" {$9="1"} 1' > tmp && mv tmp all.translated_proteins.iprscan.all

sed -i -E 's/date=|Target=|ID=|signature_desc=|Name=|status=|Dbxref=//g' all.translated_proteins.iprscan.all
 
sed -i -z 's/^/molecule\tstart\tend\tsource\tprotein_match\tstartaa\tendaa\tevalue\torientation\tframe\tdate\ttarget\tID\tsignature_desc\tName\tstatus\tDbxref\n/1' all.translated_proteins.iprscan.all

## now get he gene preditions
cat *.predicted_proteins.gff > all.predicted_proteins.all

sed -i 's/;/\t/g' all.predicted_proteins.all

awk 'BEGIN{OFS="\t"} $7 == "-" {$7="0"} $7 == "+" {$7="1"} 1' all.predicted_proteins.all > tmp && mv tmp all.predicted_proteins.all

sed -i -E 's/ID=|partial=|start_type=|rbs_motif=|rbs_spacer=|gc_cont=|conf=|score=|cscore=|sscore=|rscore=|uscore=|tscore=//g' all.predicted_proteins.all

sed -i -z 's/^/molecule\tsource\ttype\tstart\tend\tgene_score\torientation\tframe\tID\tpartial\tstart_type\trbs_motif\trbs_spacer\tgc_cont\tconf\tscore\tcscore\tsscore\trscore\tuscore\ttscore\n/1' all.predicted_proteins.all