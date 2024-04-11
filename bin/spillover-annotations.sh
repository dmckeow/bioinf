#!/bin/bash

conda activate bioinftools



#### run prodigal on fasta concatenated by binning
for Concat in *.fa; do 
  prodigal -i $Concat -a $Concat.predicted_proteins.faa -o $Concat.predicted_proteins.gff -g 11 -q -f gff -p meta; 
done

### run interproscan
rm -f run_iprscan; touch run_iprscan
Q="sbatch --time=8:00:00 --cpus-per-task=12 --mem=64GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='module load interproscan;"
R="'"
for f in *.faa; do 
  OutName=$(basename $f .faa); echo "$Q sed 's/\*/X/g' ${OutName}.faa | interproscan.sh --iprlookup -t p -f GFF3 -i - -o ${OutName}.iprscan$R" >> run_iprscan; 
done

bash run_iprscan

### get rid of the stupid fucking fasta file that interproscan appends to fucking the gff file
sed -i -e '/^##FASTA/q' -e '/^#/d' *.iprscan
