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


# example_genes
#   molecule  gene  start    end  strand orientation source
#1   Genome1  genA  15389  17299 reverse           1  
#2   Genome1  genB  17301  18161 forward           0
#3   Genome1  genC  18176  18640 reverse           1
#4   Genome1  genD  18641  18985 forward           0
#5   Genome1  genE  18999  20078 reverse           1
#6   Genome1  genF  20086  20451 forward           1
#7   Genome1 protC  22777  22989 forward           1


#### split interproscan results and get protein fastas
rm -fr iprscansplit; mkdir iprscansplit

rm -f iprscansplit/*.faa
rm -f iprscansplit/*.fai

for f in *.predicted_proteins.faa; do
  cat $f >> iprscansplit/all.faa
done

seqkit rmdup iprscansplit/all.faa > iprscansplit/all.faa.tmp && mv iprscansplit/all.faa.tmp iprscansplit/all.faa

for f in *.fa.predicted_proteins.iprscan; do
  outdir="iprscansplit/$(basename $f .fa.predicted_proteins.iprscan).iprsplitdir"
  mkdir $outdir
awk -F "\t|;" 'BEGIN{OFS=";"}; {if($2 == "Pfam") print $1,$2,$4,$5,$12,$13}' $f > $outdir/$(basename $f)
sed -E -e 's/signature_desc=//g' -e 's/Name=//g' -e 's/[^a-zA-Z0-9._;-]/_/g' $outdir/$(basename $f) | \
sort -Vu \
> $outdir/$(basename $f).tmp && mv $outdir/$(basename $f).tmp $outdir/$(basename $f)

rm -f $outdir/*.iprprotlist

for g in $outdir/*iprscan; do
  while IFS=";" read -r contig source start end signature_desc name; do
    echo -e "$contig\t$start\t$end" >> $outdir/"$source"__"$signature_desc"__"$name".iprprotlist
  done < "$g"
done

for h in $outdir/*.iprprotlist; do
  bedtools getfasta -fi iprscansplit/all.faa -bed $h >> $outdir/$(basename $h .iprprotlist).faa
done

done #######

########## remove seqs shorter than half the average length for each protein type
for f in iprscansplit/*.iprsplitdir/*.faa; do
  AVGLEN=$(seqkit stat -T $f | cut -f 7 | sed '1,1d' | sed 's/\..*//g')
  seqkit seq -m $(($AVGLEN /3)) $f > ${f}.tmp && mv ${f}.tmp ${f}
done

######### check what protein groups you have and decide what will be concatenated for phylogeny and networks
realpath iprscansplit/*.iprsplitdir/*.faa > iprscansplit/list_fastas
ls iprscansplit/*.iprsplitdir/*.faa | sed 's/.*\///g' | sort -Vu > iprscansplit/list_fastas_unique


#Pfam__Bunyavirus_RNA_dependent_RNA_polymerase__PF04196.faa
#Pfam__DNA-directed_RNA_polymerase_N-terminal__PF14700.faa
#Pfam__Mononegavirales_RNA_dependent_RNA_polymerase__PF00946.faa
#Pfam__RNA_dependent_RNA_polymerase__PF00680.faa
#Pfam__RNA_dependent_RNA_polymerase__PF00978.faa
#Pfam__Reovirus_RNA-dependent_RNA_polymerase_lambda_3__PF07925.faa
#Pfam__Viral_RNA_dependent_RNA_polymerase__PF00998.faa

#### all polymerases, all virus groups
rm -f iprscansplit/all_RdRp.faa
for f in $(grep "polymerase" iprscansplit/list_fastas); do cat $f >> iprscansplit/all_RdRp.faa; done

### lists faas to run separaely by virus group
grep "polymerase" iprscansplit/list_fastas > iprscansplit/list_fastas_RdRp


### group by Pfam group

# Define the file containing the list of paths
rm -f iprscansplit/list_fastas_RdRp.*.Pfam
# Loop through each line in the file
while IFS= read -r line; do
    # Extract the basename using the basename command
    basename=$(basename "$line")
    
    # Print the basename
    echo "$line" >> iprscansplit/list_fastas_RdRp."$basename".Pfam
done < "iprscansplit/list_fastas_RdRp"

cd iprscansplit


sbatch --time=96:00:00 --cpus-per-task=8 --mem=64GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; bioinf-phylogeny.py -i all_RdRp.faa -N -p all_RdRp'

####### run alignments for proteins
for f in $(cat list_fastas_RdRp); do
    OUT=$(echo $f | sed -e 's/.*iprscansplit\///g' -e 's/.iprsplitdir\//__/g' -e 's/\.faa//g')
  sbatch --time=96:00:00 --cpus-per-task=8 --mem=64GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; bioinf-phylogeny.py -i "'$f'" -N -p "'$OUT'"'
done


for f in list_fastas_RdRp.*.Pfam; do
  OUT=$(basename $f .faa.Pfam | sed 's/list_fastas_RdRp.//g')
  sbatch --time=96:00:00 --cpus-per-task=8 --mem=64GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; bioinf-phylogeny.py -i "'$f'" -N -p "'$OUT'"'
done

##########################################

########## generate similarity scores for netowrks

sbatch --time=96:00:00 --cpus-per-task=4 --mem=32GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; blastp -query all_RdRp.faa -subject all_RdRp.faa -outfmt 6 > all_RdRp.blastp'

for f in $(cat list_fastas_RdRp); do
    OUT=$(echo $f | sed -e 's/.*iprscansplit\///g' -e 's/.iprsplitdir\//__/g' -e 's/\.faa//g')
  sbatch --time=96:00:00 --cpus-per-task=4 --mem=32GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; blastp -query "'$f'" -subject "'$f'" -outfmt 6 > "'$OUT'".blastp'
done

rm -f *.blastp.input.tmp
for f in list_fastas_RdRp.*.Pfam; do
    OUT=$(basename $f .faa.Pfam | sed 's/list_fastas_RdRp.//g')
  for g in $(cat $f); do
    cat $g >> ${OUT}.blastp.input.tmp
  done
done

for f in list_fastas_RdRp.*.Pfam; do
    OUT=$(basename $f .faa.Pfam | sed 's/list_fastas_RdRp.//g')
  sbatch --time=96:00:00 --cpus-per-task=4 --mem=32GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; blastp -query "'${OUT}.blastp.input.tmp'" -subject "'${OUT}.blastp.input.tmp'" -outfmt 6 > "'$OUT'".blastp'
done

rm -f *.blastp.input.tmp

for f in *blastp; do
  NAME=$(basename $f .blastp)
  awk -F "\t" '!a[$1,$2]++' $f | awk -F "\t" '{if($1 != $2) print $0}' > ${f}.tmp && mv ${f}.tmp ${f}
  sed -z -i 's/^/from\tto\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n/g' $f
done

find . -type f -empty -name "*.blastp" -delete


##########################################
### tidy up the alignments 

find . -type f -empty -name "*.aln" -delete

rm -f *.trimal.aln
for f in *.aln; do
  trimal -in $f -out $(basename $f .aln).trimal.aln -keepheader -automated1 -htmlout $(basename $f .aln).trimal.html
done

for f in *.trimal.aln; do
  mv $f $(basename $f .trimal.aln).aln
done


####### run trees for proteins

sbatch --time=96:00:00 --cpus-per-task=8 --mem=64GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; bioinf-phylogeny.py -i all_RdRp.faa -P -p all_RdRp'

####### run alignments for proteins
for f in $(cat list_fastas_RdRp); do
    OUT=$(echo $f | sed -e 's/.*iprscansplit\///g' -e 's/.iprsplitdir\//__/g' -e 's/\.faa//g')
  sbatch --time=96:00:00 --cpus-per-task=8 --mem=64GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; bioinf-phylogeny.py -i "'$f'" -P -p "'$OUT'"'
done


for f in list_fastas_RdRp.*.Pfam; do
  OUT=$(basename $f .faa.Pfam | sed 's/list_fastas_RdRp.//g')
  sbatch --time=96:00:00 --cpus-per-task=8 --mem=64GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; bioinf-phylogeny.py -i "'$f'" -P -p "'$OUT'"'
done


