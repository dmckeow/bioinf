#!/bin/bash

module purge;
module load bedtools

OUTNAME=$2

##########################################
if [[ $1 == "viralrecall" ]]; then
##########################################

## running viralrecall.py (redo, bigger window) - tried -w 30 needs to be bigger window still
## window good but redo with smaller min size (2 kb) and viral hits #

###sbatch --time=96:00:00 --cpus-per-task=12 --mem=240GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; python viralrecall.py -i ALL_first_bincontigs.fna.split -p ALL_first_bincontigs -w 50 -g 1 -t 12 -b -f -r -m 2'


ls ALL_first_bincontigs2.fa.split/* > tmp.list; split -d -n l/20 tmp.list ALL_first_bincontigs2.fa.split.list.

for f in ALL_first_bincontigs2.fa.split.list.*; do
  name=$(basename $f | sed -E 's/list\.//g')
  mkdir $name
  for g in $(cat $f); do
    cp $g $name/
  done
done

rm -fr ALL_first_bincontigs2.fa.split
rm -f ALL_first_bincontigs2.fa.split.list.*

##### make fna file extension

for f in ALL_first_bincontigs2.fa.split.[0-9][0-9]; do
  cd $f
  for i in *.fa ; do mv "$i" "${i/.fa/.fna}" ; done
  cd ..
done

###### RUN version 2 - with 4951 contigs (higher pre-filter evalue threshold)
for f in $(find -type d -name "ALL_first_bincontigs2.fa.split.*"); do
  sq="'"
echo -e "sbatch --time=96:00:00 --cpus-per-task=12 --mem=240GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap=${sq}eval \"\$(conda shell.bash hook)\"; conda activate bioinftools; python viralrecall.py -i $f -p ViralRecall_$(basename $f) -w 50 -g 1 -t 12 -b -f -m 2${sq}"
done > run_viral_recall_ALL_first_bincontigs2

bash run_viral_recall_ALL_first_bincontigs2


##########################################
fi
##########################################

##########################################
if [[ $1 == "hmmsearch" ]]; then
##########################################

#### AT THE SAME TIME as the viral recall runs on the gene predictions, we need to 6 frame translation of the contigs and do another hmmsearch vs vogdb/gvog to identify potential hidden frameshifted genes

for f in $(find -type d -name "ALL_first_bincontigs2.fa.split.[0-9][0-9]"); do
  sq="'"
  rm -fr ${f}-translated
  mkdir ${f}-translated
echo -e "sbatch --time=8:00:00 --cpus-per-task=8 --mem=64GB --partition aglarge -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap=${sq}eval \"\$(conda shell.bash hook)\"; conda activate bioinftools; for f in ${f}/*.fna; do esl-translate \$f > ${f}-translated/\$(basename \$f .fna)-translated.faa; done${sq}"
done > run_translate_ALL_first_bincontigs2-translated

bash run_translate_ALL_first_bincontigs2-translated

############

ls ALL_first_bincontigs.fna.split-translated/*-translated.faa | split -d -n l/20 - fasta_list.1.tmp.

for f in $(find -type d -name "ALL_first_bincontigs2.fa.split.[0-9][0-9]-translated" | sort -V); do
  sq="'"
echo -e "sbatch --time=96:00:00 --cpus-per-task=4 --mem=64GB --partition ag2tb -o /scratch.global/dcschroe/dmckeow/slurm.%N.%j.out -e /scratch.global/dcschroe/dmckeow/slurm.%N.%j.err --wrap=${sq}eval \"\$(conda shell.bash hook)\"; conda activate bioinftools; for f in ${f}/*.faa; do output=\$(echo \$f | sed 's/\.faa//g'); hmmsearch -E 1e-10 --cpu 2 --tblout \${output}.vogout hmm/gvog.hmm \$f; hmmsearch -E 1e-10 --cpu 4 --tblout \${output}.markerout hmm/NCLDV_markers.hmm \$f; done${sq}"
done > run_hmmsearch_ALL_first_bincontigs2-gvog_marker-translated

bash run_hmmsearch_ALL_first_bincontigs2-gvog_marker-translated

## also pfam
for f in $(find -type d -name "ALL_first_bincontigs2.fa.split.[0-9][0-9]-translated" | sort -V); do
  sq="'"
echo -e "sbatch --time=96:00:00 --cpus-per-task=4 --mem=64GB --partition ag2tb -o /scratch.global/dcschroe/dmckeow/slurm.%N.%j.out -e /scratch.global/dcschroe/dmckeow/slurm.%N.%j.err --wrap=${sq}eval \"\$(conda shell.bash hook)\"; conda activate bioinftools; for f in ${f}/*.faa; do output=\$(echo \$f | sed 's/\.faa//g'); hmmsearch -E 1e-10 --cpu 4 --tblout \${output}.pfamout hmm/pfam.hmm \$f; done${sq}"
done > run_hmmsearch_ALL_first_bincontigs2-pfam-translated

##### BEFORE RUNNING, we subsample these two subsets, because they were taknig way too long for their jobs to finish within the 3 days run time limit of MSI
for f in ./ALL_first_bincontigs2.fa.split.11-translated/*.faa; do seqkit seq -m 50 $f > ${f}.tmp && mv ${f}.tmp $f; done
for f in ./ALL_first_bincontigs2.fa.split.13-translated/*.faa; do seqkit seq -m 75 $f > ${f}.tmp && mv ${f}.tmp $f; done

bash run_hmmsearch_ALL_first_bincontigs2-pfam-translated

##########################################
fi
##########################################

##########################################
if [[ $1 == "part1" ]]; then
##########################################


############### WAIT FOR ALL sampels to be done with runnign viralrecall and hmm before continuing
###### get all the viral regions fixed

rm -f ALL_summary.tsv; touch ALL_summary.tsv
for f in $(find ViralRecall_ALL_first_bincontigs2.fa.split.*/ALL_first_bincontigs2.part_* -name "*.summary.tsv"); do cat $f >> ALL_summary.tsv ; done


head -1 ALL_summary.tsv | sed 's/\t/;/g' > ALL_summary.header.tmp
sed -i '/^viral_regions\treplicon/d' ALL_summary.tsv
sed -i 's/;/,/g' ALL_summary.tsv
sed -i 's/\t/;/g' ALL_summary.tsv

awk 'BEGIN {FS=OFS=";"} {gsub("\\.[0-9]+", "", $3)} 1' ALL_summary.tsv > tmp && mv tmp ALL_summary.tsv
awk 'BEGIN {FS=OFS=";"} {gsub("\\.[0-9]+", "", $4)} 1' ALL_summary.tsv > tmp && mv tmp ALL_summary.tsv
awk 'BEGIN {FS=OFS=";"} {gsub("\\.[0-9]+", "", $5)} 1' ALL_summary.tsv > tmp && mv tmp ALL_summary.tsv
awk 'BEGIN {FS=OFS=";"} {gsub("\\.[0-9]+", "", $6)} 1' ALL_summary.tsv > tmp && mv tmp ALL_summary.tsv
awk 'BEGIN {FS=OFS=";"} {gsub("\\.[0-9]+", "", $8)} 1' ALL_summary.tsv > tmp && mv tmp ALL_summary.tsv
awk 'BEGIN {FS=OFS=";"} {gsub("\\.[0-9]+", "", $9)} 1' ALL_summary.tsv > tmp && mv tmp ALL_summary.tsv

sort -t ";" -k 2,2V -k 3,3g -k 4,4g ALL_summary.tsv -o ALL_summary.tsv

### get the inter region distances

# Read the file line by line
while IFS=';' read -r viral_regions contig start stop vregion_length contig_length score num_viralhits num_ORFs markers; do
  if [[ "$contig" == "$prev_contig" ]]; then
      echo "$viral_regions;$contig;$start;$stop;$vregion_length;$contig_length;$score;$num_viralhits;$num_ORFs;$markers;$(($start-$prev_stop))"
  elif [[ "$contig" != "$prev_contig" ]]; then
    echo "$viral_regions;$contig;$start;$stop;$vregion_length;$contig_length;$score;$num_viralhits;$num_ORFs;$markers;$(($start-1))"
  fi

  # Store the current line's values as previous for the next iteration
  prev_contig="$contig"
  prev_start="$start"
  prev_stop="$stop"

done < ALL_summary.tsv > ALL_summary.1.tmp


sort -t ";" -k 2,2V -k 3,3gr -k 4,4gr ALL_summary.1.tmp -o ALL_summary.1.tmp


# Read the file line by line
while IFS=';' read -r viral_regions contig start stop vregion_length contig_length score num_viralhits num_ORFs markers distance_from_upstream_vregion; do
  if [[ "$contig" == "$prev_contig" ]]; then
      echo "$viral_regions;$contig;$start;$stop;$vregion_length;$contig_length;$score;$num_viralhits;$num_ORFs;$markers;$distance_from_upstream_vregion"
  elif [[ "$contig" != "$prev_contig" ]]; then
    echo "$viral_regions;$contig;$start;$stop;$vregion_length;$contig_length;$score;$num_viralhits;$num_ORFs;$markers;$distance_from_upstream_vregion;$(($contig_length-$stop))"
  fi

  # Store the current line's values as previous for the next iteration
  prev_contig="$contig"
  prev_start="$start"
  prev_stop="$stop"

done < ALL_summary.1.tmp | sort -t ";" -k 2,2V -k 3,3g -k 4,4g > ALL_summary.2.tmp


rm -f tmp.viralregions.1.*
awk 'BEGIN {FS=OFS=";"} {print $11"C,"$5$1","$12"C" > "tmp.viralregions.1."$2}' ALL_summary.2.tmp

rm -f tmp.ALL.viralregions.by_contig; touch tmp.ALL.viralregions.by_contig
for f in tmp.viralregions.1.*; do
  sed 's/,C/,/g' $f | tr -d '\n' | sed 's/$/\n/1' | awk -v f=$f '{print f";"$0}' | sed 's/tmp\.viralregions\.1\.//g'
done >> tmp.ALL.viralregions.by_contig
sed -i 's/viral_region_/vr/g' tmp.ALL.viralregions.by_contig

rm -f tmp.viralregions.1.*




####################################
####################################


### prepare file to convert GVOG to NCVOG
## if a GVOG has no corresponding NCVOG, then the GVOG will be retained
cut -f 1,3 hmm/gvog_annotation.tsv | sed -e '1,1d' -e 's/\t$/\tno_annot/g' -e 's/ | /,/g' | sort -Vu | awk -F "\t" '{if($2 == "no_annot") print $1"\t"$1; else print $0}' | awk -F "\t" '{print "s/^"$1"$/"$2"/g"}' > GVOG_to_NCVOG.tmp

## parse the vog hmm hits to identify which originate from possible frameshifted or intronised proteins

### from prodigal gene predictions
rm -f ALL_first_bincontigs.vogout; touch ALL_first_bincontigs.vogout; for f in $(find ViralRecall_ALL_first_bincontigs2.fa.split.*/ALL_first_bincontigs2.part_* -name "*.vogout" | sed '/-translated.vogout/d'); do cat $f >> ALL_first_bincontigs.vogout ; done

## fix delimiters and reduce to best vog hit per query gene
sed -E -e '/^#/d' -e 's/ # |;/\t/g' -e 's/ +/\t/g' ALL_first_bincontigs.vogout | sort -t $'\t' -k 1,1V -k 5,5g -k 6,6gr | awk -F "\t" '!a[$1]++' | cut -f 1,3,19,20 | sed -E 's/(.+[0-9]+)_([0-9]+\t)(.+)/\1\t\2\3/1' > tmp && mv tmp ALL_first_bincontigs.vogout

sort -t $'\t' -k 1,1V -k 4,5n ALL_first_bincontigs.vogout -o ALL_first_bincontigs.vogout
sed -i 's/\t/;/g' ALL_first_bincontigs.vogout

### now we need to add the NCVOG because there can be multiple GVOGs for a single NCVOG
cut -d ";" -f 3 ALL_first_bincontigs.vogout > ALL_first_bincontigs.vogout.tofix.tmp
sort -Vu ALL_first_bincontigs.vogout.tofix.tmp | grep -f - GVOG_to_NCVOG.tmp | sed -f - ALL_first_bincontigs.vogout.tofix.tmp | paste -d ";" ALL_first_bincontigs.vogout - > tmp && mv tmp ALL_first_bincontigs.vogout

#####################
## NOW prepare for the -translated hits

rm -f ALL_first_bincontigs-translated.vogout; touch ALL_first_bincontigs-translated.vogout; for f in $(find ALL_first_bincontigs2.fa.split.*-translated -name "*-translated.vogout"); do cat $f >> ALL_first_bincontigs-translated.vogout ; done

### teh sort here is important because it keeps only the best vog hit per orf
sed -E -e '/^#/d' -e 's/ [a-z]+=|\.\./\t/g' -e 's/ +/\t/g' ALL_first_bincontigs-translated.vogout | sort -t $'\t' -k 19,19V -k 1,1V -k 5,5g -k 6,6gr | awk -F "\t" '!a[$1,$19]++' | awk 'BEGIN {FS=OFS="\t"} {print $19,$1,$3,$20,$21}' > tmp && mv tmp ALL_first_bincontigs-translated.vogout

sort -t $'\t' -k 1,1V -k 4,5n ALL_first_bincontigs-translated.vogout -o ALL_first_bincontigs-translated.vogout
sed -i 's/\t/;/g' ALL_first_bincontigs-translated.vogout

awk 'BEGIN {FS=OFS=";"} {if($4 > $5) print $1,$2"_REV",$3,$5,$4; else print $0}' ALL_first_bincontigs-translated.vogout > tmp && mv tmp ALL_first_bincontigs-translated.vogout

### now we need to add the NCVOG because there can be multiple GVOGs for a single NCVOG
cut -d ";" -f 3 ALL_first_bincontigs-translated.vogout > ALL_first_bincontigs-translated.vogout.tofix.tmp
grep -f ALL_first_bincontigs-translated.vogout.tofix.tmp GVOG_to_NCVOG.tmp | sed -f - ALL_first_bincontigs-translated.vogout.tofix.tmp | paste -d ";" ALL_first_bincontigs-translated.vogout - > tmp && mv tmp ALL_first_bincontigs-translated.vogout

#### NEW part
### identify overlaps between the the gene predictions and framtrans together

awk -F ";" '{print $1"\t"$4"\t"$5"\t"$2"___"$3"___"$6}' ALL_first_bincontigs.vogout > ALL_first_bincontigs.vogout.bed
awk -F ";" '{print $1"\t"$4"\t"$5"\t"$2"___"$3"___"$6}' ALL_first_bincontigs-translated.vogout > ALL_first_bincontigs-translated.vogout.bed

## get the frametrans overlaps with no genepred overlaps
bedtools intersect -v -wa -wb -b ALL_first_bincontigs.vogout.bed -a ALL_first_bincontigs-translated.vogout.bed > ALL_first_bincontigs.vogout.4.tmp

###################

## get the frametrans overlaps with genepred overlaps
bedtools intersect -wa -wb -b ALL_first_bincontigs.vogout.bed -a ALL_first_bincontigs-translated.vogout.bed > ALL_first_bincontigs.vogout.5.tmp
## now filter it to keep only the frametrans that do NOT have the same NCVOG/GVOG as the overlap
## first list the frametrans that DID overlap and have same NCVOG as a genepred, then remove them all from those that did have a different ncvog but overlapped
sed -i 's/___\|\t/;/g' ALL_first_bincontigs.vogout.5.tmp
awk -F ";" '{if($6 == $12) print $1";"$2";"$3";"$4"___"$5"___"$6}' ALL_first_bincontigs.vogout.5.tmp | sort -Vu > ALL_first_bincontigs.vogout.5.tmp.overlap_samencvog.list
awk -F ";" '{if($6 != $12) print $1";"$2";"$3";"$4"___"$5"___"$6}' ALL_first_bincontigs.vogout.5.tmp | sort -Vu > ALL_first_bincontigs.vogout.5.tmp.overlap_diffncvog.list
grep -vwF -f ALL_first_bincontigs.vogout.5.tmp.overlap_samencvog.list ALL_first_bincontigs.vogout.5.tmp.overlap_diffncvog.list | sed 's/;/\t/g' > ALL_first_bincontigs.vogout.5.tmp


### concatenate teh gene preds and the frame trans which overlap genepreds ONLY if they have a different NCVOG
cat ALL_first_bincontigs.vogout.bed ALL_first_bincontigs.vogout.4.tmp ALL_first_bincontigs.vogout.5.tmp | sed 's/___/\t/g' | sort -Vu | sort -t $'\t' -k 1,1V -k 2,3g > ALL_first_bincontigs.vogout.6.tmp

### now we need to flag consecutive hits to the same NCVOG/GVOG as FR, IN, or EX
while IFS=$'\t' read -r contig start stop gene gvog ncvog; do
  if [[ "$contig" == "$prev_contig" ]] && [[ "$ncvog" == "$prev_ncvog" ]]; then
    echo "$contig;$start;$stop;$gene;$gvog;$ncvog;$(($start-$prev_stop))"
  elif [[ "$contig" != "$prev_contig" ]] || [[ "$ncvog" != "$prev_ncvog" ]]; then
    echo "$contig;$start;$stop;$gene;$gvog;$ncvog;"
  fi

  # Store the current line's values as previous for the next iteration
  prev_contig="$contig"
  prev_start="$start"
  prev_stop="$stop"
  prev_ncvog="$ncvog"
  prev_gene="$gene"

done < ALL_first_bincontigs.vogout.6.tmp > ALL_first_bincontigs.vogout.7.tmp

### now reverse and repeat for leading genes

sort -t ';' -k 1,1V -k 2,3gr ALL_first_bincontigs.vogout.7.tmp -o ALL_first_bincontigs.vogout.7.tmp

while IFS=';' read -r contig start stop gene gvog ncvog gap; do
  if [[ ! $gap =~ [0-9] ]] && [[ $prev_gap =~ [0-9] ]]; then
    echo "$contig;$start;$stop;$gene;$gvog;$ncvog;$prev_gap"
  elif [[ $gap =~ [0-9] ]] || [[ ! $prev_gap =~ [0-9] ]]; then
    echo "$contig;$start;$stop;$gene;$gvog;$ncvog;$gap"
  fi

  # Store the current line's values as previous for the next iteration
  prev_contig="$contig"
  prev_start="$start"
  prev_stop="$stop"
  prev_ncvog="$ncvog"
  prev_gene="$gene"
  prev_gap="$gap"


done < ALL_first_bincontigs.vogout.7.tmp > tmp
mv tmp ALL_first_bincontigs.vogout.7.tmp

### now flag the chance of a gene not being monoexonic

awk 'BEGIN {FS=OFS=";"} {if($7 <= 0 && $7 != "") print $0,"FR_y__IN_n__EX_n"; else print $0}' ALL_first_bincontigs.vogout.7.tmp > tmp && mv tmp ALL_first_bincontigs.vogout.7.tmp ## overlapping and same NCVOG - likely frameshift

awk 'BEGIN {FS=OFS=";"} {if($7 >= 1 && $7 <= 100) print $0,"FR_y__IN_m__EX_n"; else print $0}' ALL_first_bincontigs.vogout.7.tmp > tmp && mv tmp ALL_first_bincontigs.vogout.7.tmp ## small intergene gap and same NCVOG - likely frameshift

awk 'BEGIN {FS=OFS=";"} {if($7 >= 101 && $7 <= 5000) print $0,"FR_m__IN_y__EX_n"; else print $0}' ALL_first_bincontigs.vogout.7.tmp > tmp && mv tmp ALL_first_bincontigs.vogout.7.tmp ## intergene gap of range around the average in brown algae or phaeoviruses and same NCVOG - likely intron

awk 'BEGIN {FS=OFS=";"} {if($7 >= 5001 && $7 <= 10000) print $0,"FR_n__IN_m__EX_y"; else print $0}' ALL_first_bincontigs.vogout.7.tmp > tmp && mv tmp ALL_first_bincontigs.vogout.7.tmp ## intergene gap beyond the range around the average in brown algae or phaeoviruses and same NCVOG - likely monoexonic gene - but 1 % of brown algae genes can have introns this large

awk 'BEGIN {FS=OFS=";"} {if($7 >= 10001) print $0,"FR_n__IN_n__EX_y"; else print $0}' ALL_first_bincontigs.vogout.7.tmp > tmp && mv tmp ALL_first_bincontigs.vogout.7.tmp ## intergene gap way beyond the range around the average in brown algae or phaeoviruses and same NCVOG - likely monoexonic gene

##########################################
fi
##########################################


##########################################
if [[ $1 == "part2" ]]; then
##########################################

##################################
####################################
#### NOW we need to flag which genes are maybe part of multiexonic genes
### get mRNA lines from gff for IVEX001 for contigs of interest only and on those with multiple exons

## get all deflines to convert back to original contig names
#rm -f tmp.all.deflines.phex; touch tmp.all.deflines.phex; for f in $(find /home/dcschroe/dmckeow/data/phaeoexplorer_main -name "*.deflinekey"); do cat $f >> tmp.all.deflines.phex; done

### moved to scratch due to storage over quota
rm -f tmp.all.deflines.phex; touch tmp.all.deflines.phex; for f in $(find /scratch.global/dcschroe/dmckeow/BINNING_* -name "*.deflinekey"); do cat $f >> tmp.all.deflines.phex; done

## lsit the contigs that have viralregions
cut -d ";" -f 1 ALL_first_bincontigs.vogout.7.tmp | sort -Vu > tmp.contigs.with.vrs
## reduce to only contigs of interest with orignial names
grep -wF -f tmp.contigs.with.vrs tmp.all.deflines.phex > tmp.contigs.with.vrs.ognames
sed -i -E 's/_part[0-9]+//g' tmp.contigs.with.vrs.ognames

##########################################
fi
##########################################


##########################################
if [[ $1 == "part2_MANUAL" ]]; then
##########################################

### THEN TRANSFER tmp.contigs.with.vrs.ognames TO SB-ROSCOFF
  for f in *.gff.gz; do zcat $f | sed '/\tGeneMark.hmm\t/d' > ${f}.tmp.vr1.vr2 ; done
  sed -E 's/(contig_)(.*)(_[0-9]+\t.+)/\2\t\1\2\3/g' tmp.contigs.with.vrs.ognames | sort -Vu | awk '{print $3 > "tmp.contigs.with.vrs.ognames."$1}'
  for f in tmp.contigs.with.vrs.ognames.*; do sort -Vu $f -o $f ; done

  ls tmp.contigs.with.vrs.ognames.* > fix.list
#### MANUALLY edit fix.list to have two columns, with correct corresponding files e.g. :
#tmp.contigs.with.vrs.ognames.Ascophyllum_nodosum_M      Ascophyllum-nodosum_MALE.gff.gz.tmp.vr1.vr2
#tmp.contigs.with.vrs.ognames.Chordaria_linearis Chordaria-linearis.gff.gz.tmp.vr1.vr2
#tmp.contigs.with.vrs.ognames.Chrysoparadoxa_australica  Chrysoparadoxa-australica.gff.gz.tmp.vr1.vr2
#tmp.contigs.with.vrs.ognames.Desmarestia_dudresnayi     Desmarestia-dudresnayi.gff.gz.tmp.vr1.vr2

### now lets reduce down to the gff lines only for the contigs of interest
## some contigs will not make it through if they dont have any gene predictions
awk -F "\t" '{print "grep -wF -f "$1" "$2" > "$2".vr3"}' fix.list | bash -

## nwo filter out only multiexonic genes THAT are phaeoexplorer
for f in *.gff.gz.tmp.vr1.vr2.vr3; do grep ';exons=[0-9]' $f | sed -E '/;exons=1;|;exons=1$/d' > ${f}.vr4; done

### now for the public genomes 
### also add exon count to the end of each line by counting CDS
############ PUBLIC genomes with CDS annotations that represent exons - the same CDS ID repeated for every CDS per gene

#genome="PUBLIC_Cladosiphon-okamuranus.gff.gz.tmp.vr1.vr2.vr3"
#genome="PUBLIC_Ectocarpus-sp7.gff.gz.tmp.vr1.vr2.vr3"
#genome="PUBLIC_Ectocarpus-subulatus.gff.gz.tmp.vr1.vr2.vr3"
#genome="PUBLIC_Nemacystus-decipiens.gff.gz.tmp.vr1.vr2.vr3"
#genome="PUBLIC_Saccharina-japonica.gff.gz.tmp.vr1.vr2.vr3"
#genome="PUBLIC_Sargassum-fusiforme.gff.gz.tmp.vr1.vr2.vr3"
#genome="PUBLIC_Tribonema-minus.gff.gz.tmp.vr1.vr2.vr3"
#genome="PUBLIC_Undaria-pinnatifida.gff.gz.tmp.vr1.vr2.vr3"

## before running genomes:
#sed -E -i 's/\tID=CDS[0-9]+\./\tID=CDS./g' PUBLIC_Undaria-pinnatifida.gff.gz.tmp.vr1.vr2.vr3
#sed -E -i 's/ |\"/_/g' PUBLIC_Sargassum-fusiforme.gff.gz.tmp.vr1.vr2.vr3
###

### RUN this via an sbatch script on sb-roscoff
genome="$1"
sed -i '/^$/d' $genome
awk '/\tCDS\t/' $genome | awk -F "\t|;" '{print $9}' | sed '/^$/d' | sort -V | uniq -d | grep -wF -f - $genome | awk '/\tCDS\t/' > ${genome}.vr4
rm -f ${genome}.tmp; touch ${genome}.tmp
for f in $(awk -F "\t|;" '{print $9}' ${genome}.vr4); do grep -cwF $f ${genome}.vr4 >> ${genome}.tmp; done
sed 's/^/exons=/g' ${genome}.tmp | paste -d ";" ${genome}.vr4 - > ${genome}.tmp2 && mv ${genome}.tmp2 ${genome}.vr4
sed -i '/;exons=1$/d' ${genome}.vr4

## when the scripts are all done:
for f in *.gff.gz.tmp.vr1.vr2.vr3.vr4; do awk 'BEGIN {FS=OFS="\t"} {print $1,$4,$5,$9}' $f > ${f}.bed; done


### now add back the anvio names for contigs. only for those that are relevent
cut -f 1 *.gff.gz.tmp.vr1.vr2.vr3.vr4.bed | sort -Vu | grep -wF -f - tmp.contigs.with.vrs.ognames | sed -E 's/(contig_)(.*)(_[0-9]+\t.+)/\2\t\1\2\3/g' | sort -Vu | awk '{print "/^"$3"\\t/ s/^"$3"\\t/"$2"\\t/g" > "tmp.contigs.with.vrs.ognames."$1}'
cp fix.list fix.list2 ## use the same genome list as before
sed -i 's/.gff.gz.tmp.vr1.vr2/.gff.gz.tmp.vr1.vr2.vr3.vr4.bed/g' fix.list2
sed -i 's/^/sed -f /g' fix.list2
sed -i -E 's/(.*)\t(.*)/\1 \2 > \2.final/g' fix.list2

### IF NEEDED change back the _part to specific genomes that had to have their contigs split for prodigal
## PUBLIC_Undaria-pinnatifida - LG0 split into 6 parts (max 30000000 bp)
## PUBLIC_Saccharnia-japonica - chr0 split into 7 parts (max 30000000 bp)

awk '/^chr0\t/' PUBLIC_Saccharina-japonica.gff.gz.tmp.vr1.vr2.vr3.vr4.bed > tmp.gff.1
awk -F "\t" '{if($2 >= 1 && $3 <= 30000000) print "contig_P_Saccharina_japonica_000000000001\t"$0}' tmp.gff.1 > tmp.gff.2
awk -F "\t" '{if($2 >= 30000001 && $3 <= 60000000) print "contig_P_Saccharina_japonica_000000000002\t"$0}' tmp.gff.1 >> tmp.gff.2
awk -F "\t" '{if($2 >= 60000001 && $3 <= 90000000) print "contig_P_Saccharina_japonica_000000000003\t"$0}' tmp.gff.1 >> tmp.gff.2
awk -F "\t" '{if($2 >= 90000001 && $3 <= 120000000) print "contig_P_Saccharina_japonica_000000000004\t"$0}' tmp.gff.1 >> tmp.gff.2
awk -F "\t" '{if($2 >= 120000001 && $3 <= 150000000) print "contig_P_Saccharina_japonica_000000000005\t"$0}' tmp.gff.1 >> tmp.gff.2
awk -F "\t" '{if($2 >= 150000001 && $3 <= 180000000) print "contig_P_Saccharina_japonica_000000000006\t"$0;}' tmp.gff.1 >> tmp.gff.2
awk -F "\t" '{if($2 >= 180000001 && $3 <= 210000000) print "contig_P_Saccharina_japonica_000000000006\t"$0;}' tmp.gff.1 >> tmp.gff.2

sed -i '/^chr0\t/d' PUBLIC_Saccharina-japonica.gff.gz.tmp.vr1.vr2.vr3.vr4.bed

awk '/^LG0\t/' PUBLIC_Undaria-pinnatifida.gff.gz.tmp.vr1.vr2.vr3.vr4.bed > tmp.gff.3
awk -F "\t" '{if($2 >= 1 && $3 <= 30000000) print "contig_P_Undaria_pinnatifida_000000000001\t"$0}' tmp.gff.3 > tmp.gff.4
awk -F "\t" '{if($2 >= 30000001 && $3 <= 60000000) print "contig_P_Undaria_pinnatifida_000000000002\t"$0}' tmp.gff.3 >> tmp.gff.4
awk -F "\t" '{if($2 >= 60000001 && $3 <= 90000000) print "contig_P_Undaria_pinnatifida_000000000003\t"$0}' tmp.gff.3 >> tmp.gff.4
awk -F "\t" '{if($2 >= 90000001 && $3 <= 120000000) print "contig_P_Undaria_pinnatifida_000000000004\t"$0}' tmp.gff.3 >> tmp.gff.4
awk -F "\t" '{if($2 >= 120000001 && $3 <= 150000000) print "contig_P_Undaria_pinnatifida_000000000005\t"$0}' tmp.gff.3 >> tmp.gff.4
awk -F "\t" '{if($2 >= 150000001 && $3 <= 180000000) print "contig_P_Undaria_pinnatifida_000000000006\t"$0;}' tmp.gff.3 >> tmp.gff.4

sed -i '/^LG0\t/d' PUBLIC_Undaria-pinnatifida.gff.gz.tmp.vr1.vr2.vr3.vr4.bed

bash fix.list2

cut --complement -f 2 tmp.gff.2 >> PUBLIC_Saccharina-japonica.gff.gz.tmp.vr1.vr2.vr3.vr4.bed.final
cut --complement -f 2 tmp.gff.4 >> PUBLIC_Undaria-pinnatifida.gff.gz.tmp.vr1.vr2.vr3.vr4.bed.final


############ NOW TRANSFER BACK TO MSI UMN
scp -r dmckeown@slurm0.sb-roscoff.fr:/shared/projects/phaeoexplorer_virus/phaeoex_screen/finalresult/IVEX002/*.gff.gz.tmp.vr1.vr2.vr3.vr4.bed.final dmckeow@agate.msi.umn.edu:/home/dcschroe/dmckeow/viralrecall/

##########################################
fi
##########################################

##########################################
if [[ $1 == "part3" ]]; then
##########################################

sed -e 's/;/\t/1' -e 's/;/\t/1' -e 's/;/\t/1' ALL_first_bincontigs.vogout.7.tmp > ALL_first_bincontigs.vogout.7.tmp.bed
sed -i 's/;$/;;/g' ALL_first_bincontigs.vogout.7.tmp.bed

### now check if our viral region genes overlap with multiexonic genes
rm -f ALL_first_bincontigs.vogout.8.tmp.bed; touch ALL_first_bincontigs.vogout.8.tmp.bed

for f in *.gff.gz.tmp.vr1.vr2.vr3.vr4.bed.final; do
  bedtools intersect -wa -wb -b $f -a ALL_first_bincontigs.vogout.7.tmp.bed | sed 's/;exons=/\toverlap_exons=/g' | cut -f 1-4,9 | sort -V | sed -E 's/\t(overlap_exons=[0-9]+).*/;\1/g' >> ALL_first_bincontigs.vogout.8.tmp.bed
done

## remove duplicate hits, keeping highest exon hit
sort -t ';' -k 1,1V -k 6,6Vr ALL_first_bincontigs.vogout.8.tmp.bed | awk -F ";" '!a[$1]++' > tmp && mv tmp ALL_first_bincontigs.vogout.8.tmp.bed

#####sed -E -i 's/(FR_.__)IN_.__EX_./\1IN_y__EX_n/g' ALL_first_bincontigs.vogout.8.tmp.bed
sed -i 's/;;overlap_exons=/;FR_n__IN_y__EX_n;overlap_exons=/g' ALL_first_bincontigs.vogout.8.tmp.bed


rm -f ALL_first_bincontigs.vogout.9.tmp.bed; touch ALL_first_bincontigs.vogout.9.tmp.bed

## get the genes with no overlaps to multiexon genes
for f in *.gff.gz.tmp.vr1.vr2.vr3.vr4.bed.final; do
  bedtools intersect -v -wa -wb -b $f -a ALL_first_bincontigs.vogout.7.tmp.bed | sort -Vu | sed 's/;$/;FR_n__IN_n__EX_y/g' >> ALL_first_bincontigs.vogout.9.tmp.bed
done
sort -Vu ALL_first_bincontigs.vogout.9.tmp.bed > tmp && mv tmp ALL_first_bincontigs.vogout.9.tmp.bed
sed -i 's/$/;overlap_exons=NONE/g' ALL_first_bincontigs.vogout.9.tmp.bed

## remove the false negative hits
cut -d ";" -f 1 ALL_first_bincontigs.vogout.8.tmp.bed | grep -vwF -f - ALL_first_bincontigs.vogout.9.tmp.bed > tmp
mv tmp ALL_first_bincontigs.vogout.9.tmp.bed

## final merge
cat ALL_first_bincontigs.vogout.8.tmp.bed ALL_first_bincontigs.vogout.9.tmp.bed | sort -t $'\t' -k 1,1V -k 2,3g > ALL_first_bincontigs.vogout.10.tmp.bed

sed -i 's/\t/;/g' ALL_first_bincontigs.vogout.10.tmp.bed

#################################
## NOW get marker gene annotations
#################################

### get ALL marker gene hmm hits
## NOTE that the marker gene is matched to gene name, so some VOG codes won't correspond to the marker gene flag
### from prodigal gene predictions
rm -f ALL_first_bincontigs.markerout; touch ALL_first_bincontigs.markerout; for f in $(find ViralRecall_ALL_first_bincontigs2.fa.split.*/ALL_first_bincontigs2.part_* -name "*.markerout" | sed '/-translated.markerout/d'); do cat $f >> ALL_first_bincontigs.markerout ; done

sed -E -e '/^#/d' -e 's/ # |;/\t/g' -e 's/ +/\t/g' ALL_first_bincontigs.markerout | cut -f 1,3,19,20 | sed -E 's/(.+[0-9]+)_([0-9]+\t)(.+)/\1\t\2\3/1' > tmp && mv tmp ALL_first_bincontigs.markerout

### from frametrans
rm -f ALL_first_bincontigs-translated.markerout; touch ALL_first_bincontigs-translated.markerout; for f in $(find ALL_first_bincontigs2.fa.split.*-translated -name "*-translated.markerout"); do cat $f >> ALL_first_bincontigs-translated.markerout ; done

sed -E -e '/^#/d' -e 's/ [a-z]+=|\.\./\t/g' -e 's/ +/\t/g' ALL_first_bincontigs-translated.markerout | awk 'BEGIN {FS=OFS="\t"} {print $19,$1,$3,$20,$21}' > tmp && mv tmp ALL_first_bincontigs-translated.markerout

## ALL_first_bincontigs.vogout.9.tmp to ALL_first_bincontigs.vogout.10.tmp.bed

awk -F ";" '{print $1"_"$4"___"}' ALL_first_bincontigs.vogout.10.tmp.bed | sed 's/_REV___/___/g' > marker_tofix.tmp

cat ALL_first_bincontigs.markerout ALL_first_bincontigs-translated.markerout | awk -F "\t" '{print "s/"$1"_"$2"___/"$3"/g"}' | grep -wF -f marker_tofix.tmp - > marker.fix.tmp

sed -f marker.fix.tmp marker_tofix.tmp > marker.fixed.tmp
sed -i 's/.*___$/na/g' marker.fixed.tmp

paste -d ";" ALL_first_bincontigs.vogout.10.tmp.bed marker.fixed.tmp > ALL_first_bincontigs.vogout.11.tmp.bed

## ALL_first_bincontigs.vogout.10.tmp to ALL_first_bincontigs.vogout.11.tmp.bed

### now we just need to connect the ALL_first_bincontigs.vogout.FINAL with the viral regions summary
## ALL_summary.tsv vs ALL_first_bincontigs.vogout.10.tmp

awk -F ";" '{print $2"\t"$3"\t"$4"\t"$1}' ALL_summary.tsv > tmp.ALL_summary.bed
awk -F ";" '{print $1"\t"$2"\t"$3}' ALL_first_bincontigs.vogout.11.tmp.bed > ALL_first_bincontigs.vogout.12.tmp.bed
bedtools intersect -wa -wb -a tmp.ALL_summary.bed -b ALL_first_bincontigs.vogout.12.tmp.bed > tmp.ALL_summary.vroverlaps

cut -d ";" -f 1-3 ALL_first_bincontigs.vogout.11.tmp.bed | sed 's/$/;/g' > tmp.tofix.vr

awk -F "\t" '{print "/^"$5";"$6";"$7";$/ s/"$5";"$6";"$7";/"$4"/g"}' tmp.ALL_summary.vroverlaps > tmp.fixer.vr

###sbatch --time=2:00:00 --cpus-per-task=4 --mem=32GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='sed -f tmp.fixer.vr tmp.tofix.vr > tmp.fixed.vr'
sed -f tmp.fixer.vr tmp.tofix.vr > tmp.fixed.vr

## add the viralregion names
awk -F "\t" '{if($0 !~ "viral_region_") print "na"; else print $0}' tmp.fixed.vr | paste -d ";" ALL_first_bincontigs.vogout.11.tmp.bed - > ALL_phexmain_viralrecall_processed.summary.genes
## add genome name
cut -d ";" -f 1 ALL_first_bincontigs.vogout.11.tmp.bed | sed -E 's/contig_(.+)_[0-9]+$/\1/g' | paste -d ";" - ALL_phexmain_viralrecall_processed.summary.genes > tmp && mv tmp ALL_phexmain_viralrecall_processed.summary.genes

### add viralregion context codes
awk -F ";" '{print "/^"$1";$/ s/"$1";/"$2"/g"}' tmp.ALL.viralregions.by_contig > tmp.ALL.viralregions.by_contig.fixer
cut -d ";" -f 2 ALL_phexmain_viralrecall_processed.summary.genes | sed 's/$/;/g' > tmp.ALL.viralregions.by_contig.tofix

##sbatch --time=2:00:00 --cpus-per-task=4 --mem=32GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='sed -f tmp.ALL.viralregions.by_contig.fixer tmp.ALL.viralregions.by_contig.tofix > tmp.ALL.viralregions.by_contig.fixed'
sed -f tmp.ALL.viralregions.by_contig.fixer tmp.ALL.viralregions.by_contig.tofix > tmp.ALL.viralregions.by_contig.fixed

awk -F "\t" '{if($0 ~ ";") print "na"; else print $0}' tmp.ALL.viralregions.by_contig.fixed | paste -d ";" ALL_phexmain_viralrecall_processed.summary.genes - > tmp && mv tmp ALL_phexmain_viralrecall_processed.summary.genes

### same for all summary info by viralregion
awk -F ";" '{print "/^"$2";"$1";$/ s/"$2";"$1";/"$3";"$4";"$5";"$6";"$7";"$8";"$9"/g"}' ALL_summary.tsv > ALL_summary.tmp.fixer
cut -d ";" -f 2,12 ALL_phexmain_viralrecall_processed.summary.genes | sed 's/$/;/g' > ALL_summary.tmp.tofix


###sbatch --time=2:00:00 --cpus-per-task=4 --mem=32GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='sed -f ALL_summary.tmp.fixer ALL_summary.tmp.tofix > ALL_summary.tmp.fixed'
sed -f ALL_summary.tmp.fixer ALL_summary.tmp.tofix > ALL_summary.tmp.fixed

awk -F "\t" '{if($0 ~ /^contig.*;$/) print "na;na;na;na;na;na;na"; else print $0}' ALL_summary.tmp.fixed | paste -d ";" ALL_phexmain_viralrecall_processed.summary.genes - > tmp && mv tmp ALL_phexmain_viralrecall_processed.summary.genes

cut --complement -d ";" -f 8 ALL_phexmain_viralrecall_processed.summary.genes > tmp && mv tmp ALL_phexmain_viralrecall_processed.summary.genes

sed -i 's/overlap_exons=//g' ALL_phexmain_viralrecall_processed.summary.genes

sed -z -i 's/^/genome__id;contig__id;gene__start;gene__stop;gene__id;gene__gvog;gene__gvog_ncvog;gene__frame_intron_monoex_flag;gene__num_exons_of_overlapping_multiexonic_gene;gene__ncldv_marker;viralregion__name;contig__viralregion_context;viralregion__start;viralregion__stop;viralregion__length;contig__length;viralregion__score;viralregion__num_viralhits;viralregion__num_orfs\n/g' ALL_phexmain_viralrecall_processed.summary.genes


### get all viralrecall hits and scores
rm -f test1; touch test1
for f in $(find ViralRecall_ALL_first_bincontigs2.fa.split.*/ALL_first_bincontigs2.part_* -name "*.vogout"); do cat $f >> test1 ; done

for f in $(find ALL_first_bincontigs2.fa.split.*-translated -name "*-translated.vogout"); do cat $f >> test1 ; done

sed -E -e '/^#/d' -e 's/ # |;/\t/g' -e 's/ +/\t/g' test1 | cut -f 1,3,5,6,19,20 | sed -E -e 's/source=|coords=//g' -e 's/\.\./\t/g' | awk 'BEGIN {FS=OFS="\t"} {if($0 ~ /orf[0-9]+\t/) print $5"_"$1,$2,$3,$4,$6,$7; else print $0}' | sed 's/\t/;/g' > tmp.gene.scores

awk 'BEGIN {FS=OFS=";"} {print "/^"$1,$2"$/ s/"$1,$2"/"$3,$4,$6-$5"/g"}' tmp.gene.scores | sed -E 's/;-([0-9]+\/g)$/;\1/g' > tmp.genes.fixer

cut -d ";" -f 2,5,6 ALL_phexmain_viralrecall_processed.summary.genes | sed 's/;/_/1' | sed '1,1d' | sed 's/_REV;/;/g' > tmp.genes.tofix

split -d -n l/100 tmp.genes.tofix tmp.split.genes.tofix.

for f in tmp.split.genes.tofix.*; do
  sbatch --time=2:00:00 --cpus-per-task=2 --mem=64GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap="grep -wF -f $f tmp.genes.fixer | sed -f - $f > fixed.${f}"
done

#####################
fi
#####################

##########################################
if [[ $1 == "part4" ]]; then
##########################################

cat fixed.tmp.split.genes.tofix.* | sed -z 's/^/gene__hmm_evalue;gene__hmm_score;gene__length\n/g' | paste -d ";" ALL_phexmain_viralrecall_processed.summary.genes - > tmp && mv tmp ALL_phexmain_viralrecall_processed.summary.genes

######### NEED to check for overlapping genes of the same GVOG/NCVOG between frametrans and separately for between genepreds

## for frametrans
bedtools intersect -wao -b ALL_first_bincontigs-translated.vogout.bed -a ALL_first_bincontigs-translated.vogout.bed | sed -E 's/\t|___/;/g' | awk 'BEGIN {FS=OFS=";"} {if($6 == $12 && $4 != $10) print $0}' > ALL_first_bincontigs.vogout.12.tmp

## for genepreds
bedtools intersect -wao -b ALL_first_bincontigs.vogout.bed -a ALL_first_bincontigs.vogout.bed | sed -E 's/\t|___/;/g' | awk 'BEGIN {FS=OFS=";"} {if($6 == $12 && $4 != $10) print $0}' > ALL_first_bincontigs.vogout.13.tmp

## merge
cat ALL_first_bincontigs.vogout.12.tmp ALL_first_bincontigs.vogout.13.tmp > ALL_first_bincontigs.vogout.14.tmp

### change frameshift flag to yes if frametrans OR gene preds have overlaps with the same NCVOG/GVOG
cut -d ";" -f 1-4 ALL_first_bincontigs.vogout.14.tmp | sed 's/$/;/g' | awk 'BEGIN {FS=OFS=";"} {print "/"$0"/ s/FR_.__IN/FR_y__IN/g"}' | sed -f - ALL_phexmain_viralrecall_processed.summary.genes > tmp

### change frameshift flag to yes if frametrans OR gene preds have overlaps with the same NCVOG/GVOG
cut -d ";" -f 1-4 ALL_first_bincontigs.vogout.14.tmp | sed 's/$/;/g' | awk 'BEGIN {FS=OFS=";"} {print "/"$0"/ s/__EX_y/__EX_n/g"}' | sed -f - tmp > tmp2
mv tmp2 ALL_phexmain_viralrecall_processed.summary.genes

#####################
fi
#####################

##########################################
if [[ $1 == "part5" ]]; then
##########################################

########### PFAM hmm hits
### get all viralrecall hits and scores
rm -f tmp.pfamout; touch tmp.pfamout

for f in $(find ViralRecall_ALL_first_bincontigs2.fa.split.*/ALL_first_bincontigs2.part_* -name "*.pfamout"); do cat $f >> tmp.pfamout ; done

for f in $(find ALL_first_bincontigs2.fa.split.*-translated -name "*-translated.pfamout"); do cat $f >> tmp.pfamout ; done

sed -E -e '/^#/d' -e 's/ # |;/\t/g' -e 's/ +/\t/g' tmp.pfamout | cut -f 1,3,5,6,19 | sed -E -e 's/source=//g' -e 's/\.\./\t/g' | awk 'BEGIN {FS=OFS="\t"} {if($0 ~ /orf[0-9]+\t/) print $5"_"$1,$2,$3,$4; else print $1,$2,$3,$4}' | sed 's/\t/;/g' > tmp.gene.scores

awk 'BEGIN {FS=OFS=";"} {print "/^"$1";$/ s/"$1";/"$2,$3,$4"/g"}' tmp.gene.scores > tmp.genes.fixer

cut -d ";" -f 2,5 ALL_phexmain_viralrecall_processed.summary.genes | sed 's/;/_/1' | sed 's/$/;/g' | sed '1,1d' > tmp.genes.tofix

split -d -n l/100 tmp.genes.tofix tmp.split.genes.tofix.

for f in tmp.split.genes.tofix.*; do
  sbatch --time=2:00:00 --cpus-per-task=2 --mem=64GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap="grep -wF -f $f tmp.genes.fixer | sed -f - $f > fixed.pfamout.${f}"
done

#####################
fi
#####################

##########################################
if [[ $1 == "part6" ]]; then
##########################################

cat fixed.pfamout.tmp.split.genes.tofix.* | sed 's/contig.*;$/na;na;na/g' | sed -z 's/^/gene__pfam_hmm_hit_name;gene__pfam_hmm_evalue;gene__pfam_hmm_score\n/g' | paste -d ";" ALL_phexmain_viralrecall_processed.summary.genes - > tmp && mv tmp ALL_phexmain_viralrecall_processed.summary.genes

rm -f fixed.pfamout.tmp.split.genes.tofix.*
rm -f tmp.split.genes.tofix.*


#####################
fi
#####################

##### marker gene AND phylogeny curation (manual):

## gene__frame_intron_monoex_flag
  # frameshifts are ok, but must be noted and counting them all must be avoided
  # genes with introns should not be included
## gene__multi_exon_overlap_flag
  # genes with introns should not be included
## gene__ncldv_marker
  # viralrecall having given a marker gene marker is not definitive
## viralregion__name
  # being within a viralregion increases marker gene probability, but is not definitive
## viralregion__length
  # larger viral region can indicate increased chance of being marker gene, but is not definitive
## viralregion__score
  # higher scoring viral region can indicate increased chance of being marker gene, but is not definitive
## gene__hmm_evalue
  # 
## gene__hmm_score
## gene__length


### 1. genes that were probably multi-exonic (including overlapping with multi-exonic genes) were excluded from being marker genes, with the exception of some genes with 2 or 3 exons, especially if they were within viral regions AND had good viral hits

### 2. genes with monoexonic annotations BUT which fell outside viral regions AND/OR had strong cellular hits, AND/OR low viral hits

### 3. Went over the genes again, considering their cellualr hmm hits results, and removing consecutive hits that were likely part of a single frameshifted gene
  # A32 - cellular hmm hit to Pox_A32 was typical of real viral hits. Many cellular domains hit against A32. DIfficult to separate cellular and viral hits
  # D5 - a typical phaeovirus genome possesses multiple genes that will have strong hits to this VOG - D5_N, DNA_primase_S, PPL4, Herpes_UL52 - we only want to count D5_N as D5, and maybe some with no clear domain hit from pfam. The presence of multiple viral fragments will make it uncertain how copies are present because an ambiguous pfam hit could be the D5 or the viral homolog from another virus, etc.
    # Some genes considered as marker genes even if they may be part of of a gene with a few introns (2 or 3?)

  # intres - diffcult, is not widely conserved in viruses and it is not clear that phaeoviruses even share the same integrase
    # these integrases clearly overlap with a host gene family - so reject multiexons as marker genes
    # there is also no reason to expect this to be a single copy gene
    # count as marker gene if in viral region and/or hits to PFAM Resolvase or MerR


  ## mcp - easily distinguished as it has few cellualr homologs - considering overlapping with 2 or 3 exon genes where the mcp has good hits, because we cannot trust misannotation
    # mcp alone might be the only reliable indicator of virus genome numbers present in a genome

  ## mRNAc - it is very unlikely that phaeoviruses have any RNA polymerase, so expect these all to be cellular homologs
    # the ectocarpus sp7 virus has no real hit to this gene - there is a single gene found in all 3 reference phaeoviruses, but it is a hypothetical protein with no conserved domains recognised - consider all hits to this VOG as false
    # no further work done on this gene, as it has a complex history tangled up with multiple gene families

  ## PolB - quite straightforward, quite a few with introns - count them anyawy

  ## RNAPL - same as before - phaeoviruses unlikely to have RNA polymerase
    # Pelvetia viralregion has one with top viral hits
  ## RNAPS - same as before - phaeoviruses unlikely to have RNA polymerase

  ## RNR - a bit problematic, because the hsot has RNR - but those with exons or good cellular pfam scores are a good indication

  ## SFII - ResIII or DEAD (EsV-1-66) are this VOG but Helicase_C (EsV-1-23) is not
    # these two families are related and share domains but are in different vogs

  ## VLTF3 - quite a virus-specific gene, quite easily identified as viral