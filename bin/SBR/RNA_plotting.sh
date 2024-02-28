#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 20
#SBATCH --mem 128GB

###### MANUAL STEPS AND PREREQUISITES ######
## this is to visualise plot RNA expression data (TPM) along a genome assembly using R
## the RNA expression can also be viewed interactively in Integrated Genome Viewer (.gct file)

####################### INPUTS #####################################
###### load software ######
module load seqkit/0.14.0; module load blast/2.9.0

##### command line arguments required #####
## $1 = IVEX002 gff.gz of interest
## $2 = a file with protein name and TPM for each condition (average if replicated), incl header e.g.:
## $3 = a file with protein name and NumReads for each condition (average if replicated), incl header e.g.:

##### example code to average TPM of 3 replicates:
## paste /shared/projects/phaeoexplorer_virus/counts/3-pseudomapping/Porterinema-fluviatile__SEAWATER_RNAseq_Replicat*/quant.sf | cut -f 1,4,9,14 | sed '1,1d' | awk '{print $1"\t"($2+$3+$4)/3}' > Porterinema-fluviatile__SEAWATER_RNAseq_AVERAGE
## merge all condition TPMs together as in example:
##NAME	CONDITION1_TPM	CONDITION2_TPM	CONDITION3_TPM
##mRNA_P-fluviatile_contig1.1.1	2.51772	4.22896	4.04597
##mRNA_P-fluviatile_contig1.10.1	-1.68467	2.36635	3.58555
##mRNA_P-fluviatile_contig1.100.1	-2.0128	-0.024116	0.155218
##mRNA_P-fluviatile_contig1.101.1	1.34185	2.58805	3.12893
##mRNA_P-fluviatile_contig1.102.1	2.17676	4.62274	4.20783
##mRNA_P-fluviatile_contig1.103.1	4.16524	5.66958	5.24017

###### do the same but for NumReads (number of reads) per gene
## paste /shared/projects/phaeoexplorer_virus/counts/3-pseudomapping/Porterinema-fluviatile__SEAWATER_RNAseq_Replicat*/quant.sf | cut -f 1,5,10,15 | sed '1,1d' | awk '{print $1"\t"($2+$3+$4)/3}' > Porterinema-fluviatile__SEAWATER_RNAseq_AVERAGE_numreads


##### variables #####
n1=$(basename "$1" | cut -d "." -f 1)
n2=$(basename "$2" | cut -d "." -f 1)
n3=$(basename "$3" | cut -d "." -f 1)
fa="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/GeneMarkS/"
FA="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/GeneMarkS/IVEX000_fixed_genomes/*.fa.gz"

####################### SCRIPT #####################################

################ STEP 001 ##############
### assign TPM value to corresponding gff info
zcat $1 | awk -F "\t|;" '{if($3 ~ "mRNA") print $0}' > tmp0
awk -F "\t|;" '{print "___"$9"___"}' tmp0 | sed 's/ID=//g' > tmp1
sed '1d' $2 | awk -F "\t" 'BEGIN{OFS="\t"};{gsub(/\./,"\\.",$1)}1' | sed -e 's/\t/;/1' -e 's/\t/\\t/g' -Ee 's/^(.+);(.+)/s\/___\1___\/\2\/g/g' > tmp2
## remove gene names with no expression or that have expression but are not in gff
sed -i '/0\\t0\\t0/d' tmp2
sed 's/\./\\./g' tmp1 | grep -wF -f - tmp2 > tmp && mv tmp tmp2

### parallel task to replace all appropriate entries quickly, then reunite results
#split -d --number=l/20 tmp1 tmp1_
#for f in tmp1_*; do
  #sed -f tmp2 $f > tmp3_"$f" &
#done
#wait

cat tmp3_tmp1_* | sed -Ee 's/___.*___|^$/0\t0\t0/g' -e 's/\t/::/g' > tmp3

### extract only relevant gff fields and convert to GCT format
paste tmp3 tmp0 | awk -F"\t|;" '{print $10"__"$4"__"$23"__"$22"\tna|@"$2":"$5"-"$6"|\t"log($1+1)/log(2)}' | sed -E 's/ID=|ncvog_abbrev=|ncvog_group=//g' | sed 's/::/\t/g' > tmp4

#### prepare header
c1=$(wc -l tmp4 | awk '{print $1-1}') ## count no. of genes with data
c2=$(head -1 $2 | sed 's/\t/\n/g' | wc -l | awk '{print $1-1}') ## count no. of sample conditions
head -1 $2 | cut --complement -f 1 | sed -E 's/(.+)/Name\tDescription\t\1/1' | sed -E "s/(.+)/#1.2\n$c1\t$c2\n\1/g" > tmp5

cat tmp5 tmp4 > "$n2".gct

### prepare for R heatmap "n"$5,"n"$6
paste tmp3 tmp0 | awk -F"\t|;" 'BEGIN{OFS="\t"};{print $10,$4,$23,$22,$2,$5,$6,$11,$14,$25,$26,$1}' | sed -Ee 's/ID=|ncvog_abbrev=|ncvog_group=|vsalltitles=|vevalue=//g' -Ee 's/salltitles=|evalue=//g' | sed 's/::/\t/g' > "$n2"_heatmap

### manually reformat so that sample is labelled on duplicate line of each gene, rather than a column of numbers e.g.:
### prepare header
sed '1,2d' tmp5 | cut -f 1 | sed 's/.*/NAME\tSOURCE\tABBREV\tGROUP\tCONTIG\tSTART\tEND\tVSALLTITLES\tVEVALUE\tSALLTITLES\tEVALUE\tTPM\tCONDITION/g' > tmp6
## change to have correct conditon name under CONDITION (in order left to right as in the headers of your input file)
cut -f 1-11,12 "$n2"_heatmap | sed 's/$/\tTPM_AVERAGE_FRESHWATER/g' > tmp_tpm_1
cut -f 1-11,13 "$n2"_heatmap | sed 's/$/\tTPM_AVERAGE_SEAWATER/g' > tmp_tpm_2
cat tmp6 tmp_tpm_1 tmp_tpm_2 > tmp && mv tmp "$n2"_heatmap
#### simplify condition names to fit figure header better
sed -i -e 's/\tTPM_AVERAGE_FRESHWATER$/\tFW/g' -e 's/\tTPM_AVERAGE_SEAWATER$/\tSW/g' "$n2"_heatmap

#### log2(TPM+1) transform the TPMs - use this file for plots and extracting contig specific plot info
awk -F "\t" 'BEGIN{OFS="\t"};{gsub(/.*/,log($12+1)/log(2),$12)}1' "$n2"_heatmap | sed 's/EVALUE\t0/EVALUE\tTPM/g' > "$n2"_LOG2p1_heatmap

######## blast expressed genes against genome to get their putative repeats within genome
#sed '1,1d' "$n2"_LOG2p1_heatmap | awk -F"\t|;" 'BEGIN{OFS="\t"};{if($12 > 0) print $5,"SOURCE","MRNA",$6,$7,".",".",".",$1}' | sort -Vu > tmp.gtf
#seqkit subseq --gtf tmp.gtf $FA > tmp.fa

##### self blast
## make database of expressed genome mRNAs
#rm -f tmp.n*
#makeblastdb -in tmp.fa -dbtype nucl -title tmp -out tmp

#### run self blasts
blastn -query tmp.fa -db tmp -out tmp.blast -outfmt "6 qseqid sseqid" -task dc-megablast -evalue 1e-40 -max_target_seqs 500 -num_threads 20 -dust no -soft_masking no

#### prep blast self hits excluding self hits
### here we want to count the number of other mRNAs each mRNA matches (max is -max_target_seqs N)
sed -E 's/:|:\.//g' tmp.blast | sort -Vu | awk -F"\t" '{if($1 != $2) print $1}' | sort -V | uniq -c | sed -Ee 's/ +/\t/g' -Ee 's/^\t//g' | awk -F"\t" '{print "/"$2"/ s/"$2"/"$1"/g"}' > tmp7

#### get lines for non-expressed genes
awk -F"\t" '{print $5"_"$6"-"$7}' "$n2"_LOG2p1_heatmap | sed '1,1d' > tmp8

## parallel jobs
split -d --number=l/20 tmp8 tmp8_
for f in tmp8_*; do
  sed -f tmp7 $f > tmp9_"$f" &
done
wait

cat tmp9_* | sed -E 's/.*[A-Za-z]+.*/0/g' | sed -z 's/^/REPEATS\n/1' | paste "$n2"_LOG2p1_heatmap - > tmp && mv tmp "$n2"_LOG2p1_heatmap

### summarise min, average, max TPMs of whole genome
cut -f 12 "$n2"_LOG2p1_heatmap | sed '1,1d' | sort -nu | head -1 > tmp_min
count=$(cut -f 12 "$n2"_LOG2p1_heatmap | sed '1,1d' | wc -l); cut -f 12 "$n2" | sed '1,1d' | awk -v count=$count '{sum+=$1;} END{print sum/count;}' > tmp_ave
cut -f 12 "$n2"_LOG2p1_heatmap | sed '1,1d' | sort -nur | head -1 > tmp_max
paste tmp_min tmp_ave tmp_max | sed -z 's/^/TPM_min\tTPM_average\tTPM_max\n/g' > "$n2"_LOG2p1_min_ave_max

##### MANUAL AREA - CHANGE TO MATCH YOUR FILES - desired contigs and regions of EVEs
awk '/NAME\t|P-fluviatile_contig7\t/' "$n2"_LOG2p1_heatmap > "$n2"_LOG2p1_contig7_heatmap
awk '/NAME\t|P-fluviatile_contig12\t/' "$n2"_LOG2p1_heatmap > "$n2"_LOG2p1_contig12_heatmap
awk '/NAME\t|P-fluviatile_contig21\t/' "$n2"_LOG2p1_heatmap > "$n2"_LOG2p1_contig21_heatmap
awk '/NAME\t|P-fluviatile_contig22\t/' "$n2"_LOG2p1_heatmap > "$n2"_LOG2p1_contig22_heatmap
awk '/NAME\t|P-fluviatile_contig40\t/' "$n2"_LOG2p1_heatmap > "$n2"_LOG2p1_contig40_heatmap
awk '/NAME\t|P-fluviatile_contig94\t/' "$n2"_LOG2p1_heatmap > "$n2"_LOG2p1_contig94_heatmap

############ GET INTERPROSCAN INFO
#### get fastas of expressed genes for interproscan
for f in *contig*_heatmap; do ff=$(basename "$f" | awk -F"__" '{print $1}'); cut -f 1 $f | sort -Vu | seqkit grep -n -f - "$fa"*"$ff".faa | sed 's/\*//g' > $(basename $f).faa; done

#### run default interproscan on protein fasta via usegalaxy.eu and then do on tab output:
for f in *contig*_heatmap; do awk -F"\t" '/\tPfam\t/ {print $1,$5,$6,$9} {OFS="\t"}' $(basename $f .faa).interproscan | awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/[^a-zA-Z0-9]/,"_",$3)}1' | sed -Ee 's/_+/_/g' -Ee 's/^_|_$//g' -Ee 's/\t_|_\t/\t/g'| sort -t $'\t' -k 1,1 -k 4,4g | awk -F"\t" '!a[$1]++' | awk -F"\t" '{print "s/___"$1"___/"$2"__"$3"/g"}' > $(basename $f .faa).interproscan_reduced; done

for f in *contig*_heatmap; do cut -f 1 $f | sed -e '1,1d' -Ee 's/(.+)/___\1___/g' | sed -f $(basename $f .faa).interproscan_reduced - | sed -E 's/___.+___/-/g' | sed -z 's/^/PFAM\n/1' | paste $f - > tmp && mv tmp $f; done

##### add in number of reads info for contigs
for f in *contig*_heatmap; do cut -f 1 $f | sed '1,1d' | grep -wF -f - $3 | awk -F"\t" '{print "s/"$1"/"$2"\\t"$3"\\t"$4"/g"}' > tmp_numr_1; cut -f 1 $f | sed '1,1d' | sed -f tmp_numr_1 - | sed -E 's/.*[Aa-Zz]+.*/0\t0\t0/g' | sed -z 's/^/NUMREADS\tNUMREADS\tNUMREADS\n/1' | paste $f - | awk -F"\t" 'BEGIN{OFS="\t"};{if($13 == "FW") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$17; else print $0}' | awk -F"\t" 'BEGIN{OFS="\t"};{if($13 == "SW") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$18; else print $0}' | cut -f 1-16 > tmp && mv tmp $f; done

##### extract viral region from contigs of interest to make heatmap of gene identities vs TPM, e.g.:

a="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig7_heatmap"; b="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_EVE7_heatmap"; awk -F "\t" '{if($6 >= 3498000 && $7 <= 3815000) print $0}' $a > $b; cut -f 6 $b | sed 's/$/;/g' | tr -d '\n' | sed -e 's/;$/\n/g' -Ee 's/([0-9]+);.+;([0-9]+)/\1\n\2/g' > $(basename "$b" _heatmap)_coords; grep -nwF -f $(basename "$b" _heatmap)_coords $a | cut -d ":" -f 1 | head -2 | paste $(basename "$b" _heatmap)_coords - | sed -z 's/^/START_END_LOCI\tSTART_END_LINES\n/1' > tmp && mv tmp $(basename "$b" _heatmap)_coords

a="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig12_heatmap"; b="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_EVE12_heatmap"; awk -F "\t" '{if($6 >= 22500 && $7 <= 345000) print $0}' $a > $b; cut -f 6 $b | sed 's/$/;/g' | tr -d '\n' | sed -e 's/;$/\n/g' -Ee 's/([0-9]+);.+;([0-9]+)/\1\n\2/g' > $(basename "$b" _heatmap)_coords; grep -nwF -f $(basename "$b" _heatmap)_coords $a | cut -d ":" -f 1 | head -2 | paste $(basename "$b" _heatmap)_coords - | sed -z 's/^/START_END_LOCI\tSTART_END_LINES\n/1' > tmp && mv tmp $(basename "$b" _heatmap)_coords

a="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig21_heatmap"; b="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_EVE21_heatmap"; awk -F "\t" '{if($6 >= 297000 && $7 <= 414000) print $0}' $a > $b; cut -f 6 $b | sed 's/$/;/g' | tr -d '\n' | sed -e 's/;$/\n/g' -Ee 's/([0-9]+);.+;([0-9]+)/\1\n\2/g' > $(basename "$b" _heatmap)_coords; grep -nwF -f $(basename "$b" _heatmap)_coords $a | cut -d ":" -f 1 | head -2 | paste $(basename "$b" _heatmap)_coords - | sed -z 's/^/START_END_LOCI\tSTART_END_LINES\n/1' > tmp && mv tmp $(basename "$b" _heatmap)_coords

a="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig22_heatmap"; b="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_EVE22_heatmap"; awk -F "\t" '{if($6 >= 1171000 && $7 <= 1454000) print $0}' $a > $b; cut -f 6 $b | sed 's/$/;/g' | tr -d '\n' | sed -e 's/;$/\n/g' -Ee 's/([0-9]+);.+;([0-9]+)/\1\n\2/g' > $(basename "$b" _heatmap)_coords; grep -nwF -f $(basename "$b" _heatmap)_coords $a | cut -d ":" -f 1 | head -2 | paste $(basename "$b" _heatmap)_coords - | sed -z 's/^/START_END_LOCI\tSTART_END_LINES\n/1' > tmp && mv tmp $(basename "$b" _heatmap)_coords

a="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig40_heatmap"; b="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_EVE40-1_heatmap"; awk -F "\t" '{if($6 >= 1557000 && $7 <= 1695800) print $0}' $a > $b; cut -f 6 $b | sed 's/$/;/g' | tr -d '\n' | sed -e 's/;$/\n/g' -Ee 's/([0-9]+);.+;([0-9]+)/\1\n\2/g' > $(basename "$b" _heatmap)_coords; grep -nwF -f $(basename "$b" _heatmap)_coords $a | cut -d ":" -f 1 | head -2 | paste $(basename "$b" _heatmap)_coords - | sed -z 's/^/START_END_LOCI\tSTART_END_LINES\n/1' > tmp && mv tmp $(basename "$b" _heatmap)_coords

a="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig40_heatmap"; b="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_EVE40-2_heatmap"; awk -F "\t" '{if($6 >= 220000 && $7 <= 248000) print $0}' $a > $b; cut -f 6 $b | sed 's/$/;/g' | tr -d '\n' | sed -e 's/;$/\n/g' -Ee 's/([0-9]+);.+;([0-9]+)/\1\n\2/g' > $(basename "$b" _heatmap)_coords; grep -nwF -f $(basename "$b" _heatmap)_coords $a | cut -d ":" -f 1 | head -2 | paste $(basename "$b" _heatmap)_coords - | sed -z 's/^/START_END_LOCI\tSTART_END_LINES\n/1' > tmp && mv tmp $(basename "$b" _heatmap)_coords

a="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_contig94_heatmap"; b="Porterinema-fluviatile__RNAseq_TPM_LOG2p1_EVE94_heatmap"; awk -F "\t" '{if($6 >= 43000 && $7 <= 85000) print $0}' $a > $b; cut -f 6 $b | sed 's/$/;/g' | tr -d '\n' | sed -e 's/;$/\n/g' -Ee 's/([0-9]+);.+;([0-9]+)/\1\n\2/g' > $(basename "$b" _heatmap)_coords; grep -nwF -f $(basename "$b" _heatmap)_coords $a | cut -d ":" -f 1 | head -2 | paste $(basename "$b" _heatmap)_coords - | sed -z 's/^/START_END_LOCI\tSTART_END_LINES\n/1' > tmp && mv tmp $(basename "$b" _heatmap)_coords

#### get info of all potentially new NON-NCVOG proteins with putative functions
for f in *EVE*_heatmap; do awk -F "\t" '{if($3 == "-" && $14 != "-") print $0}' $f | cut -f 1-11,14,15 | sort -Vu; done > "$n2"_nNCVOG_PFAM

#### reduce to only expressed genes in viral regions and limit length of hit names, and sort by TPM
for f in *EVE*_heatmap; do awk -F "\t" '{if($12 > 0) print $0}' $f | cut -f 1-11 | sort -Vu | grep -wF -f - $f > tmp && mv tmp $f; done

for f in *EVE*_heatmap; do awk -F"\t" 'BEGIN{OFS="\t"}; {gsub(/_/,";",$10)}1' $f | sed 's/;/;;/10' | awk -F"\t" 'BEGIN{OFS="\t"}; {gsub(/;;.*/,"",$10)}1' | sed 's/;/_/g' | awk -F"\t" 'BEGIN{OFS="\t"}; {gsub(/_/,";",$8)}1' | sed 's/;/;;/10' | awk -F"\t" 'BEGIN{OFS="\t"}; {gsub(/;;.*/,"",$8)}1' | sed 's/;/_/g' > tmp && mv tmp $f; done

for f in *EVE*_heatmap; do sort -t $'\t' -k 12,12g $f | sed -z 's/^/NAME\tSOURCE\tABBREV\tGROUP\tCONTIG\tSTART\tEND\tVSALLTITLES\tVEVALUE\tSALLTITLES\tEVALUE\tTPM\tCONDITION\tREPEATS\tPFAM\tNUMREADS\n/g' > tmp && mv tmp $f; done
