#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 4
#SBATCH --mem 8GB

###### the purpose of this script is to visualise the total percent similarity of a promer alignment between

### command line
## $1 = gff file with NCVOG info
## $2 = dataframe number (record which genome is which df)
## $3 = coords file WITHOUT k option

##### NCVOG defintion files #####
NCVKEY="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/IVEX_NCVOG_006"

##### prep gff to match extracted fasta region (so that loci start from 1 in both)
#s="3498000"; e="3815000"; awk '/^P-fluviatile_contig7\t/' IVEX004_Porterinema-fluviatile.gff | awk -F"\t" -v s=$s -v e=$e 'BEGIN{OFS="\t"};{if($3 ~ "mRNA" && $4 >= s && $5 <= e) print $1,$2,$3,$4-s,$5-s,$6,$7,$8,$9,$10}'

##### prep coords summary of how similar genomes are


awk -F "\t" 'BEGIN{OFS="\t"};{print $1,$2,$12"___"$13}' ../IVEX_WGA/Porterinema-fluviatile_EVE_ALL_EsV-1_FsV-158_Porterinema-fluviatile_EVE_ALL_EsV-1_FsV-158_promer_FULL.coords > tmp1; awk -F "\t" 'BEGIN{OFS="\t"};{if($1 > $2) print $2,$1,$3; else print $0}' tmp1 | sort -t $'\t' -k 3,3V -k 1,2n | uniq > tmp2


##### repeat steps until line count does not change
## remove alignments completely overlapped by another
be=`wc -l tmp2 | sed 's/ .*//g'`
af=0
until [[ "$be" -eq "$af" ]]
  do
      be=`wc -l tmp2 | sed 's/ .*//g'`
      awk -F "\t" 'BEGIN{OFS="\t"};{if(a >= $1 && c == $3) print $1"X",$2,$3; else print $0} {a=$2} {c=$3}' tmp2 | awk -F "\t" 'BEGIN{OFS="\t"};{if(a >= $2 && c == $3) print $1,$2"X",$3; else print $0} {a=$2} {c=$3}' | sed -e '/.*X\t.*X/d' -e 's/X//g' | sort -t $'\t' -k 3,3V -k 1,2n > tmp && mv tmp tmp2
      af=`wc -l tmp2 | sed 's/ .*//g'`
    sleep 1
done
cp tmp2 tmp3

## identify partially overlapping regions and their size (to be subtracted from total)
awk -F "\t" 'BEGIN{OFS="\t"};{if(a >= $1 && c == $3) print $1,$2,$3,$2-$1,a-$1; else print $0,$2-$1,"0"} {a=$2} {c=$3}' tmp3 | sort -t $'\t' -k 3,3V -k 1,2n | awk -F "\t" 'BEGIN{OFS="\t"};{print $1,$2,$4,$5 > $3".tmp4"}'
for f in *.tmp4; do awk '{sum1+=$3; sum2+=$4} END {print sum1-sum2}' $f > tmp && mv tmp $f; done
for f in *.tmp4; do awk -F "\t" '{print FILENAME"\t"$0}' $f; done | sed -E 's/(.+)___(.+)\.tmp4(\t.*)/\1\t\2\3/g' > tmp5
awk -F "\t" 'BEGIN{OFS="\t"};{print $0 > $1".tmp6"}' tmp5
for f in *.tmp6; do max=$(awk -F"\t" '{if($1 ~ $2) print $3}' $f); awk -F "\t"  -v max=$max '{print $0"\t"($3/max)*100}' $f; done | sed -z 's/^/REF\tQUERY\tPROMER_ALIGNMENT_LENGTH_BP\tPERCENT_REF_ALIGNED_TO_QUERY\n/g' > IVEX_genoplotr_EVE_similarity_$(basename $3)

###### R code to make heatmap
## data <- read.delim("IVEX005_EVE_similarity_...", sep = "\t", header = T)

## hm <- ggplot(data, aes(REF, QUERY)) + geom_raster(aes(fill=PERCENT_REF_ALIGNED_TO_QUERY), hjust=0.5, vjust=0.5) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.title = element_text(size=10, hjust = 0.5)) + scale_y_discrete(expand=c(0,0)) + labs(x="reference",y="query") + scale_fill_viridis_c(option = "magma") + geom_text(aes(label = round(PERCENT_REF_ALIGNED_TO_QUERY, 1)), color="grey") + theme(panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.background = element_rect(fill = "transparent"), legend.box.background = element_rect(fill = "transparent"), axis.ticks = element_blank(), plot.title = element_text(color="black", size=16, hjust=0.5, face="bold"))

## ggsave(hm, filename="IVEX005_EVE_similarity_..._heatmap.pdf",height=125,width=200,units="mm",dpi=300,bg="transparent")
