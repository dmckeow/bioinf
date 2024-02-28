#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 8
#SBATCH --mem 24GB

module load seqkit/0.14.0;
gff="$1"; con="$2"; fas="$3"; ## UBE05_B005_gff, contig for visualisation, genome assembly fasta
abb="/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/NCVOG_prot_name_abbrev";
consize=$(seqkit grep -n -p $con $fas | seqkit stat - | sed -e '1d' -Ee 's/ +/\t/g' -e 's/,//g' | cut -f 8);

cut -f 1,2 $abb | sed -E 's/(ncvog_prot_name=.+)\t/\1_____\t/g' > UBE05_cgview_001;
awk -v con=$con -F"\t|;" 'BEGIN{OFS="\t"}; {if($1 ==con && $3 =="mRNA") print $7,"2",$4,$5,"1","1",$3" "$30,$33; else if ($1 ==con && $3 =="CDS") print $7,"1",$4,$5,"1","0.5",$3" "$30,$33}' $gff | sed -E 's/(ncvog_prot_name=.+)/\1_____/g' | awk 'FNR==NR { array[$1]=$2; next } { for (i in array) gsub(i, array[i]) }1' UBE05_cgview_001 - | sed -e 's/^+\t/forward\t/g' -e 's/^-\t/reverse\t/g' -Ee 's/..._sskingdoms=|ncvog_abbrev=//g' -e 's/\tCDS \t/\tscore_II\t/g' -e 's/\tmRNA Cellular\t/\tpredicted_gene\t/g' -e 's/\tmRNA \t/\tpredicted_gene\t/g' -e 's/\tmRNA Viruses\t/\tpromoter\t/g' -e 's/\tmRNA none\t/\topen_reading_frame\t/g' -e 's/\t$/\t-/g' -e 's/\tncvog_prot_name=NA_____/\t-/g' | sed -z "s/^/\#$con\n%$consize\n!strand\tslot\tstart\tstop\topacity\tthickness\ttype\tlabel\n/g" > UBE05_cgview_$(basename $con).tab;

cut -f 8 UBE05_cgview_$(basename $con).tab | sort -u | sed -E 's/(.+)/ncvog_abbrev=\1_____/g' > UBE05_cgview_$(basename $con .tab)_002;
cut -f 2,3 /projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/NCVOG_prot_name_abbrev | sed -E 's/(\tncvog_abbrev_key=)/_____\1/g' | grep -F -f UBE05_cgview_$(basename $con .tab)_002 - | sort -u | sed -Ee 's/ncvog_abbrev=(.+)_____\tncvog_abbrev_key=(.+)/\1=\2/g' -e 's/_/ /g' -e 's/$/_/g' | tr -d '\n' | sed -e 's/_$/.\n/g' -e 's/_/, /g' > UBE05_cgview_$(basename $con .tab)_003; ## cgview_003 is gene name abbreviation key for figure legend, also details colours to be relabelled
sed -z -i 's/^/### COLOUR KEY SUBSTITUTIONS: ###\nForward gene(red) - REMOVE\nReverse gene(blue) - REMOVE\nscore_II(grey) - exons\npromoter(green) - viral gene\npredicted_gene(yellow) - cellular gene\nopen_reading_frame(pink) - ORFan gene\n\n### GENE ABBREVIATION LEGEND: ###\n/g' UBE05_cgview_$(basename $con .tab)_003; ## add colour details to be changed
###TO VISUALISE: 
### normal visualisation, no labels or legend
#java -jar cgview.jar -i UBE05_cgview_Porterinema-fluviatile_Contig22.tab -o UBE05_cgview_Porterinema-fluviatile_Contig22.svg -f svg -R T -r T
### zoomed visualisation
#java -jar cgview.jar -i UBE05_cgview_Porterinema-fluviatile_Contig22.tab -o UBE05_cgview_Porterinema-fluviatile_Contig22_zoom.svg -f svg -I T -c 1310000 -z 8 -A 10 -D 10 -U 9
