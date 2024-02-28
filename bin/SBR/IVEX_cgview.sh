#!/bin/bash

### the purpose of this script is to generate circular genome/contig visualisations using cgview
### CGview is not available on the cluster, so this script must be edited to your specific PC directory setup / CGview installation
### this script does not require editing - the output can be changed using the command line arguments below
### if you wish to edit the colours and appearance, it is possible by editing the script, but this will likely require studying the CGview tutorials

##### COMMAND LINE INPUTS REQUIRED
## $1 = contig name to visualise (must match contig name in gff and fasta)
## $2 = genome annotation .gff.gz (IVEX002 or IVEX004)
## $3 = genome assembly .fa.gz
## $4 = zoom multiplier (use 1 is unsure/not needed) OR start loci of region of interest
## $5 = centre loci for zoom multiplier (use 1 if not needed) OR end loci of region of interest
## $6 = legend off T or F
## $7 = gene labels off T or F
## $8 = show GC content and skew T or F
## $9 = zoom style A or B - A is zoomed on area of circular contig, B extracts region of interest and circularises it alone
## e.g. /mnt/c/Users/Dean\ Mckeown/Documents/PC_iSSD/SBR_backup/script/IVEX_cgview.sh P-fluviatile_contig22 /mnt/c/Users/Dean\ Mckeown/Documents/PC_iSSD/SBR_backup/finalresult/IVEX002/Porterinema-fluviatile.gff.gz /mnt/c/Users/Dean\ Mckeown/Documents/PC_iSSD/SBR_backup/finalresult/IVEX002/Porterinema-fluviatile.fa.gz 1200000 1460000 F F F B

## name variable (incase of different file extensions)
n2=$(basename "$2" | cut -d "." -f 1)
n3=$(basename "$3" | cut -d "." -f 1)

##### reduce gff to contig of interest, reformat feature labels for CGview
if grep -q ".gz" "$2"; then
  zcat "$2" | grep -wF "$1" - | awk -F "\t" 'BEGIN{OFS="\t"};{if($3 =="CDS" || $3 =="mRNA_cellular" || $3 =="mRNA_cellular_or_viral" || $3 =="mRNA_viral" || $3 =="mRNA_ORFan") print $0}' | sed -e 's/;/\t/g' -e 's/\tCDS\t/\tother\t/g' -e 's/\tmRNA_cellular_or_viral\t/\trRNA\t/g' -e 's/\tmRNA_cellular\t/\ttRNA\t/g' -e 's/\tmRNA_viral\t/\tCDS\t/g' -e 's/\tmRNA_ORFan\t/\t\t/g' | awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/.+=/,$3"_",$9)}1' | awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/.+=|.+=-/,"",$22)}1' | sort -t $'\t' -k 4,5n > tmp1
fi

if ! grep -q ".gz" "$2"; then
  grep -wF "$1" "$2" | awk -F "\t" 'BEGIN{OFS="\t"};{if($3 =="CDS" || $3 =="mRNA_cellular" || $3 =="mRNA_cellular_or_viral" || $3 =="mRNA_viral" || $3 =="mRNA_ORFan") print $0}' | sed -e 's/;/\t/g' -e 's/\tCDS\t/\tother\t/g' -e 's/\tmRNA_cellular_or_viral\t/\trRNA\t/g' -e 's/\tmRNA_cellular\t/\ttRNA\t/g' -e 's/\tmRNA_viral\t/\tCDS\t/g' -e 's/\tmRNA_ORFan\t/\t\t/g' | awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/.+=/,$3"_",$9)}1' | awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/.+=|.+=-/,"",$22)}1' | sort -t $'\t' -k 4,5n > tmp1
fi

##### prepare gff-style input for CGview -genes option
awk -F "\t" 'BEGIN{OFS="\t"};{print $9,".",$3,$4,$5,".",$7,"."}' tmp1 > tmp2
### list lines for other(true CDS) for multiexonic genes only
awk '/^other_/' tmp2 | cut -f 1 | sort -V | uniq -d | grep -wF -f - tmp2 > tmp3
### remove all monoexonic other(true CDS) lines and add back multiexonic ones
sed '/^other_/d' tmp2 | cat - tmp3 | sort -t $'\t' -k 4,5n | sed -z 's/^/seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\n/g' > tmp_IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9"_genes
### reduce all CDS sizes down to x bp (so they appear as bands)
awk -F"\t" 'BEGIN{OFS="\t"; OFMT="%.0f"};{if($3 =="other") print $1,$2,$3,$4+($5-$4)/2,$4+200+($5-$4)/2,$6,$7,$8; else print $0}' tmp_IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9"_genes > tmp && mv tmp tmp_IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9"_genes

##### make reduced fastas and prepare list input for CGview -labels_to_show option (to label NCVOG abbreviations) if($22 !="")
### option A
if [[ "$9" =~ "A" ]]; then
  seqkit grep -n -p "$1" "$3" > tmp_"$n2"_"$1".fa
  awk -F "\t" 'BEGIN{OFS="\t"};{if($22 ~ /.+/ && $22 !="u") print $9,$22}' tmp1 > tmp_IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9"_labels_to_show
fi

### option B
if [[ "$9" =~ "B" ]]; then
  seqkit grep -n -p "$1" "$3" | seqkit subseq -r $4:$5 > tmp_"$n2"_"$1".fa
  awk -F "\t" -v a="$5" 'BEGIN{OFS="\t"};{if($4 <= a && $5 <= a) print $0}' tmp_IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9"_genes | awk -F "\t" -v a="$4" 'BEGIN{OFS="\t"};{print $1,$2,$3,$4-a,$5-a,$6,$7,$8}' | awk -F "\t" 'BEGIN{OFS="\t"};{if($4 >= 1 && $5 >= 1) print $0}' | sed -z 's/^/seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\n/g' > tmp && mv tmp tmp_IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9"_genes
  cut -f 1 tmp_IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9"_genes | grep -wF -f - tmp1 | awk -F "\t" 'BEGIN{OFS="\t"};{if($22 ~ /.+/ && $22 !="u") print $9,$22}' > tmp_IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9"_labels_to_show
fi

### define colours
ora="255,140,0" ## orange
mgr="102,102,102" ## medium grey
bla="0,0,0" ## black
ver="128,0,0" ## vermilion
red="255,0,0" ## red
oli="0,100,0" ## olive green

lgr="204,204,204" ## light grey
mau="128,0,128" ## mauve
blu="65,105,225" ## blue
tur="33,144,141" ## turquoise
gre="90,200,101" ## green
yel="249,231,33" ## yellow

##### cgview_xml_builder to create XML input
perl /home/dmckeown/cgview/scripts/cgview_xml_builder/cgview_xml_builder.pl -sequence tmp_"$n2"_"$1".fa -genes tmp_IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9"_genes -labels_to_show tmp_IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9"_labels_to_show -output IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml -linear T -draw_divider_rings T -gene_labels T -gc_content $8 -gc_skew $8 -log IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".log -custom rRNAColor="rgb($ora)" tRNAColor="rgb($mgr)" otherColor="rgb($bla)" proteinColor="rgb($ver)" gcSkewColorNeg="rgb($red)" gcSkewColorPos="rgb($oli)"

##### XML edits
### remove shading on gene arrows and linear topology text
sed -i 's/showShading=\"true\"/showShading=\"false\"/g' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml
sed -i '/legendItem text=\"Topology: linear\"/d' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml
### change legend text, positions
sed -i 's/legendItem text=\"CDS\"/legendItem text=\"viral\"/g' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml
sed -i 's/legendItem text=\"Length: /legendItem text=\" /g' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml
sed -i 's/legendItem text=\"rRNA\"/legendItem text=\"viral or cellular\"/g' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml
sed -i 's/legendItem text=\"tRNA\"/legendItem text=\"cellular\"/g' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml
sed -i 's/legendItem text=\"Other\"/legendItem text=\"ORFan\"/g' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml
sed -i 's/legend position=\"upper-left\"/legend position=\"upper-center\"/g' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml
sed -i 's/legend position=\"upper-right\"/legend position=\"lower-right\"/g' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml

########## change colours
### replace ORFan colour so that it is not assigned other colour
awk -v co="$lgr" '{if($0 ~ /^<feature color="rgb\(0,0,0\)"/ && $0 ~ /label="_/) gsub(/feature color="rgb\([0-9]+,[0-9]+,[0-9]+\)"/,"feature color=\"rgb("co")\"")}1' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml > tmp && mv tmp IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml

### group x NCVOGs to viridis colour scheme
awk -v co="$mau" '{if($0 ~ /^<feature color=/ && $0 ~ /\(1\)/) gsub(/feature color="rgb\([0-9]+,[0-9]+,[0-9]+\)"/,"feature color=\"rgb("co")\"") gsub(/\(1\)/,"")}1' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml > tmp && mv tmp IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml

awk -v co="$blu" '{if($0 ~ /^<feature color=/ && $0 ~ /\(2\)/) gsub(/feature color="rgb\([0-9]+,[0-9]+,[0-9]+\)"/,"feature color=\"rgb("co")\"") gsub(/\(2\)/,"")}1' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml > tmp && mv tmp IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml

awk -v co="$tur" '{if($0 ~ /^<feature color=/ && $0 ~ /\(3\)/) gsub(/feature color="rgb\([0-9]+,[0-9]+,[0-9]+\)"/,"feature color=\"rgb("co")\"") gsub(/\(3\)/,"")}1' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml > tmp && mv tmp IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml

awk -v co="$gre" '{if($0 ~ /^<feature color=/ && $0 ~ /\(4\)/) gsub(/feature color="rgb\([0-9]+,[0-9]+,[0-9]+\)"/,"feature color=\"rgb("co")\"") gsub(/\(4\)/,"")}1' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml > tmp && mv tmp IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml

awk -v co="$yel" '{if($0 ~ /^<feature color=/ && $0 ~ /\(5\)/) gsub(/feature color="rgb\([0-9]+,[0-9]+,[0-9]+\)"/,"feature color=\"rgb("co")\"") gsub(/\(5\)/,"")}1' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml > tmp && mv tmp IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml

########### add new colors to legend and fix opacity
### ORFan colour legend
sed -i -E 's/(^<legendItem text="ORFan".+swatchColor=)"rgb\([0-9]+,[0-9]+,[0-9]+\)" \/>/\1"rgb\('"$lgr"'\)" \/>/g' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml

### NCVOG colour legend
sed -i -E 's/(^<legend position="lower-right".+)/\1\n<legendItem text=\"NCVOG group 1\" drawSwatch=\"true\" swatchOpacity="0.5" swatchColor=\"rgb\('"$mau"'\)\" \/>/g' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml
sed -i -E 's/(^<legend position="lower-right".+)/\1\n<legendItem text=\"NCVOG group 2\" drawSwatch=\"true\" swatchOpacity="0.5" swatchColor=\"rgb\('"$blu"'\)\" \/>/g' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml
sed -i -E 's/(^<legend position="lower-right".+)/\1\n<legendItem text=\"NCVOG group 3\" drawSwatch=\"true\" swatchOpacity="0.5" swatchColor=\"rgb\('"$tur"'\)\" \/>/g' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml
sed -i -E 's/(^<legend position="lower-right".+)/\1\n<legendItem text=\"NCVOG group 4\" drawSwatch=\"true\" swatchOpacity="0.5" swatchColor=\"rgb\('"$gre"'\)\" \/>/g' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml
sed -i -E 's/(^<legend position="lower-right".+)/\1\n<legendItem text=\"NCVOG Phaeovirus\" drawSwatch=\"true\" swatchOpacity="0.5" swatchColor=\"rgb\('"$yel"'\)\" \/>/g' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml

### opacity of colours arrows and legend
sed -i -e 's/opacity=\"0\.5\"/opacity="0.9"/g' -e 's/swatchOpacity=\"0\.5\"/swatchOpacity="0.9"/g' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml

##### run CGview to generate image (svg)
### option A To make circular map of a whole contig OR a zoomed in region of contig
if [[ "$9" =~ "A" ]]; then
  java -jar -Xmx1500m /home/dmckeown/cgview/bin/cgview.jar -i IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml -o IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".svg -f svg -U 32 -A 48 -D 48 -I T -z $4 -c $5 -r $6 -R $7
fi

### option B to make a circular map of extracted contig region
if [[ "$9" =~ "B" ]]; then
  java -jar -Xmx1500m /home/dmckeown/cgview/bin/cgview.jar -i IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".xml -o IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".svg -f svg -U 32 -A 48 -D 48 -I T -z 1 -c 1 -r $6 -R $7
fi

##### prepare gene abbreviation key for legend
cut -f 1 tmp_IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9"_labels_to_show | grep -wF -f - tmp1 | awk -F "\t" 'BEGIN{OFS="\t"};{if($22 ~ /.+/ && $22 !="u") print $22,$23}' | sort -Vu | sed -e 's/ncvog_legend=//g' -e '/^\t/d' -e '/^$/d' -e 's/\t/: /g' -e 's/_/ /g' -e 's/$/, /g' | tr -d '\n' | sed -e 's/, $/\n/g' > IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9"_legend.txt

#### extract infor on GC values (if plotted)
if [[ "$8" =~ "T" ]]; then
  awk '/GC skew | GC content/' IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9".log | awk '/value is|window/' | sed -Ee 's/^Plotting.+ ([0-9]+).+ ([0-9]+)/GC values plotted with \1 bp window and \2 step/1' -e 's/ value is /;/g' -e '5,5d' -e 's/ GC content| GC skew//g' -e 's/The maximum;/max /g' -e 's/The minimum;/min /g' -e 's/The average;/average /g' -e 's/\.$/, /g' | tr -d '\n' | sed -e 's/, $/.\n/g' -e 's/max/GC content: max/1' -e 's/, max/; GC skew: max/1' | cat IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9"_legend.txt - > tmp && mv tmp IVEX_cgview_"$n2"_"$1"_"$4"_"$5"_"$9"_legend.txt
fi

mv IVEX_cgview_* /mnt/c/Users/Dean\ Mckeown/Documents/PC_iSSD/SBR_backup/finalresult/IVEX_cgview/
rm -f tmp*
