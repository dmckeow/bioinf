#!/bin/bash

### script function: DNA whole genome alignment using mummer, specific contig input, self comparison to ID repeat sequences
#### COMMAND LINE inputs required
## $1 = input genomes assembly
## $2 = contig name to align with itself
## $3 = contig prefix to remove (leaving only identifiable number on y axis)
## $4 = title (at top of figure)
## $5 = minimum nucmer alignment length (1500 good start)

####################### INPUTS #####################################

#### software directory
M="/usr/bin/mummer/bin/"
n1=$(basename "$1" | cut -d "." -f 1 | sed 's/^/tmp_/g')

### RUN in a subdirectory in home directory
seqkit grep -n -p "$2" "$1" > "$n1"__$(basename "$2").fa

#### reduce contig names to identifable numbers
echo "s/$3//g" > tmp1
sed -i -f tmp1 "$n1"__$(basename "$2").fa

########### NUCMER
#### perform NUCMER alignment
### $5 sets the minimum length of alignments - to low and there will be much noise/complexity
"$M"nucmer --maxmatch --nosimplify -L $5 -p "$n1"__$(basename "$2" .fa)_nucmer "$n1"__$(basename "$2").fa "$n1"__$(basename "$2").fa

#### for FULL (= all contigs), generate image, png, order contigs, with title, similarity colour scale
"$M"mummerplot --postscript -title $4 --color -p "$n1"__$(basename "$2" .fa)_nucmer "$n1"__$(basename "$2" .fa)_nucmer.delta

##### get repeat coords and ref
### get coords and remove all in-place self-hits
"$M"show-coords -HTd "$n1"__$(basename "$2" .fa)_nucmer.delta | awk '{if($1 == $3) gsub(/.+/,"")}1' | sed '/^$/d' > tmp2
### remove duplicated entries forward/reverse
awk -F "\t" 'BEGIN{OFS="\t"};{if($3 > $4) print $1,$2,$4,$3,$5,$6,$7,$8,$9,$10,$11; else print $0}' tmp2 | awk -F "\t" 'BEGIN{OFS="\t"};{if($1 > $3) print $3,$4,$1,$2,$5,$6,$7,$8,$9,$10,$11; else print $0}' | sort -u -t $'\t' -k 1,4n > tmp3

### assign broad categories to repeats (direct/inverted, interspersed/tandem)
awk -F "\t" 'BEGIN{OFS="\t"};{if($8 =="1" && $9 =="-1") print $0,"INV_"; else print $0,"DIR_"}' tmp3 | awk -F "\t" 'BEGIN{OFS="\t"};{if($2 > $3) print $0"TAN"; else print $0"INT"}' > tmp4

### rearrange and reformat into gff style
awk -F "\t" -v a=$2 'BEGIN{OFS="\t"};{print a,"nucmer",$12,$1,$2,".",$8,".","ID="a"_repeat"NR";Name="NR";length="$5";pid="$7}' tmp4 | sed -e 's/\.\t1\t\./\.\t+\t\./g' -e 's/\.\t-1\t\./\.\t-\t\./g' > tmp5
awk -F "\t" -v a=$2 'BEGIN{OFS="\t"};{print a,"nucmer",$12,$3,$4,".",$9,".","ID="a"_repeat"NR";Name="NR";length="$6";pid="$7}' tmp4 | sed -e 's/\.\t1\t\./\.\t+\t\./g' -e 's/\.\t-1\t\./\.\t-\t\./g' >> tmp5

mv tmp5 IVEX_nucmerRep_"$n1"__$(basename "$2" .fa).gff

### get fasta of repeats
awk -F"\t|;" 'BEGIN{OFS="\t"};{gsub(/$/,$10,$3)}1' IVEX_nucmerRep_"$n1"__$(basename "$2" .fa).gff | sed -e 's/\t/;/9g' -e 's/Name=/_/1' | sed -f tmp1 - | bedtools getfasta -fi "$n1"__$(basename "$2").fa -bed - -fo IVEX_nucmerRep_"$n1"__$(basename "$2").fa -name

############ CGview prep
##### reformat feature labels for CGview
sed -e 's/Name=\|ID=//g' -e 's/;/\t/g' -e 's/\tDIR_TAN\t/\trRNA\t/g' -e 's/\tINV_TAN\t/\t\t/g' -e 's/\tDIR_INT\t/\tCDS\t/g' -e 's/\tINV_INT\t/\ttRNA\t/g' IVEX_nucmerRep_"$n1"__$(basename "$2" .fa).gff | sort -t $'\t' -k 4,5n > tmp1

##### prepare gff-style input for CGview -genes option
awk -F "\t" 'BEGIN{OFS="\t"};{print $9,".",$3,$4,$5,".",$7,"."}' tmp1 | sed -z 's/^/seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\n/g' > tmp_genes

##### prepare list input for CGview -labels_to_show option (to label NCVOG abbreviations) if($22 !="")
awk -F "\t" 'BEGIN{OFS="\t"};{print $9,$10}' tmp1 > tmp_labels_to_show

##### cgview_xml_builder to create XML input
perl /home/dmckeown/cgview/scripts/cgview_xml_builder/cgview_xml_builder.pl -sequence "$n1"__$(basename "$2").fa -genes tmp_genes -labels_to_show tmp_labels_to_show -output tmp.xml -linear T -draw_divider_rings T -gene_labels T -custom rRNAColor="rgb(128,0,128)" tRNAColor="rgb(65,105,225)" otherColor="rgb(33,144,141)" proteinColor="rgb(90,200,101)" gcSkewColorNeg="rgb(255,0,0)" gcSkewColorPos="rgb(0,100,0)"

##### XML edits
### remove shading on gene arrows and linear topology text
sed -i 's/showShading=\"true\"/showShading=\"false\"/g' tmp.xml
sed -i '/legendItem text=\"Topology: linear\"/d' tmp.xml
### change legend text, positions
sed -i 's/legendItem text=\"CDS\"/legendItem text=\"direct interspersed\"/g' tmp.xml
sed -i 's/legendItem text=\"Length: /legendItem text=\" /g' tmp.xml
sed -i 's/legendItem text=\"tRNA\"/legendItem text=\"inverted interspersed\"/g' tmp.xml
sed -i 's/legendItem text=\"rRNA\"/legendItem text=\"direct tandem\"/g' tmp.xml
sed -i 's/legendItem text=\"\"/legendItem text=\"inverted tandem\"/g' tmp.xml
sed -i 's/legend position=\"upper-left\"/legend position=\"upper-center\"/g' tmp.xml
sed -i 's/legend position=\"upper-right\"/legend position=\"lower-right\"/g' tmp.xml

########## change colours

########### add new colors to legend and fix opacity
### ORFan colour legend
sed -i -E 's/(^<legend position="lower-right".+)/\1\n<legendItem text=\"inverted tandem\" drawSwatch=\"true\" swatchOpacity="0.5" swatchColor=\"rgb\(33,144,141\)\" \/>/g' tmp.xml

### opacity of colours arrows and legend
sed -i -e 's/opacity=\"0\.5\"/opacity="0.9"/g' -e 's/swatchOpacity=\"0\.5\"/swatchOpacity="0.9"/g' tmp.xml

##### run CGview to generate image (svg)
java -jar -Xmx1500m /home/dmckeown/cgview/bin/cgview.jar -i tmp.xml -o IVEX_nucmerRep_"$n1"__$(basename "$2" .fa).svg -f svg -A 32 -D 32 -I F

#### final rename, move to results, and delete tmp file
for f in "$n1"__$(basename "$2" .fa)_nucmer.ps; do mv "$f" "${f//tmp_/IVEX_nucmerRep_}"; done
for f in IVEX_nucmerRep_tmp*; do mv "$f" "${f//IVEX_nucmerRep_tmp_/IVEX_nucmerRep_}"; done
mv IVEX_nucmerRep_* /mnt/c/Users/Dean\ Mckeown/Documents/PC_iSSD/SBR_backup/finalresult/IVEX_WGA/
#rm -f tmp*
