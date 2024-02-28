#!/bin/bash

### mummer is not available on the cluster, so this script must be edited to your specific PC directory setup / installation

#### COMMAND LINE inputs required
## $1 = query genome assembly (can be multifasta, ensure fasta deflines are what you want on figure - no spaces, not gzipped)
## $2 = reference genome assemlby (ensure fasta deflines are what you want on figure - no spaces, not gzipped)
## $3 = genome name for title (at top of figure)


#### software directory
M="/usr/bin/mummer/bin/"
n1=$(basename "$1" | cut -d "." -f 1)
n2=$(basename "$2" | cut -d "." -f 1)

########### NUCMER
#### for FULL perform NUCMER alignment
"$M"nucmer --maxmatch --nosimplify -p "$n1"_"$n2"_nucmer_FULL $2 $1

#### for FULL (= all contigs), generate image, png, order contigs, with title, similarity colour scale
"$M"mummerplot --color --postscript -title $3 -p "$n1"_"$n2"_nucmer_FULL "$n1"_"$n2"_nucmer_FULL.delta

### reduce to only well-aligned contigs
#grep ">" "$n1"_"$n2"_nucmer_FULL.delta | sed 's/ /\t/g' | cut -f 2 | sort -Vu | seqkit grep -n -f - $1 > PART_"$n1".fa

#### for PART perform NUCMER alignment
#"$M"nucmer -p "$n1"_"$n2"_nucmer_PART $(basename "$2") PART_"$n1".fa

#### for PART (= well-aligned contigs), generate image, png, order contigs, with title, similarity colour scale
#"$M"mummerplot --postscript -title $3 -p "$n1"_"$n2"_nucmer_PART "$n1"_"$n2"_nucmer_FULL.delta

#### move to results
mv *nucmer*.ps /mnt/c/Users/Dean\ Mckeown/Documents/PC_iSSD/SBR_backup/finalresult/IVEX_WGA/
mv *nucmer*.png /mnt/c/Users/Dean\ Mckeown/Documents/PC_iSSD/SBR_backup/finalresult/IVEX_WGA/
rm -f *nucmer*
