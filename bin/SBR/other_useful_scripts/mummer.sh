#!/bin/bash

### script function: DNA whole genome alignment using mummer, specific contig input, self comparison to ID repeat sequences
#### COMMAND LINE inputs required
## $1 = input fasta for alignment
## $2 = reference genome for alignment

####################### INPUTS #####################################

#### software directory
M="/usr/bin/mummer/bin/"
n1=$(basename "$1" | cut -d "." -f 1)
n2=$(basename "$2" | cut -d "." -f 1)

########### NUCMER
#### perform NUCMER alignment
### $5 sets the minimum length of alignments - to low and there will be much noise/complexity
#"$M"nucmer --maxmatch -p "$n1"_"$n2"_nucmer $2 $1

#### for FULL (= all contigs), generate image, png, order contigs, with title, add --color for similarity colour scale
#"$M"mummerplot --postscript -title "$n1"_"$n2" -p "$n1"_"$n2"_nucmer "$n1"_"$n2"_nucmer.delta

########### PROMER
#### perform PROMER alignment
### $5 sets the minimum length of alignments - to low and there will be much noise/complexity
"$M"promer --maxmatch -p "$n1"_"$n2"_promer $2 $1

#### for FULL (= all contigs), generate image, png, order contigs, with title, add --color for similarity colour scale
"$M"mummerplot --postscript -title "$n1"_"$n2" -p "$n1"_"$n2"_promer "$n1"_"$n2"_promer.delta

#### move outputs somewhere
cp *mer.ps /mnt/c/Users/Dean\ Mckeown/Downloads/
