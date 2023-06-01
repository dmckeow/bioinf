#!/bin/bash

### mummer is not available on the cluster, so this script must be edited to your specific PC directory setup / installation
### running this script is recommended, as it can identify good viral regions that may not stand out amongst the blast result summary
### the coords output of this script are also requried for the IVEX_genoplotr scripts

#### COMMAND LINE inputs required
## $1 = query genome assembly (ensure fasta deflines are what you want on figure - no spaces, not gzipped)
## $2 = reference genome assembly (ensure fasta deflines are what you want on figure - no spaces, not gzipped)
## $3 = genome name for title (at top of figure)


#### software directory
M="/usr/bin/mummer/bin/"
n1=$(basename "$1" | cut -d "." -f 1)
n2=$(basename "$2" | cut -d "." -f 1)

############ STEP_001 #####################
########### PROMER
#### for TMP perform PROMER alignment try maxmatch or --mum
"$M"promer --maxmatch -c 100 -p "$n1"_"$n2"_promer_TMP $2 $1

#### for TMP (= all contigs), generate image, png, order contigs, with title, similarity colour scale
"$M"mummerplot --color --postscript -title $3 -p "$n1"_"$n2"_promer_TMP "$n1"_"$n2"_promer_TMP.delta

### reduce to only aligned contigs
grep ">" "$n1"_"$n2"_promer_TMP.delta | sed 's/ /\t/g' | cut -f 2 | sort -Vu | seqkit grep -n -f - $1 > FINAL_"$n1".fa

#### for FINAL perform PROMER alignment
"$M"promer -p "$n1"_"$n2"_promer_FINAL $(basename "$2") FINAL_"$n1".fa

#### for FINAL (= well-aligned contigs), generate image, png, order contigs, with title, similarity colour scale
"$M"mummerplot --layout --color --postscript -title $3 -p "$n1"_"$n2"_promer_FINAL "$n1"_"$n2"_promer_FINAL.delta

"$M"show-coords -HTd "$n1"_"$n2"_promer_FINAL.delta > "$n1"_"$n2"_promer_FINAL.coords

############ STEP_002 #####################

######## SUMMARISE MUMMER PROMER RESULTS FOR SUMMARY
#### get contigs as listed in ../IVEX002/IVEX002_final_005_contig_summary and flag with y if aligned to reference (repeat for whatever reference aligned in promer)
#### run these steps interactively once STEP_001 ran on all genomes
#### modify these steps to match your reference genomes and repeat for each reference genome (ensure FINAL output file name is changed to specify reference genome)
##for f in *_Ectocarpus_siliculosus_virus_1_promer_FINAL.coords; do awk -F"\t" -v f=$f '{print "s/"f";"$13"___/Y/g"}' $f | sort -Vu; done | sed 's/_Ectocarpus_siliculosus_virus_1_promer_FINAL.coords//g' > IVEX_promer_TMP_fix1
##awk -F ";" '{print $1";"$2"___"}' ../IVEX002/IVEX002_final_005_contig_summary | sed '1,1d' > IVEX_promer_TMP_fix2
##sed -f IVEX_promer_TMP_fix1 IVEX_promer_TMP_fix2 | sed 's/.*___$/n/g' > IVEX_promer_TMP_fix3;
##sed 's/___$//g' IVEX_promer_TMP_fix2 | paste -d ";" - IVEX_promer_TMP_fix3 > IVEX_promer_FINAL_esv1

#### when done, as in you have generated all summary results for genome, move to results, e.g.
mv *promer_FINAL* /mnt/c/Users/Dean\ Mckeown/Documents/PC_iSSD/SBR_backup/finalresult/IVEX_WGA/
rm -f *_promer_TMP*
