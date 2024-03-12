#!/bin/bash

conda activate bioinftools

step="a4"

##########################################
if [[ "${step}" =~ "A" ]] || [[ "$step" =~ "a1" ]]; then
###########################################

rm -fr tmp.CATBAT; mkdir tmp.CATBAT

for f in $(cat bin_list | sed '1,1d' | cut -f 1  | sort -Vu); do
	mkdir tmp.CATBAT/$f
done


#### prepare contigs fasta in bin folders for CATBAT
while IFS=$'\t' read -r AnvioBin ContigNameBin ContigPath ReadPath SampleNameBin ContigNameCanu ContigNameDerep; do
	echo "touch tmp.CATBAT/$AnvioBin/$(basename $ContigPath .contigs.fasta).fna"
done < bin_list | sed '1,1d' | sort -Vu > tmp.CATBAT/RunFolderPrep
bash tmp.CATBAT/RunFolderPrep

#### prepare contigs fasta in bin folders for CATBAT
while IFS=$'\t' read -r AnvioBin ContigNameBin ContigPath ReadPath SampleNameBin ContigNameCanu ContigNameDerep; do
	echo "seqkit grep -p $ContigNameCanu $ContigPath | sed 's/>/>$(basename $ContigPath)______/g' >> tmp.CATBAT/$AnvioBin/$(basename $ContigPath .contigs.fasta).fna"
done < bin_list | sed '1,1d' | sort -Vu > tmp.CATBAT/RunFastaCopy
bash tmp.CATBAT/RunFastaCopy

find tmp.CATBAT -type f -empty -print -delete

ls -d tmp.CATBAT/BINNING_* > tmp.CATBAT/CATBATBinList

while read -r BinFolder; do
	rm -f $BinFolder/$(basename $BinFolder).concatenated.fasta; touch $BinFolder/$(basename $BinFolder).concatenated.fasta
	cat $BinFolder/*.fna >> $BinFolder/$(basename $BinFolder).concatenated.fasta
done < tmp.CATBAT/CATBATBinList

#####################################
fi
#####################################


##########################################
if [[ "${step}" =~ "A" ]] || [[ "$step" =~ "a2" ]]; then
###########################################


ls -d tmp.CATBAT/BINNING_*/*.concatenated.fasta > tmp.CATBAT/CATBATConcatList
while read -r Concat; do
	prodigal -i $Concat -a $Concat.predicted_proteins.faa -o $Concat.predicted_proteins.gff -p meta -g 11 -q -f gff
done < tmp.CATBAT/CATBATConcatList


#####################################
fi
#####################################


##########################################
if [[ "${step}" =~ "A" ]] || [[ "$step" =~ "a3" ]]; then
###########################################

rm -f tmp.CATBAT/AllBins.predicted_proteins.faa; touch tmp.CATBAT/AllBins.predicted_proteins.faa

while read -r BinFolder; do
	sed "s/>/>$(basename $BinFolder)_________/g" $BinFolder/*.predicted_proteins.faa >> tmp.CATBAT/AllBins.predicted_proteins.faa
done < tmp.CATBAT/CATBATBinList

diamond blastp --fast -p 24 -d /panfs/jay/groups/27/dcschroe/shared/20231120_CAT_nr/db/2023-11-21_CAT.dmnd -q tmp.CATBAT/AllBins.predicted_proteins.faa -o tmp.CATBAT/AllBins.predicted_proteins.diamond

#####################################
fi
#####################################

##########################################
if [[ "${step}" =~ "A" ]] || [[ "$step" =~ "a4" ]]; then
###########################################

#while read -r BinFolder; do
#	rm -f $BinFolder/$(basename $BinFolder).predicted_proteins.diamond; touch $BinFolder/$(basename $BinFolder).predicted_proteins.diamond
#done < tmp.CATBAT/CATBATBinList

#while read -r bin qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore; do
#	echo -e "$qseqid\t$sseqid\t$pident\t$length\t$mismatch\t$gapopen\t$qstart\t$qend\t$sstart\t$send\t$evalue\t$bitscore" >> tmp.CATBAT/$bin/$bin.predicted_proteins.diamond
#done < <(awk -F "\t" 'BEGIN{OFS=FS} { gsub(/_________/,"\t", $1) }1' tmp.CATBAT/AllBins.predicted_proteins.diamond)

while read -r BinFolder; do
	/panfs/jay/groups/27/dcschroe/dmckeow/CAT_pack/CAT_pack/CAT_pack bins -b $BinFolder \
	-d /panfs/jay/groups/27/dcschroe/shared/20231120_CAT_nr/db \
	-t /panfs/jay/groups/27/dcschroe/shared/20231120_CAT_nr/tax \
	-o CATBAT_$(basename $BinFolder) --no_log --no_stars --force \
	-p $BinFolder/$(basename $BinFolder).concatenated.fasta.predicted_proteins.faa \
	-a $BinFolder/$(basename $BinFolder).predicted_proteins.diamond
done < tmp.CATBAT/CATBATBinList

#####################################
fi
#####################################
