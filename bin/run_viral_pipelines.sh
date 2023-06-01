#!/bin/bash -l

find . -name "barcode*" | awk '/fastq_pass\/barcode[0-9]+$/' > tmp
for f in $(cat tmp); do realpath $f; done > tmp2

cp tmp2 run_viral_pipeline0
sed -E 's/(.+)/sbatch \/home\/dcschroe\/dmckeow\/projects\/DWV\/script\/viral_pipeline1.sh \1 10000/g' tmp2 > run_viral_pipeline1
cp tmp2 run_viral_pipeline2
sed -E 's/(.+)/sbatch \/home\/dcschroe\/dmckeow\/projects\/DWV\/script\/viral_pipeline3.sh \1 10000/g' tmp2 > run_viral_pipeline3

awk '{print "## COPY THE CONTENTS OF THIS FILE TO STEP001 OF viral_pipeline0 ##\nGSIZE=10000\nOUTNAME=$(basename \""$0"\")\n""cd "$0"\n""rm -f tmp_allreads.fastq; cat *.fastq > tmp_allreads.fastq\n""rm -fr DENOVO_\"$OUTNAME\"; mkdir DENOVO_\"$OUTNAME\"\n""cd DENOVO_\"$OUTNAME\"\n""$CANU -p $OUTNAME -d . genomeSize=$GSIZE maxInputCoverage=all corOutCoverage=all corMhapSensitivity=high corMinCoverage=0 redMemory=32 oeaMemory=32 batMemory=64 -nanopore ../tmp_allreads.fastq\nrm -f tmp_checkrun\nuntil [ -s tmp_checkrun ]\ndo\nsqueue -u dmckeow > tmp_checkrun\ngrep -l " canu_" tmp_checkrun | xargs rm -f\nsleep 30\ndone\n"}' run_viral_pipeline0 > tmp && mv tmp run_viral_pipeline0


"OUTNAME=$(basename $0)\ncd "$0"\nrm -f tmp_allreads.fastq; cat *.fastq > tmp_allreads.fastq\nrm -fr DENOVO_"$0"; mkdir DENOVO_"$0"\ncd DENOVO_"$0"\n$CANU -p $(basename $1) -d . genomeSize=$2 maxInputCoverage=all corOutCoverage=all corMhapSensitivity=high corMinCoverage=0 redMemory=32 oeaMemory=32 batMemory=64 -nanopore ../tmp_allreads.fastq




awk '{print $0"\nrm -f tmp_checkrun\nuntil [ -s tmp_checkrun ]\ndo\nsqueue -u dmckeow > tmp_checkrun\ngrep -l " canu_" tmp_checkrun | xargs rm -f\nsleep 30\ndone\n"}' run_viral_pipeline2 > tmp && mv tmp run_viral_pipeline2


#rm -f tmp_checkrun &
#until [ -s tmp_checkrun ]
#do
#  squeue -u dmckeow > tmp_checkrun
#  grep -l " canu_" tmp_checkrun | xargs rm -f
#    sleep 30
#done &
