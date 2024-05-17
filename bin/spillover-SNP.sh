#!/bin/bash

conda activate medaka

#
## map reads vs ref
## already done in spillover-summarise2.sh - bams in:
BAMS="/panfs/jay/groups/27/dcschroe/dmckeow/data/Spillover_FINAL/tmp.mapping.self"

#RUN="01FEB24DM1"
#RUN="15FEB24DM1"
#RUN="18OCTDM23_1"
#RUN="21APR23DM1"
#RUN="26JAN24DM1"
#RUN="29JAN24DM1"
#RUN="30JAN24DM1"
#RUN="060723_mixedlib_PHB_DM"
RUN="091523_BIP23plus"


## model is dna_r10.4.1_e8.2_400bps_hac@v4.2.0
for f in $(ls $BAMS/"$RUN"*.bam); do 
    HDF=$(echo SNP/$(basename $f .bam).hdf)
	FA=$(echo tmp.derep/$(basename $f .bam).contigs.fasta)
	VCF=$(echo SNP/$(basename $f .bam).vcf)
	AVCF=$(echo SNP/$(basename $f .bam).annot.vcf)

	medaka consensus $f $HDF --threads $SLURM_CPUS_PER_TASK --model r1041_e82_400bps_hac_v4.2.0

	medaka variant $FA $HDF $VCF

	medaka tools annotate $VCF $FA $f $AVCF
done

