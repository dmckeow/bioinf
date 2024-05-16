conda activate medaka


## map reads vs ref
## already done in spillover-summarise2.sh - bams in:
BAMS="/panfs/jay/groups/27/dcschroe/dmckeow/data/Spillover_FINAL/tmp.mapping.self"

## model is dna_r10.4.1_e8.2_400bps_hac@v4.2.0
for f in $(ls $BAMS/*.bam); do 
    HDF=$(ls SNP/$(basename $f .bam).hdf)
	FA=$(ls tmp.derep/$(basename $f .bam).contigs.fasta)
	VCF=$(ls SNP/$(basename $f .bam).vcf)

	medaka consensus $f $HDF --threads $SLURM_CPUS_PER_TASK --model r1041_e82_400bps_hac_v4.2.0

	medaka variant $FA $HDF $VCF
done

