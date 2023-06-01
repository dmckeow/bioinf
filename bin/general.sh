#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=600GB
#SBATCH --partition ag2tb
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

module purge
eval "$(conda shell.bash hook)"

####################### PREREQUISITES #####################################

####################### SCRIPT PURPOSE #####################################
####################### HOW TO RUN SCRIPT #####################################

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################

####################### OUTPUT FILES #####################################

####################### LOAD SOFTWARE #####################################

####################### SET VARIABLES/ARGUMENTS #####################################

################################################################################################


##### gzip parallel stuff
#module load pigz; for f in $(cat $1); do pigz -p $SLURM_CPUS_PER_TASK -d $f ; done

##### split sequences in a multifasta
#conda activate gaas; gaas_fasta_splitter.pl --size_seq 30000000 --overlap 0 --nb_chunks 1 -f $1 -o split_$1


#sed -i -E 's/^(VOG[0-9]+)\t/______\1______\t/g' genes.txt

#zcat genes.hmm.gz | sed -E 's/^(NAME +)(VOG[0-9]+)$/\1______\2______/g' > genes.hmm2


#split -d -l 1000 genes.txt REFORMAT_GENES_

#echo $SLURM_MEM_PER_NODE $SLURM_NTASKS

#for f in REFORMAT_GENES_*; do
 #   srun --ntasks=1 --nodes=1 --cpus-per-task=$SLURM_CPUS_PER_TASK --mem-per-cpu=$(($SLURM_MEM_PER_NODE/$SLURM_NTASKS)) sed -f /panfs/jay/groups/27/dcschroe/dmckeow/VOGDB/vog_annotations.lca.reformatted_for_replacement.txt $f > tmp2_${f} &
#done
#wait

#cat tmp2_REFORMAT_GENES_* > genes.txt3
#rm -fr REFORMAT_GENES_* tmp2_REFORMAT_GENES_*


#sed -f /panfs/jay/groups/27/dcschroe/dmckeow/VOGDB/vog_annotations.lca.reformatted_for_replacement.txt genes.txt > genes.txt2

#sed -f /panfs/jay/groups/27/dcschroe/dmckeow/VOGDB/vog_annotations.lca.reformatted_for_replacement.txt genes.hmm2 > genes.hmm3
#gzip genes.hmm3


#sed -f replace_ncvog genes.hmm2 > genes.hmm3



######### KAIJU DB
##conda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/kaiju
##kaiju-makedb -s nr_euk -t $SLURM_CPUS_PER_TASK
##kaiju-makedb -s rvdb -t $SLURM_CPUS_PER_TASK


cp /scratch.global/dcschroe/kaijudb_nr/nr_euk/kaiju_db_nr_euk.fmi /home/dcschroe/shared/tools/kaijudb/nr_euk

