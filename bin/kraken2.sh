#!/bin/bash -l
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=130GB
#SBATCH --tmp=200GB
#SBATCH --partition agsmall
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

####################### PREREQUISITES #####################################
####################### SCRIPT PURPOSE #####################################

####################### HOW TO RUN SCRIPT #####################################

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################

####################### OUTPUT FILES #####################################

####################### LOAD SOFTWARE #####################################
KRAKEN="/panfs/roc/msisoft/kraken/2.0.7beta/kraken2"
KDB="/panfs/roc/msisoft/kraken/kraken_db"

####################### SET VARIABLES/ARGUMENTS #####################################

################################################################################################
####################### RUNNING THE SCRIPT #####################################################
################################################################################################

#################### STEP 000 #################################
kraken2 --db $KDB Albert_BC03_06_gene_calls.fa \
  --use-mpa-style --threads 8 --use-names \
  --output Albert_BC03_06.results.kraken2 \
  --report Albert_BC03_06.report.kraken2 \
  --unclassified-out Albert_BC03_06.unclassified.kraken2 \
  --classified-out Albert_BC03_06.classified.kraken2 \


##### reformat for anvio import
### gene_callers_id 	t_domain 	t_phylum 	t_class 	t_order 	t_family 	t_genus 	t_species

cut -f 2,3 Albert_BC03_06.results.kraken2 | sed -E 's/(.*)\(taxid ([0-9]+)\)/\1\t\2/g' > tmp1
taxonkit reformat -F -I 3 tmp1 | cut -f 1,4 | sed 's/;/\t/g' | sed -z 's/^/gene_callers_id\tt_domain\tt_phylum\tt_class\tt_order\tt_family\tt_genus\tt_species\n/1' > tmp2
