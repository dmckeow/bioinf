#!/bin/bash -l
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=600GB
#SBATCH --tmp=300GB
#SBATCH --partition ag2tb
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

####################### PREREQUISITES #####################################

####################### SCRIPT PURPOSE #####################################
## Make database for Kraken2
####################### HOW TO RUN SCRIPT #####################################

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################

####################### OUTPUT FILES #####################################

####################### LOAD SOFTWARE #####################################
#KRAKEN="/panfs/roc/msisoft/kraken/2.0.7beta/kraken2"
#KDB2="/home/dcschroe/dmckeow/kraken2db_nt/"

#SCRATCH="/scratch.global/dcschroe"

## BLAST is INSTALLED LOCALLY AND IN PATH
## KRAKEN2 is INSTALLED LOCALLY AND IN PATH
####################### SET VARIABLES/ARGUMENTS #####################################

################################################################################################
####################### RUNNING THE SCRIPT #####################################################
################################################################################################

#conda activate anvio-7.1
#anvi-script-reformat-fasta Phex_SEP22.fa -o Phex_SEP22.fa2 --seq-type NT


########### KRAKEN2 ####################
########## TO BUILD CUSTOM DATABASE, KRAKEN SCRIPTS MUST BE CHANGED AS FOLLOWS:
#download_genomic_library.sh (line 17)
#FTP_SERVER="ftp://$NCBI_SERVER"
#to
#FTP_SERVER="https://$NCBI_SERVER"

#rsync_from_ncbi.pl (line 46)
#if (! ($full_path =~ s#^ftp://${qm_server}${qm_server_path}/##)) {
#to
#if (! ($full_path =~ s#^https://${qm_server}${qm_server_path}/##)) {

#################### STEP 000 #################################

########## build nt DB
#if [[ "$1" =~ "0" ]] || [[ "$1" =~ "2" ]]; then

### get taxonomy
#kraken2-build --download-taxonomy --db $KDB2
### get specific databases (one line per db)
#kraken2-build --download-library nt --db $KDB2
## make the finalised db
#kraken2-build --build --fast-build --threads 12 --db $KDB2
#fi

########### DRAM ####################
###### setup DRAM database
#conda activate DRAM
#DRAM-setup.py prepare_databases --skip_uniref --output_dir /home/dcschroe/dmckeow/DRAMDB2/
#### test DRAM with its own gene calls
#DRAM.py annotate --threads 12 --min_contig_size 1000 --use_vogdb -i Albert_BC03_06.fa -o DRAM_test1
#DRAM-v.py annotate --threads 12 --min_contig_size 1000 --use_vogdb -i Albert_BC03_06.fa -o DRAM_test1

### test on gene calls done by anvio/prodigal
## issue is that there is no --use_vogdb option, so must be entered as a custom database
#DRAM.py annotate_genes --threads 12 -i Albert_BC03_06_gene_calls.faa -o DRAM_test2 --custom_hmm_name VOGDB --custom_hmm_loc /home/dcschroe/dmckeow/DRAMDB/vog_latest_hmms.txt

###### once DRAM is installed, run the following in terminal:
#### and copy the $DRAMDB
#conda activate DRAM
#DRAMDB=$(DRAM-setup.py print_config | grep "VOGDB db" | sed -e 's/.* \//\//g' -Ee 's/(\/.+)\/.*/\1/g')
#cd $DRAMDB
#rm -f vog.lca.tsv.gz; rm -f vog.lca.tsv;
#wget http://fileshare.csb.univie.ac.at/vog/latest/vog.lca.tsv.gz;
#gunzip vog.lca.tsv.gz;
#cut -f 4 vog.lca.tsv > tmp && mv tmp vog.lca.tsv
#zcat vog_annotations_latest.tsv.gz | paste - vog.lca.tsv | sed -e 's/|/\t\t/g' -e 's/GroupName/VOGDB_id/g' | sed -E -e 's/[^Aa-Zz0-9\t]/_/g' -e 's/_+/_/g' -e 's/^_|_$//g' -e 's/\t_|_\t/\t/g' | sed -e 's/\t\t/__/g' -e 's/\t/;/g' | awk -F ";" '{print "s/___"$1"___/"$0"/g"}' > vog_annotations_latest.lca.reformatted_for_replacement.txt
#conda deactivate


#DRAM.py annotate_genes --threads 12 -i Albert_BC03_06_gene_calls.faa -o DRAM_test2 --custom_hmm_name VOGDB --custom_hmm_loc /home/dcschroe/dmckeow/DRAMDB/vog_latest_hmms.txt

#awk -F "\t" '{print "___"$5"___"}' DRAM_test2/annotations.tsv | sed -e 's/______//g' > tmp_vogid
#sed -i -f $DRAMDB/vog_annotations_latest.lca.reformatted_for_replacement.txt tmp_vogid












DRAMDB="/home/dcschroe/dmckeow/DRAMDB/"

VOGDB="/home/dcschroe/dmckeow/VOGDB"

VOGHMM="$VOGDB/anvioHMMs"; ### the VOG db hmmer profile, downloaded and then additional files created to match the anvio requirements for external hmmer profiles

VOGHMMcutoff="1e-30" #### this is the evalue cut off for VOG HMM scan with anvio


############ SETUP HMMs for VOGDB for Anvi'o

echo "Virus_orthologous_groups" > kind.txt

echo -e "-E $VOGHMMcutoff" > noise_cutoff_terms.txt

echo "https://vogdb.org/" > reference.txt

echo -e "AA:GENE" > target.txt
