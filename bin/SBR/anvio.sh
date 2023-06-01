#!/bin/bash -l
#SBATCH --cpus-per-task=24
#SBATCH --mem=128GB
#SBATCH --partition long
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

####################### PREREQUISITES #####################################
## nanopore seq reads data

####################### SCRIPT PURPOSE #####################################

####################### HOW TO RUN SCRIPT #####################################

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################

## if you intend to merge all samples together for analysis and visualisation in avi'o, then you must use either a co-assembly (multilpe sequencing runs and one contig fasta assembly) OR a concatenated fasta of independently assembled contigs

####################### OUTPUT FILES #####################################

####################### LOAD SOFTWARE #####################################
### LOAD software available via shared environment on server:
#module load conda
module purge
module load minimap2/2.17
module load samtools
module load seqkit

### set paths for DATABASES that must be installed locally:
#KDB1="/scratch.global/dcschroe/kraken2db_ABVFPP"; ### a DB comprised of reference genomes for Archaea, Bacteria, Viruses, Fungi, Plants, and Protists: see Kraken2 database bulding for instructions. NOTE - this will be very large (~> 1 TB), so do it in scratch
#KDB2="/scratch.global/dcschroe/kraken2db_nt"; ### a DB comprised of NCBI's nt: see Kraken2 database bulding for instructions. NOTE - this will be very large (~> 1 TB), so do it in scratch
#VOGHMM="/home/dcschroe/dmckeow/VOGDB_HMMs_anvio"; ### the VOG db hmmer profile, downloaded and then additional files created to match the anvio requirements for external hmmer profiles
#NCLDVHMM="/home/dcschroe/dmckeow/HMMs_NCLDVs/HMM_NCLDV_149_core_genes" ### a hmmer profile made by the anvio team for NCLDVs

######### SOFTWARE/OTHER SCRIPTS THAT NEED SETUP BEFORE RUNNING:

#VS2_ANVIO="/home/dcschroe/dmckeow/VirSorter2_to_Anvio/virsorter_to_anvio.py" ### a custom script made by the VS2 developers to reformat VS2 results for import into anvio
##### Seqkit:
### must be installed locally, e.g.:
#SEQKIT="/home/dcschroe/dmckeow/seqkit"

##### KRAKEN2:
## KRAKEN2 must be installed locally, and be in your $PATH

##### ANVIO 7.1 or later:
## ANVIO is accessed throughout script using conda activate anvio-7.1 - it must be installed locally via conda

##### DRAM:
## DRAM is accessed throughout script using conda activate DRAM - it must be installed locally via conda
### before runing script, and after DRAM installation, the following needs to be run once (it adds more meaningful info to DRAM results):
#conda activate DRAM
#DRAMDB=$(DRAM-setup.py print_config | grep "VOGDB db" | sed -e 's/.* \//\//g' -Ee 's/(\/.+)\/.*/\1/g')
#cd ${DRAMDB}
#rm -f vog.lca.tsv.gz; rm -f vog.lca.tsv;
#wget http://fileshare.csb.univie.ac.at/vog/latest/vog.lca.tsv.gz;
#gunzip vog.lca.tsv.gz;
#cut -f 4 vog.lca.tsv > tmp && mv tmp vog.lca.tsv
#zcat vog_annotations_latest.tsv.gz | paste - vog.lca.tsv | sed -e 's/|/\t\t/g' -e 's/GroupName/VOGDB_id/g' | sed -E -e 's/[^Aa-Zz0-9\t]/_/g' -e 's/_+/_/g' -e 's/^_|_$//g' -e 's/\t_|_\t/\t/g' | sed -e 's/\t\t/__/g' -e 's/\t/;/g' | awk -F ";" '{print "s/___"$1"___/"$0"/g"}' > vog_annotations_latest.lca.reformatted_for_replacement.txt
#conda deactivate

##### Virsorter2:
## Virsorter2 is accessed throughout script using conda activate vs2 - it must be installed locally via conda

####################### SET VARIABLES/ARGUMENTS #####################################

################################################################################################
####################### RUNNING THE SCRIPT #####################################################
################################################################################################
### STEP 001 must be finished before any other STEPS
## note that visualistaion is not run as part of a batch job because it is only run interactively (visualisation step)
## flag 1 must be $SAMPLES1 file (see details in STEP 001)
## flag 2 can be a number corresponding to which step will be run, e.g. A0 means run every step, A1 means run step 001 only
## sbatch /home/dcschroe/dmckeow/projects/DWV/script/anvio.sh /home/dcschroe/dmckeow/data/test_anvio/Albert_BC03_06_SAMPLES1.txt A0
## then do bash ${OUTPUTNAME1}_run_anvi_interative or copy paste the contents of this file into the terminal to run the interactive binning and visualisation

#################### TEST STEPS #################################
if [[ "$1" =~ "TEST" ]] || [[ "$1" =~ "TEST" ]]; then
  conda activate anvio-7.1
export TMPDIR=/scratch2/fr2424/sib/dmckeown/tmpdir
export TMP=/scratch2/fr2424/sib/dmckeown/tmp
### generate the contigs database
### anvi-gen-contigs-database does not have overwrite option (-W), so rm -f db before is necessary
rm -f $(basename $2 .fa).db
anvi-gen-contigs-database -f $2 -o $(basename $2 .fa).db -n "$(basename $2 .fa)" -T $SLURM_CPUS_PER_TASK ; \

fi




#################### STEP 001 #################################
if [[ "$2" =~ "A0" ]] || [[ "$2" =~ "A1" ]]; then
  SAMPLES1="$1"
  conda activate anvio-7.1

########## This step reformats fasta deflines to be compatible with anvi'o and concatenates your fastas
#### SAMPLES1 - ; delimited, 3 fields: contigs fasta of individual assembly (e.g. per barcode per sequencing run) [;] trimmed read fasta gzip for that individual assembly [;] NAME for subgrouping of co-assemblies OR separate assemblies to be concatenated
## e.g.:
##### ==> Albert_BC03_06_SAMPLES1.txt <==
## /home/dcschroe/dmckeow/projects/DWV/finalresult/assembly_Albert/04_24_2022_Canada_honeybee_barcode03.contigs.fasta;/home/dcschroe/dmckeow/projects/DWV/finalresult/assembly_Albert/04_24_2022_Canada_honeybee_barcode03.trimmedReads.fasta.gz;Albert_BC03_06
## /home/dcschroe/dmckeow/projects/DWV/finalresult/assembly_Albert/04_24_2022_Canada_honeybee_barcode04.contigs.fasta;/home/dcschroe/dmckeow/projects/DWV/finalresult/assembly_Albert/04_24_2022_Canada_honeybee_barcode04.trimmedReads.fasta.gz;Albert_BC03_06
## /home/dcschroe/dmckeow/projects/DWV/finalresult/assembly_Albert/04_24_2022_Canada_honeybee_barcode06.contigs.fasta;/home/dcschroe/dmckeow/projects/DWV/finalresult/assembly_Albert/04_24_2022_Canada_honeybee_barcode06.trimmedReads.fasta.gz;Albert_BC03_06

### fix fasta deflines and save "deflinekey" files to show corresponding original format and new format deflines
for SETS in $(cat $SAMPLES1); \
do \
    if [ "$SETS" == "SETS" ]; then continue; fi; \
    OUTPUTNAME=$(echo $SETS | cut -d ";" -f 1 | cut -d "." -f 1 | sed -E 's/.*\/(.+)/\1/g'); \
    CONTIGS=$(echo $SETS | cut -d ";" -f 1); \
    DEFLINEPREFIX=$(echo $SETS | cut -d ";" -f 1 | cut -d "." -f 1 | sed -E 's/.*\/(.+)/\1/g' | sed -E -e 's/[^A-Za-z0-9_]/_/g' -e 's/_+/_/g' -e 's/^_|_$//g'); \
    anvi-script-reformat-fasta $CONTIGS -o $(basename $CONTIGS) --simplify-names --report-file ${OUTPUTNAME}.deflinekey -l 2500 --seq-type NT --prefix "contig_${DEFLINEPREFIX}"; \

done; \

##### concatenate fasta files per GROUP
#################
for SETS in $(cat $SAMPLES1); \
do \
    if [ "$SETS" == "SETS" ]; then continue; fi; \
    CONTIGS=$(echo $SETS | cut -d ";" -f 1); \
    GROUP=$(echo $SETS | cut -d ";" -f 3); \
    rm -f ${GROUP}.fa; \
    touch ${GROUP}.fa; \
  done; \

for SETS in $(cat $SAMPLES1); \
do \
    if [ "$SETS" == "SETS" ]; then continue; fi; \
    CONTIGS=$(echo $SETS | cut -d ";" -f 1); \
    GROUP=$(echo $SETS | cut -d ";" -f 3); \
    cat $(basename $CONTIGS) >> ${GROUP}.fa; \
  done; \

#### create second running file for concatenated fastas, post defline cleanup
  ##### ==> Albert_BC03_06_SAMPLES2.txt <==
  ## /home/dcschroe/dmckeow/data/test_anvio/Albert_BC03_06.fa;/home/dcschroe/dmckeow/projects/DWV/finalresult/assembly_Albert/04_24_2022_Canada_honeybee_barcode03.trimmedReads.fasta.gz
  ## /home/dcschroe/dmckeow/data/test_anvio/Albert_BC03_06.fa;/home/dcschroe/dmckeow/projects/DWV/finalresult/assembly_Albert/04_24_2022_Canada_honeybee_barcode04.trimmedReads.fasta.gz
  ## /home/dcschroe/dmckeow/data/test_anvio/Albert_BC03_06.fa;/home/dcschroe/dmckeow/projects/DWV/finalresult/assembly_Albert/04_24_2022_Canada_honeybee_barcode06.trimmedReads.fasta.gz

awk -F ";" '{print $3".fa;"$1}' $SAMPLES1 | sed 's/;.*\//;/g' > SAMPLES2_$(basename $SAMPLES1)

fi


#################### STEP 002 #################################
if [[ "$2" =~ "A0" ]] || [[ "$2" =~ "A2" ]]; then
  SAMPLES1="$1"
  SAMPLES2="SAMPLES2_$(basename $SAMPLES1)"
  conda activate anvio-7.1
export TMPDIR=/scratch2/fr2424/sib/dmckeown/tmpdir
export TMP=/scratch2/fr2424/sib/dmckeown/tmp
### generate the contigs database
### anvi-gen-contigs-database does not have overwrite option (-W), so rm -f db before is necessary
for SETS in $(cat $SAMPLES2 | cut -d ";" -f 1 | sort -V | uniq); \
do \
    if [ "$SETS" == "SETS" ]; then continue; fi; \
    OUTPUTNAME=$(echo $SETS | cut -d "." -f 1 | sed -E 's/.*\/(.+)/\1/g'); \
    CONTIGS=$(echo $SETS); \
    rm -f ${OUTPUTNAME}.db; \
    anvi-gen-contigs-database -f $(basename $CONTIGS) -o ${OUTPUTNAME}.db -n "${OUTPUTNAME}" -T $SLURM_CPUS_PER_TASK ; \
done; \

fi

#################### STEP 003 #################################
if [[ "$2" =~ "A0" ]] || [[ "$2" =~ "A3" ]]; then
  SAMPLES1="$1"
  SAMPLES2="SAMPLES2_$(basename $SAMPLES1)"
  conda activate anvio-7.1


#### run hmm profile, summarise contig stats, and identify COGs
### !!! both anvi-run-hmms and anvi-run-ncbi-cogs do NOT include viruses and most eukaryotes
### !!! as an alternative to the default hmm profiles of anvio, we will later run and import the results of Virsorter2
### !!! as an alternative to the default COGs, we will import results from DRAM

for SETS in $(cat $SAMPLES2 | cut -d ";" -f 1 | sort -V | uniq); \
do \
    if [ "$SETS" == "SETS" ]; then continue; fi; \
    OUTPUTNAME=$(echo $SETS | cut -d "." -f 1 | sed -E 's/.*\/(.+)/\1/g'); \
    anvi-run-hmms -c ${OUTPUTNAME}.db -T $SLURM_CPUS_PER_TASK; \
    anvi-run-hmms -c ${OUTPUTNAME}.db -T $SLURM_CPUS_PER_TASK -H $VOGHMM; \
    anvi-run-hmms -c ${OUTPUTNAME}.db -T $SLURM_CPUS_PER_TASK -H $NCLDVHMM; \
    anvi-run-ncbi-cogs -c ${OUTPUTNAME}.db -T $SLURM_CPUS_PER_TASK; \
    anvi-run-kegg-kofams -c ${OUTPUTNAME}.db -T $SLURM_CPUS_PER_TASK; \
    anvi-display-contigs-stats ${OUTPUTNAME}.db --report-as-text -o ${OUTPUTNAME}.contigstats; \
done; \

### run single core gene taxonomy - this is built into anvio and cant be customised currently. It provides quick taxonomy for prokaryotes.
### !!! the alternative is to use Kraken2 and then import the results, which will provide gene level taxonomy
for SETS in $(cat $SAMPLES2 | cut -d ";" -f 1 | sort -V | uniq); \
do \
if [ "$SETS" == "SETS" ]; then continue; fi; \
  OUTPUTNAME=$(echo $SETS | cut -d "." -f 1 | sed -E 's/.*\/(.+)/\1/g'); \
  anvi-run-scg-taxonomy -c ${OUTPUTNAME}.db -T $SLURM_CPUS_PER_TASK; \
  anvi-estimate-scg-taxonomy -c ${OUTPUTNAME}.db --metagenome-mode -T $SLURM_CPUS_PER_TASK; \
done; \

fi

#################### STEP 004 #################################
if [[ "$2" =~ "A0" ]] || [[ "$2" =~ "A4" ]]; then
  SAMPLES1="$1"
  SAMPLES2="SAMPLES2_$(basename $SAMPLES1)"
  conda activate anvio-7.1

### get taxonomy for genecalls using Kraken2 and reformat:
for SETS in $(cat $SAMPLES2 | cut -d ";" -f 1 | sort -V | uniq); \
do \
if [ "$SETS" == "SETS" ]; then continue; fi; \
  OUTPUTNAME=$(echo $SETS | cut -d "." -f 1 | sed -E 's/.*\/(.+)/\1/g'); \
  anvi-get-sequences-for-gene-calls -c ${OUTPUTNAME}.db -o ${OUTPUTNAME}_gene_calls.fa; \
##### first run Kraken2 for microbial and some eukaryote databases
  kraken2 --db $KDB1 ${OUTPUTNAME}_gene_calls.fa \
    --use-mpa-style --threads $SLURM_CPUS_PER_TASK --use-names \
    --output ${OUTPUTNAME}.results.kraken2 \
    --report ${OUTPUTNAME}.report.kraken2 \
    --unclassified-out ${OUTPUTNAME}.unclassified.kraken2 \
    --classified-out ${OUTPUTNAME}.classified.kraken2; \

##### now run unclassified reads vs nt, to get hits for eukaryotes like animals, etc
  kraken2 --db $KDB2 ${OUTPUTNAME}.unclassified.kraken2 \
    --use-mpa-style --threads $SLURM_CPUS_PER_TASK --use-names \
    --output tmp_${OUTPUTNAME}.results.kraken2 \
    --report tmp_${OUTPUTNAME}.report.kraken2 \
    --unclassified-out tmp_${OUTPUTNAME}.unclassified.kraken2 \
    --classified-out tmp_${OUTPUTNAME}.classified.kraken2; \

##### merge the kraken results
sed -i '/unclassified (taxid 0)/d' ${OUTPUTNAME}.results.kraken2; \
cat ${OUTPUTNAME}.results.kraken2 tmp_${OUTPUTNAME}.results.kraken2 | sort -n -t $'\t' -k 2,2 > tmp_${OUTPUTNAME} && mv tmp_${OUTPUTNAME} ${OUTPUTNAME}.results.kraken2; \
cat ${OUTPUTNAME}.report.kraken2 tmp_${OUTPUTNAME}.report.kraken2 | sort -V > tmp_${OUTPUTNAME} && mv tmp_${OUTPUTNAME} ${OUTPUTNAME}.report.kraken2; \
mv tmp_${OUTPUTNAME}.unclassified.kraken2 ${OUTPUTNAME}.unclassified.kraken2
cat ${OUTPUTNAME}.classified.kraken2 tmp_${OUTPUTNAME}.classified.kraken2 | $SEQKIT sort -Nn - | $SEQKIT rmdup -n - > tmp_${OUTPUTNAME} && mv tmp_${OUTPUTNAME} ${OUTPUTNAME}.classified.kraken2; \
rm -f tmp_${OUTPUTNAME}.*kraken2

##### reformat kraken results

cut -f 2,3 ${OUTPUTNAME}.results.kraken2 | sed -E 's/(.*)\(taxid ([0-9]+)\)/\1\t\2/g' > tmp1_${OUTPUTNAME}; \

taxonkit reformat --threads $SLURM_CPUS_PER_TASK -F -I 3 tmp1_${OUTPUTNAME} | awk -F "\t" 'BEGIN{OFS="\t"};{if($2 ~ "unclassified") print $1,$2,$3,";;;;;;"; else print $0}' | cut -f 1,4 | sed 's/;/\t/g' | sed -z 's/^/gene_callers_id\tt_domain\tt_phylum\tt_class\tt_order\tt_family\tt_genus\tt_species\n/1' > ${OUTPUTNAME}_input_matrix.txt; \

anvi-import-taxonomy-for-genes -c ${OUTPUTNAME}.db -i ${OUTPUTNAME}_input_matrix.txt -p default_matrix; \
done; \

fi

#################### STEP 005 #################################
if [[ "$2" =~ "A0" ]] || [[ "$2" =~ "A5" ]]; then
  SAMPLES1="$1"
  SAMPLES2="SAMPLES2_$(basename $SAMPLES1)"

### run Virsorter2 and import results into Anvio2 (uses custom scripts made by VS2 developer - virsorter_to_anvio.py - anvi'o does not support this officially for VS2 yet)
for SETS in $(cat $SAMPLES2 | cut -d ";" -f 1 | sort -V | uniq); \
do \
if [ "$SETS" == "SETS" ]; then continue; fi; \
  OUTPUTNAME1=$(echo $SETS | cut -d "." -f 1 | sed -E 's/.*\/(.+)/\1/g'); \
  CONTIGS=$(echo $SETS); \
  conda activate anvio-7.1; \
  anvi-export-table --table splits_basic_info ${OUTPUTNAME1}.db -o ${OUTPUTNAME1}_splits_basic_info.txt; \
  anvi-export-gene-calls --gene-caller prodigal -c ${OUTPUTNAME1}.db -o ${OUTPUTNAME1}_all_gene_calls.txt; \
  conda activate vs2; \
  rm -fr ${OUTPUTNAME1}_VS2_results; \
  virsorter run -i $(basename $CONTIGS) -w ${OUTPUTNAME1}_VS2_results/ --keep-original-seq --prep-for-dramv --hallmark-required-on-short -j $SLURM_CPUS_PER_TASK --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae; \
  $VS2_ANVIO -i ${OUTPUTNAME1}_VS2_results/ -s ${OUTPUTNAME1}_splits_basic_info.txt -n ${OUTPUTNAME1}_all_gene_calls.txt -d ~/vs2/db/ -A ${OUTPUTNAME1}_virsorter_additional_info.txt -C ${OUTPUTNAME1}_virsorter_collection.txt -F ${OUTPUTNAME1}_virsorter_annotations.txt; \
done; \

fi

#################### STEP 006 #################################
if [[ "$2" =~ "A0" ]] || [[ "$2" =~ "A6" ]]; then
  SAMPLES1="$1"
  SAMPLES2="SAMPLES2_$(basename $SAMPLES1)"

### run DRAM and import results into Anvi'o
## if you get DRAM error about database codes, try: DRAM-setup.py update_description_db
for SETS in $(cat $SAMPLES2 | cut -d ";" -f 1 | sort -V | uniq); \
do \
if [ "$SETS" == "SETS" ]; then continue; fi; \
  OUTPUTNAME=$(echo $SETS | cut -d "." -f 1 | sed -E 's/.*\/(.+)/\1/g'); \
  conda activate anvio-7.1; \
  anvi-get-sequences-for-gene-calls --report-extended-deflines --get-aa-sequences -c ${OUTPUTNAME}.db -o ${OUTPUTNAME}_gene_calls.faa; \
  sed -i -e 's/;/;;;/1' -Ee 's/>([0-9]+) contig:(.+);;;.*/>\2__gene__\1/g' -e 's/>.*_gene_calls_/>/g' ${OUTPUTNAME}_gene_calls.faa; \
  conda activate DRAM; \
  DRAMDB=$(DRAM-setup.py print_config | grep "VOGDB db" | sed -e 's/.* \//\//g' -Ee 's/(\/.+)\/.*/\1/g'); \
  rm -fr ${OUTPUTNAME}_DRAM_results; \
  DRAM.py annotate_genes --threads $SLURM_CPUS_PER_TASK -i ${OUTPUTNAME}_gene_calls.faa -o ${OUTPUTNAME}_DRAM_results --custom_hmm_name VOGDB --custom_hmm_loc ${DRAMDB}/vog_latest_hmms.txt; \

  mkdir ${OUTPUTNAME}_transpose
  split -d -l 1 ${OUTPUTNAME}_DRAM_results/annotations.tsv ${OUTPUTNAME}_transpose/tmp1_
  sed -i 's/\t/\n/g' ${OUTPUTNAME}_transpose/tmp1_*
  find ${OUTPUTNAME}_transpose/ -maxdepth 1 -type f -name tmp1_\* | sort -V | split -l 1000 -d - ${OUTPUTNAME}_transpose/tmp_lists_
  for f in ${OUTPUTNAME}_transpose/tmp_lists_*; do paste $(cat $f) > ${OUTPUTNAME}_transpose/tmp_merge_$(basename $f); done
  paste ${OUTPUTNAME}_transpose/tmp_merge_* > ${OUTPUTNAME}_DRAM_results/tmp_annotations.txt

  grep "VOGDB_id" ${OUTPUTNAME}_DRAM_results/tmp_annotations.txt | sed 's/\t/\n/g' | sed -Ee 's/(.*)/___\1___/g' -e 's/______//g' > ${OUTPUTNAME}_DRAM_results/tmp_vogdb; \
  sed -i -f ${DRAMDB}/vog_annotations_latest.lca.reformatted_for_replacement.txt ${OUTPUTNAME}_DRAM_results/tmp_vogdb; \
  sed -i -e '1,1d' -e 's/;/__/g' -e 's/$/\t0/g' ${OUTPUTNAME}_DRAM_results/tmp_vogdb; \
  sed -i -z 's/^/function\te_value\n/1' ${OUTPUTNAME}_DRAM_results/tmp_vogdb; \
  awk -F "\t" '{print $1"\tVOGDB\t"$4}' ${OUTPUTNAME}_DRAM_results/annotations.tsv | sed -e '1,1d' -Ee 's/^.+__gene__([0-9]+\t)/\1/g' -e 's/;|; /__/g' | sed -z 's/^/gene_callers_id\tsource\taccession\n/1' | paste - ${OUTPUTNAME}_DRAM_results/tmp_vogdb | awk -F "\t" '{if($4 ~ /VOG|function/) print $0}' > ${OUTPUTNAME}_DRAM_annotations.txt; \

  grep "pfam_hits" ${OUTPUTNAME}_DRAM_results/tmp_annotations.txt | sed 's/\t/\n/g' > ${OUTPUTNAME}_DRAM_results/tmp_pfam; \
  sed -i -E -e 's/[^Aa-Zz0-9;]/_/g' -e 's/_+/_/g' -e 's/^_|_$//g' -e 's/_;_/__/g' ${OUTPUTNAME}_DRAM_results/tmp_pfam
  sed -i -e '1,1d' -e 's/;/__/g' -e 's/$/\t0/g' ${OUTPUTNAME}_DRAM_results/tmp_pfam; \
  awk -F "\t" '{print $1"\tPFAM\t"$4}' ${OUTPUTNAME}_DRAM_results/annotations.tsv | sed -e '1,1d' -Ee 's/^.+__gene__([0-9]+\t)/\1/g' -e 's/;|; /__/g' | paste - ${OUTPUTNAME}_DRAM_results/tmp_pfam | awk -F "\t" '{if($4 ~ /[0-9]|function/) print $0}' >> ${OUTPUTNAME}_DRAM_annotations.txt; \


done; \

fi


#################### STEP 007 #################################
if [[ "$2" =~ "A0" ]] || [[ "$2" =~ "A7" ]]; then
  SAMPLES1="$1"
  SAMPLES2="SAMPLES2_$(basename $SAMPLES1)"
  conda activate anvio-7.1

#### make BAM file of mapping reads vs contigs and sort and index it
### the same concatenate contig fasta must be used for all read files

for SETS in $(cat $SAMPLES1); \
do \
    if [ "$SETS" == "SETS" ]; then continue; fi; \
    READS=$(echo $SETS | cut -d ";" -f 2); \
    CONTIGS=$(echo $SETS | cut -d ";" -f 3); \
    OUTPUTNAME1=$(echo $SETS | cut -d ";" -f 2 | cut -d "." -f 1 | sed -E 's/.*\/(.+)/\1/g'); \
    minimap2 -t $SLURM_CPUS_PER_TASK -ax map-ont ${CONTIGS}.fa $READS | samtools sort -O BAM - > ${OUTPUTNAME1}.bam; \
    samtools index -@ $SLURM_CPUS_PER_TASK -b ${OUTPUTNAME1}.bam; \
done; \

fi

#################### STEP 008 #################################
#### profile the contig databases WITH CONCATENATED FASTA OR CO_ASSEMBLIES ONLY - MULTIPLE SEQUENCING SAMPLES (SET OF READS) AND ONE CONTIG FASTA
if [[ "$2" =~ "A0" ]] || [[ "$2" =~ "A8" ]]; then
  ##SAMPLES1="$1"
  ##SAMPLES2="SAMPLES2_$(basename $SAMPLES1)"
  conda activate anvio-7.1

##for SETS in $(cat $SAMPLES2); \
##do \
    ##if [ "$SETS" == "SETS" ]; then continue; fi; \
    ##OUTPUTNAME1=$(echo $SETS | cut -d ";" -f 2 | cut -d "." -f 1 | sed -E 's/.*\/(.+)/\1/g'); \
    ##OUTPUTNAME2=$(echo $SETS | cut -d ";" -f 1 | cut -d "." -f 1 | sed -E 's/.*\/(.+)/\1/g'); \
    N1=$(echo $3 | cut -d "." -f 1 | sed -E 's/.*\/(.+)/\1/g' | sed 's/-/_/g')
    samtools index -@ $SLURM_CPUS_PER_TASK -b $3
    anvi-profile -i $3 -c ${4} --sample-name "profile_${N1}" --min-coverage-for-variability 10 -T $SLURM_CPUS_PER_TASK -W
##done; \

fi

#################### STEP 009 #################################
if [[ "$2" =~ "A0" ]] || [[ "$2" =~ "A9" ]]; then
  SAMPLES1="$1"
  SAMPLES2="SAMPLES2_$(basename $SAMPLES1)"
  SAMPLES3="SAMPLES3_$(basename $SAMPLES1)"

  conda activate anvio-7.1

#### IMPORT Virsorter2 AND DRAM results (before merging)
if [ -f "${OUTPUTNAME2}_virsorter_annotations.txt" ]; then
for SETS in $(cat $SAMPLES2); \
do \
    OUTPUTNAME1=$(echo $SETS | cut -d ";" -f 2 | cut -d "." -f 1 | sed -E 's/.*\/(.+)/\1/g'); \
    OUTPUTNAME2=$(echo $SETS | cut -d ";" -f 1 | cut -d "." -f 1 | sed -E 's/.*\/(.+)/\1/g'); \
    anvi-import-functions -c ${OUTPUTNAME2}.db -i ${OUTPUTNAME2}_virsorter_annotations.txt; \
done; \
fi

if [ -f "${OUTPUTNAME2}_DRAM_annotations.txt" ]; then
for SETS in $(cat $SAMPLES2); \
do \
    OUTPUTNAME1=$(echo $SETS | cut -d ";" -f 2 | cut -d "." -f 1 | sed -E 's/.*\/(.+)/\1/g'); \
    OUTPUTNAME2=$(echo $SETS | cut -d ";" -f 1 | cut -d "." -f 1 | sed -E 's/.*\/(.+)/\1/g'); \
    anvi-import-functions -c ${OUTPUTNAME2}.db -i ${OUTPUTNAME2}_DRAM_annotations.txt; \
done; \
fi


#### merge the database profiles
awk -F ";" -v S="$SAMPLES3" 'BEGIN{OFS=";"};{print $2 > "tmp_"S"_"$1".1"}' $SAMPLES2; \
for f in tmp_${SAMPLES3}_*.1; do sed -E 's/(.+)\.contigs\.fasta/profile_\1\/PROFILE.db /g' $f | tr -d '\n' | sed 's/$/\n/g' > $(basename $f .1).2 ; done; \
for f in tmp_${SAMPLES3}_*.2; do awk -v f=$f '{print f";"$0}' $f | sed -E "s/tmp_${SAMPLES3}_(.+)\.fa\.2;/\1;/g" > $(basename $f .2).3; done; \
cat tmp_${SAMPLES3}_*.3 | sed 's/ /__/g' > $SAMPLES3; \
rm -f tmp_${SAMPLES3}_*; \

for SETS in $(cat $SAMPLES3); \
do \
MERGENAME=$(echo $SETS | cut -d ";" -f 1); \
PROFILELIST=$(echo $SETS | cut -d ";" -f 2 | sed 's/__/ /g'); \
anvi-merge ${PROFILELIST} -o profile_${MERGENAME}-MERGED -c ${MERGENAME}.db -W; \

#### post merge Virsorter2 import
if [ -f "${MERGENAME}_virsorter_collection.txt" ]; then
anvi-import-collection ${MERGENAME}_virsorter_collection.txt -c ${MERGENAME}.db -p profile_${MERGENAME}-MERGED/PROFILE.db -C ${MERGENAME}_VirSorter2; \
fi

if [ -f "${MERGENAME}_virsorter_additional_info.txt" ]; then
anvi-import-misc-data ${MERGENAME}_virsorter_additional_info.txt -p profile_${MERGENAME}-MERGED/PROFILE.db --target-data-table items; \
fi

#### prepare command for interactive interface and visualise
### NOT for job submission or batch; interactive only
## this generates a file that can be copy pasted from or use bash command on to run visualisation step in terminal
echo "conda activate anvio-7.1; anvi-interactive -p profile_${MERGENAME}-MERGED/PROFILE.db -c ${MERGENAME}.db --browser-path /usr/bin/chromium-browser" > ${MERGENAME}_run_anvi_interactive; \
done; \
fi
