#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=600GB
#SBATCH --partition ag2tb
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err


Help()
{
echo -e "########### REQUIRED arguments ###########\n-i, --input\tfull path to file listing input file locations, semi-colon delimited, 4 fields:\n/path/barcode01.contigs.fasta;/path/barcode01.fastq.gz;group_location1;BC01"
echo -e "\t\tFIELD1=contig fasta, separate line per assembly (full path). Can be individual assemblies or co-assemblies."
echo -e "\t\tFIELD2=reads fq.gz, corresponding to each assembly contig fasta (full path)"
echo -e "\t\tFIELD3=group name; samples sharing the same group name will be concatenated together will appear together in the same anvio figures. To prevent an assembly from being concatenated with any other (for example, if you had multilpe co-assemblies for which you wanted a separate set of results and figures), simply give it a unique group name."
echo -e "\t\tFIELD4=profile name; a short name for each profile which will appear in all outputs and figures."
echo -e "!!! group and profile names MUST NOT begin with a digit, and MUST ONLY include letters, digits, and underscores !!!\n"
echo -e "-s, --step\tWhich script steps to run. A0 = run the whole script. A1-A9, B1-B4 for specific steps"
echo -e "-m, --mincontigsize\tinteger to set minimum size of contigs in bp to include in analyses (don't go below 1000)"
echo -e "-o, --output\tfull path for output directory\n"
echo -e "########### OPTIONAL arguments ###########\n-h, --help\tshow this help message and exit\n"
echo -e "-r, --rerun\tto rerun certain steps without starting all over again; A6"
echo -e "-L, --layers\ta tsv file to add additional info to the layers, visible in the layer plots, step B4. Include info for each group of samples that you want to add extra info for:\n\t\tsamples density varroa_per_100_bees\n\t\tSU_TX_22_6083   0.0000000935    1.333333333\n\t\tSU_TX_22_6084   0.0000000935    6.333333333\n\t\tEQ_TX_22_5997   0.000000194     0\n!!! FIELD1 - "samples" must match exactly your profile names !!!\n"
echo -e "-O, --order\ta tsv file to add additional info to re-order the layers by, step B4. For each variable, you must include a comma-separated list of samples sorted in the desired order:\n\t\titem_name   data_type   data_value\n\t\tphylogeny   newick  (sample2:0.0370199,(sample1:0.0227268,sample3:0.0227268)Int3:0.0370199)\n\t\tdensity   basic   sample3,sample2,sample1\n!!! FIELD1 - "samples" must match exactly your profile names !!!\n"
echo -e "-h, --help\tshow this help message and exit\n"
}

VOGDBmismatch()
{
echo -e "VOGDB versions used to build databases for DRAMDB versus anvio/diamond are NOT the same; rebuild them!"
}

while getopts i:s:m:o:rhL:O: option
do 
    case "${option}" in
        i)input=${OPTARG};;
        s)step=${OPTARG};;
    m)mincontigsize=${OPTARG};;
    o)output=${OPTARG};;
    r)rerun=${OPTARG};;
    h)Help; exit;;
    L)layers=${OPTARG};;
    O)order=${OPTARG};;
    esac
done

if [[ -z "${input}" ]]; then echo "-i, --input REQUIRED"; Help; exit; fi
if [[ -z "${step}" ]]; then echo "-s, --step REQUIRED"; Help; exit; fi
if [[ -z "${mincontigsize}" ]]; then echo "-m, --mincontigsize REQUIRED"; Help; exit; fi
if [[ -z "${output}" ]]; then echo "-o, --output REQUIRED"; Help; exit; fi




####################### PREREQUISITES #####################################
## nanopore seq reads data

####################### SCRIPT PURPOSE #####################################

####################### HOW TO RUN SCRIPT #####################################

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################

## if you intend to merge all samples together for analysis and visualisation in avi'o, then you must use either a co-assembly (multilpe sequencing runs and one contig fasta assembly) OR a concatenated fasta of independently assembled contigs

####################### OUTPUT FILES #####################################

####################### SOFTWARE, DATABSES, ETC #####################################
### LOAD software available via shared environment on server:
module purge
eval "$(conda shell.bash hook)"

### set paths for DATABASES that must be installed locally:
###### KRAKEN2
touch /scratch.global/dcschroe/kraken2db_nt/*
KDBNT="/scratch.global/dcschroe/kraken2db_nt"; ### a DB comprised of NCBI's nt: see Kraken2 database bulding for instructions. NOTE - this will be very large (~> 1 TB), so do it in scratch
KDBNT_s3="s3://dmckeow/kraken2db_nt/" ### location for krakendb to be stored on 2nd tier storage

NCLDVHMM="/home/dcschroe/dmckeow/HMMs_NCLDVs/HMM_NCLDV_149_core_genes" ### a hmmer profile made by the anvio team for NCLDVs - available online
VS2DB="/home/dcschroe/dmckeow/vs2/db"


VOGDB="/home/dcschroe/dmckeow/VOGDB"
DRAMDB="/home/dcschroe/dmckeow/DRAMDB/"
DIAMONDDB="/home/dcschroe/dmckeow/VOGDB/DBVOG.dmnd"
VOGHMM="$VOGDB/anvioHMMs"; ### the VOG db hmmer profile, downloaded and then additional files created to match the anvio requirements for external hmmer profiles
VOGHMMcutoff="1e-30" #### this is the evalue cut off for VOG HMM scan with anvio

########################################################################
########################################################################

######### SOFTWARE/OTHER SCRIPTS THAT NEED SETUP BEFORE RUNNING:
module load pigz

VS2_ANVIO="/home/dcschroe/dmckeow/VirSorter2_to_Anvio/virsorter_to_anvio.py" ### a custom script made by the VS2 developers to reformat VS2 results for import into anvio
COVERAGE2RPKM="/home/dcschroe/dmckeow/projects/DWV/script/coverage_to_RPKM.py" ## custom script that converts coverage data into RPKM

##### Seqkit - must be installed locally, e.g.:
SEQKIT="/home/dcschroe/dmckeow/seqkit"

##### KRAKEN2 - must be installed locally, and be in your $PATH

##### ANVIO 7.1 or later - accessed using conda activate anvio-7.1 - it must be installed locally via conda

##### Virsorter2 - accessed using conda activate vs2 - it must be installed locally via conda - see instructions to setup VS2 database

##### Diamond - installed locally
DIAMOND="/home/dcschroe/dmckeow/diamond"

##### datamash - installed locally
DATAMASH="~/bin/bin/datamash"

######################################################################
######################################################################

########### KRAKEN2 ####################
########## TO BUILD CUSTOM DATABASE, KRAKEN SCRIPTS MUST BE CHANGED AS FOLLOWS:
#in download_genomic_library.sh (line 17)
#replace
#FTP_SERVER="ftp://$NCBI_SERVER"
#with
#FTP_SERVER="https://$NCBI_SERVER"

#in rsync_from_ncbi.pl (line 46)
#replace
#if (! ($full_path =~ s#^ftp://${qm_server}${qm_server_path}/##)) {
#with
#if (! ($full_path =~ s#^https://${qm_server}${qm_server_path}/##)) {

if [[ "$step" == "krakenbuild" ]]; then

### get taxonomy
kraken2-build --threads $SLURM_CPUS_PER_TASK --download-taxonomy --db $KDBNT
### get specific databases (one line per db)
kraken2-build --threads $SLURM_CPUS_PER_TASK --download-library nt --db $KDBNT
## make the finalised db
kraken2-build --build --fast-build --threads $SLURM_CPUS_PER_TASK --db $KDBNT

s3cmd rm -fr $KDBNT_s3
s3cmd put -r $KDBNT_s3
fi

############## to re-download KRAKEN2 DATABASE FROM 2ND TIER STORAGE TO SCRATCH (FASTER THAN REBUILDING IT)
if [[ "$step" == "krakenrenew" ]]; then

s3cmd get -fr $KDBNT_s3 $KDBNT

fi


######################################################################
######################################################################
###### VIRSORTER2 DATABASE SETUP
if [[ "$step" == "vs2build" ]]; then

conda activate vs2
virsorter setup -d $VS2DB -j $SLURM_CPUS_PER_TASK
conda deactivate

fi

######################################################################
######################################################################

############## SETUP DBs for VOGDB, DRAM VOGDB, and VOGDB for Diamond, AND VOG HMMs for Anvi'o
#### VOGDB and DRAM MUST be SETUP AT THE SAME TIME TO ensure the same version of VOGDB is used.
##### ALL of these databases MUST be built again from scratch if you upgrade the saved version of VOGDB

if [[ "$step" == "VOGbuild" ]]; then

rm -fr $VOGDB; mkdir $VOGDB; cd $VOGDB

##### get latest VOGDB - a database of curated core genes for a wide range of virus groups
## also creates a file that will add more all useful info to VOGs

wget --cut-dirs 2 -r -np -nH -R "index.html*" https://fileshare.csb.univie.ac.at/vog/latest/
zcat vog.lca.tsv | cut -f 4 > tmp_vog.lca.tsv
zcat vog.annotations.tsv.gz | paste - tmp_vog.lca.tsv | sed -e 's/|/\t\t/g' -e 's/GroupName/VOGDB_id/g' | sed -E -e 's/[^Aa-Zz0-9\t]/_/g' -e 's/_+/_/g' -e 's/^_|_$//g' -e 's/\t_|_\t/\t/g' | sed -e 's/\t\t/__/g' -e 's/\t/;/g' | awk -F ";" '{print "s/___"$1"___/"$0"/g"}' > vog_annotations.lca.reformatted_for_replacement.txt
rm -f tmp_vog.lca.tsv

##### DRAM - prepare DRAM databases, including VOGDB:
conda activate DRAM
DRAM-setup.py prepare_databases --skip_uniref --output_dir $DRAMDB
conda deactivate


#### check if DRAM and VOGDB are the same version


diff1=$(diff $DRAMDB/vog_annotations_latest.tsv.gz $VOGDB/vog.annotations.tsv.gz | wc -l)

if [[ "$diff1" == '1' ]]; then

VOGDBmismatch
exit

fi




########## DIAMOND - Setup Diamond database for VOG
gunzip vog.members.tsv.gz
rm -fr tmp_VOG_faa; mkdir tmp_VOG_faa; tar --directory tmp_VOG_faa -xzvf vog.faa.tar.gz

for f in tmp_VOG_faa/VOG*.faa; do
  V=$(echo $f | sed 's/\.faa/./g'); sed -i "s/>/>$V/g" $f
done

cat tmp_VOG_faa/VOG*.faa > DBVOG.faa
$DIAMOND makedb --threads $SLURM_CPUS_PER_TASK --in DBVOG.faa -d $VOGDB
rm -f DBVOG.faa; rm -fr tmp_VOG_faa
gzip vog.members.tsv

############ SETUP HMMs for VOGDB for Anvi'o

rm -fr $VOGHMM; mkdir $VOGHMM; cd $VOGHMM

cp $DRAMDB/vog_latest_hmms.txt genes.hmm
##pigz -p $SLURM_CPUS_PER_TASK -dc $VOGDB/vog.hmm.tar.gz | tar --directory . -xf -
pigz -p $SLURM_CPUS_PER_TASK genes.hmm

zcat genes.hmm.gz | awk '/^NAME *VOG[0-9]+/' | sed 's/^NAME *//g' | sort -Vu | awk '{print $0"\tnone""\tVOGDB"}' | sed -z 's/^/gene\taccession\thmmsource\n/1' > genes.txt

echo "Virus_orthologous_groups" > kind.txt

echo -e "-E $VOGHMMcutoff" > noise_cutoff_terms.txt

echo "https://vogdb.org/" > reference.txt

echo -e "AA:GENE" > target.txt

fi

####################### THE SCRIPT #####################################

################################################################################################
################################################################################################
########################################################################
########################################################################

#### GO TO OUTPUT DIRECTORY
cd $output

#################### STEP A1 #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A1" ]]; then
  conda activate anvio-7.1

### fix fasta deflines and save "deflinekey" files to show corresponding original format and new format deflines
for SETS in $(cat $input | sort -Vu -t ";" -k 1,1)
do
  CONTIGNAME=$(echo $SETS | cut -d ";" -f 1 | cut -d "." -f 1)
  CONTIGS=$(echo $SETS | cut -d ";" -f 1)
  PROFILENAME=$(echo $SETS | cut -d ";" -f 4)
  anvi-script-reformat-fasta $CONTIGS -o $(basename $CONTIGS) --simplify-names --report-file $(basename $CONTIGNAME).deflinekey -l $mincontigsize --seq-type NT --prefix "contig_${PROFILENAME}"
done

##### concatenate fasta files per GROUP1
#################
for SETS in $(cat $input | sort -Vu -t ";" -k 1,1)
do
    CONTIGS=$(echo $SETS | cut -d ";" -f 1)
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    rm -f ${GROUP1}.fa
    touch ${GROUP1}.fa
    rm -f ${GROUP1}.deflinekey
    touch ${GROUP1}.deflinekey
  done

for SETS in $(cat $input | sort -Vu -t ";" -k 1,1)
do
    CONTIGS=$(echo $SETS | cut -d ";" -f 1)
    CONTIGNAME=$(echo $SETS | cut -d ";" -f 1 | cut -d "." -f 1)
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    cat $(basename $CONTIGS) >> ${GROUP1}.fa
    cat $(basename $CONTIGNAME).deflinekey >> ${GROUP1}.deflinekey
    rm -f $(basename $CONTIGS)
    rm -f $(basename $CONTIGNAME).deflinekey
  done

fi


#################### STEP A2 #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A2" ]]; then
  conda activate anvio-7.1

### generate the contigs database
### anvi-gen-contigs-database does not have overwrite option (-W), so rm -f db before is necessary
for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    rm -f ${GROUP1}.db
    anvi-gen-contigs-database -f ${GROUP1}.fa -o ${GROUP1}.db -n "${GROUP1}" -T $SLURM_CPUS_PER_TASK
done

fi

#################### STEP A3 #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A3" ]]; then
  conda activate anvio-7.1


#### run hmm profile, summarise contig stats, and identify COGs
### !!! both anvi-run-hmms and anvi-run-ncbi-cogs do NOT include viruses and most eukaryotes
### !!! as an alternative to the default hmm profiles of anvio, we will later run and import the results of Virsorter2
### !!! as an alternative to the default COGs, we will import results from DRAM

for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    anvi-run-hmms -c ${GROUP1}.db -T $SLURM_CPUS_PER_TASK --just-do-it
    anvi-run-hmms -c ${GROUP1}.db -T $SLURM_CPUS_PER_TASK -H $VOGHMM --hmmer-program hmmsearch --just-do-it
    anvi-run-hmms -c ${GROUP1}.db -T $SLURM_CPUS_PER_TASK -H $NCLDVHMM --hmmer-program hmmsearch --just-do-it
    anvi-run-ncbi-cogs -c ${GROUP1}.db -T $SLURM_CPUS_PER_TASK --just-do-it
    anvi-run-kegg-kofams -c ${GROUP1}.db -T $SLURM_CPUS_PER_TASK --just-do-it
    anvi-display-contigs-stats ${GROUP1}.db --report-as-text -o ${GROUP1}.contigstats
done

### run single core gene taxonomy - this is built into anvio and cant be customised currently. It provides quick taxonomy for prokaryotes.
### !!! the alternative is to use Kraken2 and then import the results, which will provide gene level taxonomy
for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
  GROUP1=$(echo $SETS | cut -d ";" -f 3)
  anvi-run-scg-taxonomy -c ${GROUP1}.db -T $SLURM_CPUS_PER_TASK
  anvi-estimate-scg-taxonomy -c ${GROUP1}.db --metagenome-mode -T $SLURM_CPUS_PER_TASK
done

fi

#################### STEP A4 #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A4" ]]; then
  conda activate anvio-7.1

### get taxonomy for genecalls using Kraken2 and reformat:
for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
  GROUP1=$(echo $SETS | cut -d ";" -f 3)
  anvi-get-sequences-for-gene-calls -c ${GROUP1}.db -o ${GROUP1}_gene_calls.fa
##### first run Kraken2 for microbial and some eukaryote databases
  kraken2 --db $KDBNT ${GROUP1}_gene_calls.fa --use-mpa-style --threads $SLURM_CPUS_PER_TASK --use-names --output ${GROUP1}.results.kraken2 --report ${GROUP1}.report.kraken2 --unclassified-out ${GROUP1}.unclassified.kraken2 --classified-out ${GROUP1}.classified.kraken2

##### reformat kraken results for anvio import
cut -f 2,3 ${GROUP1}.results.kraken2 | sed -E 's/(.*)\(taxid ([0-9]+)\)/\1\t\2/g' > tmp1_${GROUP1}

taxonkit reformat --threads $SLURM_CPUS_PER_TASK -F -I 3 tmp1_${GROUP1} | awk -F "\t" 'BEGIN{OFS="\t"};{if($2 ~ "unclassified" || $4 =="") print $1,$2,$3,";;;;;;"; else print $0}' | cut -f 1,4 | sed 's/;/\t/g' | sed -z 's/^/gene_callers_id\tt_domain\tt_phylum\tt_class\tt_order\tt_family\tt_genus\tt_species\n/1' > ${GROUP1}_input_matrix.txt

anvi-import-taxonomy-for-genes -c ${GROUP1}.db -i ${GROUP1}_input_matrix.txt -p default_matrix
done

fi

#################### STEP A5 #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A5" ]]; then

### run Virsorter2 and import results into Anvio2 (uses custom scripts made by VS2 developer - virsorter_to_anvio.py - anvi'o does not support this officially for VS2 yet)
for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
  GROUP1=$(echo $SETS | cut -d ";" -f 3)
  conda activate anvio-7.1
  anvi-export-table --table splits_basic_info ${GROUP1}.db -o ${GROUP1}_splits_basic_info.txt
  anvi-export-gene-calls --gene-caller prodigal -c ${GROUP1}.db -o ${GROUP1}_all_gene_calls.txt
  conda activate vs2
  rm -fr ${GROUP1}_VS2_results
  virsorter run -i ${GROUP1}.fa -w ${GROUP1}_VS2_results/ --keep-original-seq --prep-for-dramv --hallmark-required-on-short --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae -j $SLURM_CPUS_PER_TASK
  $VS2_ANVIO -i ${GROUP1}_VS2_results/ -s ${GROUP1}_splits_basic_info.txt -n ${GROUP1}_all_gene_calls.txt -d ${VS2DB} -A ${GROUP1}_virsorter_additional_info.txt -C ${GROUP1}_virsorter_collection.txt -F ${GROUP1}_virsorter_annotations.txt
done

fi


#################### STEP A6 #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A6" ]] && [[ -z "${rerun}" ]]; then


### run DRAM and import results into Anvi'o
## if you get DRAM error about database codes, try: DRAM-setup.py update_description_db
for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
  GROUP1=$(echo $SETS | cut -d ";" -f 3)
  conda activate anvio-7.1
  anvi-get-sequences-for-gene-calls --report-extended-deflines --get-aa-sequences -c ${GROUP1}.db -o ${GROUP1}_gene_calls.faa
  sed -i -e 's/;/;;;/1' -Ee 's/>([0-9]+) contig:(.+);;;.*/>\2__gene__\1/g' -e 's/>.*_gene_calls_/>/g' ${GROUP1}_gene_calls.faa

### split the AA gene calls
  ##rm -fr ${GROUP1}_gene_calls.part_*.faa
  ##$SEQKIT split ${GROUP1}_gene_calls.faa -s 1000000 -O .

rm -fr ${GROUP1}_DRAM_results

  conda activate DRAM

srun --ntasks=1 --cpus-per-task=6 DRAM.py annotate_genes --threads $SLURM_CPUS_PER_TASK -i ${GROUP1}_gene_calls.faa -o ${GROUP1}_DRAM_results --custom_hmm_name VOGDB --custom_hmm_loc ${DRAMDB}/vog_latest_hmms.txt

conda deactivate

#############
########## rerun flag skips first part and restarts the DRAM step ONLY on
############ NEEDS FIXED - DRAM merge does not work - the transpose part next could be changed to allow merging of multiple DRAM runs with different fields
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A6" ]] && [[ ! -z "${rerun}" ]]; then

find ${GROUP1}_gene_calls.faa.split -name "annotations.tsv" | sed 's/_DRAM_results\/annotations.tsv//g' > tmp_DRAM_${GROUP1}_done;
ls ${GROUP1}_gene_calls.faa.split/*.faa | grep -vw -f tmp_DRAM_${GROUP1}_done - > tmp_DRAM_${GROUP1}_NOTdone


for f in $(cat tmp_DRAM_${GROUP1}_NOTdone); do
  rm -fr ${f}_DRAM_results
  DRAM.py annotate_genes --threads $SLURM_CPUS_PER_TASK -i $f -o ${f}_DRAM_results --custom_hmm_name VOGDB --custom_hmm_loc ${DRAMDB}/vog_latest_hmms.txt
done

#### merge all DRAM results together
rm -fr ${GROUP1}_DRAM_results
DRAM.py merge_annotations -i ${GROUP1}_gene_calls.faa.split/*_DRAM_results -o ${GROUP1}_DRAM_results




fi

#########################

  $DATAMASH transpose < ${GROUP1}_DRAM_results/annotations.tsv > ${GROUP1}_DRAM_results/tmp_annotations.txt

### rm -f ${GROUP1}_DRAM_results/tmp_annotation_headers; touch ${GROUP1}_DRAM_results/tmp_annotation_headers
### cut -f 1 ${GROUP1}_DRAM_results/tmp_annotations.txt | sort -Vu >> ${GROUP1}_DRAM_results/tmp_annotation_headers
### grep -w -f ${GROUP1}_DRAM_results/tmp_annotation_headers ${GROUP1}_DRAM_results/tmp_annotations.txt > ?_headername
  grep "VOGDB_id" ${GROUP1}_DRAM_results/tmp_annotations.txt | sed 's/\t/\n/g' | sed -Ee 's/(.*)/___\1___/g' -e 's/______//g' > ${GROUP1}_DRAM_results/tmp_vogdb
  sed -i -f ${VOGDB}/vog_annotations.lca.reformatted_for_replacement.txt ${GROUP1}_DRAM_results/tmp_vogdb
  sed -i -e '1,1d' -e 's/;/__/g' -e 's/$/\t0/g' ${GROUP1}_DRAM_results/tmp_vogdb
  sed -i -z 's/^/function\te_value\n/1' ${GROUP1}_DRAM_results/tmp_vogdb
  ##### IF ${GROUP1}_DRAM_results/tmp_vogdb
  awk -F "\t" '{print $1"\tVOGDB\t"$4}' ${GROUP1}_DRAM_results/annotations.tsv | sed -e '1,1d' -Ee 's/^.+__gene__([0-9]+\t)/\1/g' -e 's/;|; /__/g' | sed -z 's/^/gene_callers_id\tsource\taccession\n/1' | paste - ${GROUP1}_DRAM_results/tmp_vogdb | awk -F "\t" '{if($4 ~ /VOG|function/) print $0}' > ${GROUP1}_DRAM_annotations.txt

  grep "pfam_hits" ${GROUP1}_DRAM_results/tmp_annotations.txt | sed 's/\t/\n/g' > ${GROUP1}_DRAM_results/tmp_pfam
  sed -i -E -e 's/[^Aa-Zz0-9;]/_/g' -e 's/_+/_/g' -e 's/^_|_$//g' -e 's/_;_/__/g' ${GROUP1}_DRAM_results/tmp_pfam
  sed -i -e '1,1d' -e 's/;/__/g' -e 's/$/\t0/g' ${GROUP1}_DRAM_results/tmp_pfam
  awk -F "\t" '{print $1"\tPFAM\t"$4}' ${GROUP1}_DRAM_results/annotations.tsv | sed -e '1,1d' -Ee 's/^.+__gene__([0-9]+\t)/\1/g' -e 's/;|; /__/g' | paste - ${GROUP1}_DRAM_results/tmp_pfam | awk -F "\t" '{if($4 ~ /[0-9]|function/) print $0}' >> ${GROUP1}_DRAM_annotations.txt

  ### add colummn for absent headers
  ### merge together once all headers present so columns are the same

done

fi

#################### STEP A7 #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A7" ]]; then

  conda activate anvio-7.1
  module load minimap2/2.17
  module load samtools

#### make BAM file of mapping reads vs contigs and sort and index it
### the same concatenate contig fasta must be used for all read files

for SETS in $(cat $input)
do
    READS=$(echo $SETS | cut -d ";" -f 2)
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    PROFILENAME=$(echo $SETS | cut -d ";" -f 4)
    minimap2 -t $SLURM_CPUS_PER_TASK -ax map-ont ${GROUP1}.fa $READS | samtools sort -O BAM - > ${PROFILENAME}.bam
    samtools index -@ $SLURM_CPUS_PER_TASK -b ${PROFILENAME}.bam
done

fi

#################### STEP A8 #################################
#### profile the contig databases WITH CONCATENATED FASTA OR CO_ASSEMBLIES ONLY - MULTIPLE SEQUENCING SAMPLES (SET OF READS) AND ONE CONTIG FASTA
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A8" ]]; then

  conda activate anvio-7.1

for SETS in $(cat $input)
do
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    PROFILENAME=$(echo $SETS | cut -d ";" -f 4)
    anvi-profile -i ${PROFILENAME}.bam -c ${GROUP1}.db -o profile_${PROFILENAME} -S "${PROFILENAME}" --min-coverage-for-variability 10 -T $SLURM_CPUS_PER_TASK -W
done

fi

#################### STEP A9 #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A9" ]]; then

  conda activate anvio-7.1

#### IMPORT Virsorter2 (before merging)

for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    anvi-import-functions -c ${GROUP1}.db -i ${GROUP1}_virsorter_annotations.txt
done

#### IMPORT DRAM results (before merging)

for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    anvi-import-functions -c ${GROUP1}.db -i ${GROUP1}_DRAM_annotations.txt
done

fi


#################### STEP B1 #################################
#### merge the database profiles

if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "B1" ]]; then

  conda activate anvio-7.1

#### prepare input lines, grouping profile.dbs in their groups as defined by SAMPLES1.txt
rm -f tmp_SAMPLES2_*
awk -F ";" 'BEGIN{OFS=";"};{print $4 > "tmp_SAMPLES2_"$3".1"}' $input
for f in tmp_SAMPLES2_*.1; do sed -E 's/(.+)/profile_\1\/PROFILE.db /g' $f | tr -d '\n' | sed 's/$/\n/g' > $(basename $f .1).2 ; done
for f in tmp_SAMPLES2_*.2; do awk -v f=$f '{print f";"$0}' $f | sed -E "s/tmp_SAMPLES2_(.+)\.2;/\1;/g" > $(basename $f .2).3; done
cat tmp_SAMPLES2_*.3 | sed 's/ /__/g' > SAMPLES2.txt
rm -f tmp_SAMPLES2_*

for SETS in $(cat SAMPLES2.txt)
do
GROUP2=$(echo $SETS | cut -d ";" -f 1)
PROFILELIST=$(echo $SETS | cut -d ";" -f 2 | sed 's/__/ /g')
anvi-merge --enforce-hierarchical-clustering ${PROFILELIST} -o profile_${GROUP2}-MERGED -c ${GROUP2}.db -W

#### post merge Virsorter2 import
anvi-import-collection ${GROUP2}_virsorter_collection.txt -c ${GROUP2}.db -p profile_${GROUP2}-MERGED/PROFILE.db -C VS2_${GROUP2}
anvi-import-misc-data ${GROUP2}_virsorter_additional_info.txt -p profile_${GROUP2}-MERGED/PROFILE.db --target-data-table items

#### prepare command for interactive interface and visualise
### NOT for job submission or batch; interactive only
## this generates a file that can be copy pasted from or use bash command on to run visualisation step in terminal
echo -e "conda activate anvio-7.1\n##############\nanvi-interactive -p $output/profile_${GROUP2}-MERGED/PROFILE.db -c $output/${GROUP2}.db --browser-path /usr/bin/chromium-browser" -V $output/${GROUP2}-RPKM.txt > ${GROUP2}_anvio_interactive_guide

echo -e "\n##############\nAfter MANUAL binning done, save bins as collection (named manual) via interactive interface, then do:\nconda activate anvio-7.1; anvi-summarize -p $output/profile_${GROUP2}-MERGED/PROFILE.db -c $output/${GROUP2}.db -o $output/${GROUP2}-SUMMARY -C manual\nIt would also be good to save an .svg of this figure with all bins visible\nAny images saved from interactive interface will be sent to ~/Downloads" >> ${GROUP2}_anvio_interactive_guide

echo -e "\n##############\nFor an interactive summary of results:\n/usr/bin/chromium-browser $output/${GROUP2}-SUMMARY/index.html" >> ${GROUP2}_anvio_interactive_guide

echo -e "\n##############\nAfter SUMMARY done, note which bins need refined (i.e. which bins are contaminated?), then for each bad bin, do:\nanvi-refine -p $output/profile_${GROUP2}-MERGED/PROFILE.db -c $output/${GROUP2}.db --browser-path /usr/bin/chromium-browser -C manual -b [BIN_NAME]" >> ${GROUP2}_anvio_interactive_guide


echo -e "\n###########\nTo add additional info to the layers, visible in the layer plots:\ncreate a tsv file e.g.:\nsamples density varroa_per_100_bees\nSU_TX_22_6083   0.0000000935    1.333333333\nSU_TX_22_6084   0.0000000935    6.333333333\nEQ_TX_22_5997   0.000000194     0\n\nthen do\nanvi-import-misc-data ${GROUP2}-additional_info.txt -p $output/profile_${GROUP2}-MERGED/PROFILE.db --target-data-table layers\n" >> ${GROUP2}_anvio_interactive_guide

echo -e "\n###########\nTo reorder the layers:\n create a tsv file with the data name, type, and sample names listed in the desired order e.g.:\nitem_name   data_type   data_value\nphylogeny   newick  (c2:0.0370199,(c1:0.0227268,c3:0.0227268)Int3:0.0370199)\ndensity   basic   c3,c2,c1\n\nthen do\nanvi-import-misc-data $output/${GROUP2}-layer_orders.txt -p $output/profile_${GROUP2}-MERGED/PROFILE.db --target-data-table layer_orders\n" >> ${GROUP2}_anvio_interactive_guide


done


fi


##################################################################
#################### STEP B2 #################################

#### prepare matrix of split coverage for normalisation using python script
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "B2" ]]; then

  for SETS in $(cat SAMPLES2.txt)
  do
  GROUP2=$(echo $SETS | cut -d ";" -f 1)
  PROFILELIST=$(echo $SETS | cut -d ";" -f 2 | sed 's/__/ /g')
  conda activate anvio-7.1

anvi-export-splits-and-coverages --splits-mode -p $output/profile_${GROUP2}-MERGED/PROFILE.db -c $output/${GROUP2}.db -o $output -O ${GROUP2}

$SEQKIT fx2tab -nl $output/${GROUP2}-SPLITS.fa | cut -f 2 | sed -z 's/^/length\n/g' | paste $output/${GROUP2}-COVs.txt - > $output/tmp_${GROUP2} && mv $output/tmp_${GROUP2} $output/${GROUP2}-COVs.txt

$COVERAGE2RPKM -i $output/${GROUP2}-COVs.txt -o $output/${GROUP2}-RPKM.txt

cut -f 1 $output/${GROUP2}-RPKM.txt | sed 's/_split.*//g' | sed 's/^contig$/__parent__/g' | paste $output/${GROUP2}-RPKM.txt - > $output/tmp_${GROUP2} && mv $output/tmp_${GROUP2} $output/${GROUP2}-RPKM.txt

done
fi

##################################################################



#################### STEP B3 #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "B3" ]]; then
  for SETS in $(cat SAMPLES2.txt)
  do
  GROUP2=$(echo $SETS | cut -d ";" -f 1)
  PROFILELIST=$(echo $SETS | cut -d ";" -f 2 | sed 's/__/ /g')

### Diamond blastx whole contigs against VOGDB
  $DIAMOND blastx -p $SLURM_CPUS_PER_TASK -d $DIAMONDDB -q ${GROUP2}-SPLITS.fa -f 6 -o ${GROUP2}.dmnd.blastx

conda activate DRAM

sed -E 's/\t(VOG[0-9]+)\./\t___\1___\t/1' ${GROUP2}.dmnd.blastx | awk -F "\t" '!a[$1,$2]++' | awk -F "\t" -v M=${GROUP2} '{print $0 > "tmp1_"M"_"$1}' ### keep only best hit per VOG per splot, then split results by split

## add extra VOG info
for f in tmp1_${GROUP2}_*; do
  awk -F";" '{print $1";"$5"/g"}' ${VOGDB}/vog_annotations.lca.reformatted_for_replacement.txt | sed -f - $f > tmp2_${f} && mv tmp2_${f} $f
done

### reduce info for each split to split name, number of vog hits, evalues comma separated, VOG codes comma separated
for f in tmp1_${GROUP2}_*; do
S=$(cut -f 1 $f | sort -Vu); W=$(wc -l $f | sed -E 's/^([0-9]+) .*/\1/g'); E=$(cut -f 12 $f | sed 's/$/,/g' | tr -d '\n' | sed 's/,$//g'); V=$(cut -f 2 $f | sed 's/$/,/g' | tr -d '\n' | sed 's/,$//g'); echo -e "$S\t$W\t$E\t"$V > tmp2_${f} && mv tmp2_${f} $f
done

cat tmp1_${GROUP2}_* | sed -z 's/^/split\tnum_VOG_hits\tVOG_Evalues\tVOG_info\n/1' > ${GROUP2}_dmndblastx_additional_info.txt

conda activate anvio-7.1
anvi-import-misc-data ${GROUP2}_dmndblastx_additional_info.txt -p profile_${GROUP2}-MERGED/PROFILE.db --target-data-table items

rm -f tmp1_${GROUP2}_*

done

fi





#################### STEP B4 #################################
#### OPTIONAL STEP To add additional info to the layers, visible in the layer plots
## requries a manually-made file

if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "B4" ]] && [[ ! -z "${layers}" ]]; then


  conda activate anvio-7.1

for SETS in $(cat SAMPLES2.txt)
do
GROUP2=$(echo $SETS | cut -d ";" -f 1)
PROFILELIST=$(echo $SETS | cut -d ";" -f 2 | sed 's/__/ /g')

anvi-import-misc-data ${layers} -p $output/profile_${GROUP2}-MERGED/PROFILE.db --target-data-table layers

done 
fi

#### OPTIONAL STEP To add additional order info for the layers
## requries a manually-made file

if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "B4" ]] && [[ ! -z "${order}" ]]; then

  conda activate anvio-7.1

for SETS in $(cat SAMPLES2.txt)
do
GROUP2=$(echo $SETS | cut -d ";" -f 1)
PROFILELIST=$(echo $SETS | cut -d ";" -f 2 | sed 's/__/ /g')

anvi-import-misc-data ${order} -p $output/profile_${GROUP2}-MERGED/PROFILE.db --target-data-table layer_orders

done
fi


############# enable permission for all in group ##################
chmod 777 -R $output