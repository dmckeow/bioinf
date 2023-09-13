#!/bin/bash

umask 002
green='\033[32m';red='\033[31m';cyan='\033[36m';nocolor='\033[m'
echo -e "${green}RUNNING bioinf-binning with the following parameters: $@ ${nocolor}"

Help()
{
echo -e "Bioinf-binning is a pipeline that uses anvi'o to perform binning on sequence datasets using custom databases setup by the user. It is intended to run on contigs from an assembly, and to maps the reads against this assembly. It can also be run without reads. Please note that currently the default mapping uses minimap and expects Oxford Nanopore reads. This script is designed to run on a single project per job submission, i.e. every sample in your input will be analysed together and end up in a single anvi'o figure/analyses\n"
echo -e "${green}########### REQUIRED parameters ###########${nocolor}\n-i --input\tAbsolute path to a plain text file with 3 columns separated by semicolons or tabs:\n"
echo -e "\tCOLUMN 1: /absolute_path/to/each/contig_fasta - separate line per assembly"
echo -e "\tCOLUMN 2: /absolute_path/to/each/reads_file - separate line per assembly, on same line with corresponding contigs"
echo -e "\tCOLUMN 3: profile name; a short name for each profile which will appear in all outputs and figures ${red} MUST NOT begin with a digit, MUST ONLY include letters, digits, and underscores, and MUST be unique${nocolor}"
echo -e "\t--input file example:\n\t\t${green}/home/data/19AUG22_1_DM_barcode01.contigs.fasta;/home/data/19AUG22_1_DM_barcode01-trimmedbcs.fastq.gz;location_1\n\t\t/home/data/19AUG22_1_DM_barcode03.contigs.fasta;/home/data/19AUG22_1_DM_barcode03-trimmedbcs.fastq.gz;location_2\n\t\t/home/data/19AUG22_1_DM_barcode06.contigs.fasta;/home/data/19AUG22_1_DM_barcode06-trimmedbcs.fastq.gz;location_3${nocolor}"
echo -e "\n-p --project\tA meaningful and unique name for your project. Output files will contain this name\n"
echo -e "${cyan}########### OPTIONAL parameters ###########${nocolor}"
echo -e "\n-m --mincontigsize\tinteger to set minimum size of contigs in bp to include in analyses (don't go below 1000)"
echo -e "-s --step\tWhich script steps to run. If not provided, the whole script is run. \nMultiple steps can be specified, e.g.: A1_A2_A3 will only run steps A1-A3.\n\t${green}STEPS AVAILABLE:${nocolor}"
awk '/^###### STEP-/ {print "\t\t"$0}' $(which bioinf-binning.sh)
echo -e "\n-n, --nomap [optional]\trun without mapping reads to contigs - allows you to run pipeline if you only have an assembly but no reads\n"
echo -e "-S --splitlength\tinteger to set minimum size of split in bp (default is 20000)\n"
echo -e "\n${green}example job submissions:${nocolor}\n"
echo -e "\tSUBMIT TO SLURM (see README for more info): ${cyan}sbatch --time=96:00:00 --cpus-per-task=24 --mem=240GB --partition long -o slurm.%N.%j.out -e slurm.%N.%j.err bioinf-binning.sh -i /path/to/samplelistfile -p project_name${nocolor}"
echo -e "\tSUBMIT TO SLURM, running only steps A1 to A3, and without mapping: ${cyan}sbatch --time=96:00:00 --cpus-per-task=24 --mem=240GB --partition long -o slurm.%N.%j.out -e slurm.%N.%j.err bioinf-binning.sh -i /path/to/samplelistfile -p project_name -s A1_A2_A3 -n${nocolor}"
echo -e "\tRUN LOCALLY (not recommended): ${cyan}bioinf-binning.sh -i /path/to/samplelistfile -p project_name --threads 24${nocolor}"
echo -e "\n-C, --concoct\trun auto binning with CONCOCT. Requires this script be run with mapping. Requires concoct setup through anvi'o as detailed in anvi'o installation\n"
echo -e "-t, --threads\tnumber of threads for job (SLURM --cpus-per-task does the same thing); default 1\n"
echo -e "-h, --help\tshow this help message and exit"
}

### set euclidean options to default setting here
nomap="false"
concoct="false"
##

while getopts i:p:s:m:S:nCt:h option
do 
    case "${option}" in
        i)input=${OPTARG};;
        p)project=${OPTARG};;
        ### optional:
        s)step=${OPTARG};;
        m)mincontigsize=${OPTARG};;
        S)splitlength=${OPTARG};;
        n)nomap="true";;
        C)concoct="true";;
        t)threads=${OPTARG};;
        h)Help; exit;;
    esac

done

####################### SET AND CHECK VARIABLES/ARGUMENTS #####################################
OUTDIR="${PWD}/BINNING_${project}"
mkdir -p ${bioinftmp} ## make tmp data storage in scratch

######### from commandline arugments ######
if [[ -z "${input}" ]]; then echo -e "${red}-i, --input missing"; Help; exit; fi
if [[ -z "${project}" ]]; then echo -e "${red}-p, --project missing"; Help; exit; fi

### set defaults for optional arguments
if [[ -z "${splitlength}" ]]; then splitlength="20000"; echo -e "${green}Bioinf-binning is using default split length 20000 bp ${nocolor}"; fi
if [[ -z "${mincontigsize}" ]]; then mincontigsize="2000"; echo -e "${green}Bioinf-binning is using default minimum contig length 2000 bp ${nocolor}"; fi
if [[ "${nomap}" == "false" ]]; then echo -e "${green}Bioinf-binning running WITH mapping${nocolor}"; fi
if [[ "${nomap}" == "true" ]]; then echo -e "${green}Bioinf-binning in nomap mode - no read mapping will be done${nocolor}"; fi

### THREADS

if [[ ! -z "${SLURM_CPUS_PER_TASK}" ]]; then THREADS="${SLURM_CPUS_PER_TASK}"; echo -e "${green}Using threads set by SLURM_CPUS_PER_TASK:${nocolor}"; fi

if [[ ! -z "${threads}" ]]; then THREADS="${threads}"; echo -e "${green}Using threads set by --threads flag:${nocolor}"; fi

if [[ -z "${threads}" ]] && [[ -z "${SLURM_CPUS_PER_TASK}" ]]; then THREADS="1"; echo -e "${green}No SLURM_CPUS_PER_TASK or --threads set, using default threads:${nocolor}"; fi
echo -e "\t${THREADS} threads"

#######
if [[ -z $(conda env list | grep "bioinftools") ]]; then echo "NO conda environment for bioinftools found - see README"; else echo "conda environment for bioinftools FOUND"; fi

if [ -z ${bioinfdb+x} ]; then echo -e "${red}bioinf database is unset - see README and bioinf-setup.sh in bioinf directory ${nocolor}"; exit; fi
echo -e "${green}bioinf database is set to '$bioinfdb' ${nocolor}"


if [ -z ${bioinftmp+x} ]; then echo -e "${red}bioinf temporary working space is unset - see README and bioinf-setup.sh in bioinf directory ${nocolor}"; exit; fi
echo -e "${green}bioinf temporary working space is set to '$bioinftmp' ${nocolor}"

if [[ ! -d "$bioinfdb" ]]; then echo -e "${red}NO database directory for bioinftools found - see README and bioinf-setup.sh in bioinf directory ${nocolor}"; exit; fi

if [[ -z "$(ls -A "$bioinfdb"/DMND)" ]]; then echo -e "${red}NO databases for DIAMOND found - see README and bioinf-setup.sh in bioinf directory ${nocolor}"; fi

if [[ -z "$(ls -A "$bioinfdb"/HMM)" ]]; then echo -e "${red}NO databases for HMMS found - see README and bioinf-setup.sh in bioinf directory ${nocolor}"; fi

if [[ -z "$(ls -A "$bioinfdb"/KAIJU/*)" ]]; then echo -e "${red}NO databases for Kaiju found - see README and bioinf-setup.sh in bioinf directory ${nocolor}"; fi

if [[ -z "$(ls -A "$bioinfdb"/VOGDB)" ]]; then echo -e "${red}NO databases for VOGDB found - see README and bioinf-setup.sh in bioinf directory ${nocolor}"; fi


####################### SOFTWARE #####################################
### LOAD software available via shared environment on server:
module purge
eval "$(conda shell.bash hook)"
conda activate bioinftools

####################### DATABASES #####################################

########################################################################
########################################################################

####################### THE SCRIPT #####################################

################################################################################################
################################################################################################
########################################################################
########################################################################

###### STEP-A1 - pre-process your contig fasta(s)
if [[ -z "${step}" ]] || [[ "$step" =~ "A1" ]]; then
echo -e "${cyan}\t\tRUNNING STEP A1${nocolor}"
### make output parent directories if they do not already exist

dos2unix $input
sed -i 's/\t/;/g' $input

##### make the output directory
rm -fr $OUTDIR
mkdir $OUTDIR
cd $OUTDIR

rm -fr ${project}-reformat
mkdir ${project}-reformat
cd ${project}-reformat

### fix fasta deflines and save "deflinekey" files to show corresponding original format and new format deflines
for SETS in $(cat $input)
do
  CONTIGS=$(echo $SETS | cut -d ";" -f 1)
  PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
  anvi-script-reformat-fasta $CONTIGS -o ${PROFILENAME}.fa --simplify-names --report-file ${PROFILENAME}.deflinekey -l $mincontigsize --seq-type NT --prefix "contig_${PROFILENAME}"
done

cd $OUTDIR
rm -f ${project}.fa
rm -f ${project}.deflinekey
touch ${project}.fa
touch ${project}.deflinekey

cat ${project}-reformat/*.fa >> ${project}.fa
cat ${project}-reformat/*.deflinekey >> ${project}.deflinekey
rm -fr ${project}-reformat

##############################################################
fi
##############################################################


###### STEP-A2 - generate the contigs database
if [[ -z "${step}" ]] || [[ "$step" =~ "A2" ]]; then
echo -e "${cyan}\t\tRUNNING STEP A2${nocolor}"

cd $OUTDIR

### generate the contigs database
##### uses user provided split length from command line if provided

    if [ -z "$splitlength" ]; then
        anvi-gen-contigs-database -f ${project}.fa -o ${project}.db -n "${project}" -T $THREADS
    fi

    if [ ! -z "$splitlength" ]; then
        anvi-gen-contigs-database -f ${project}.fa -o ${project}.db -n "${project}" -T $THREADS --split-length ${splitlength}
    fi

##############################################################
fi
##############################################################

###### STEP-A3 - perform hmm searches vs default databases [optional]
if [[ -z "${step}" ]] || [[ "$step" =~ "A3" ]]; then
echo -e "${cyan}\t\tRUNNING STEP A3${nocolor}"

cd $OUTDIR

#### run hmm profile, summarise contig stats, and identify COGs
### !!! both anvi-run-hmms and anvi-run-ncbi-cogs do NOT include viruses and most eukaryotes
### !!! as an alternative to the default hmm profiles of anvio, we will later run and import the results of Virsorter2

### default anvio HMM profile dbs
    anvi-run-hmms -c ${project}.db -T $THREADS --just-do-it
    anvi-run-ncbi-cogs -c ${project}.db -T $THREADS
    anvi-run-kegg-kofams -c ${project}.db -T $THREADS --just-do-it --skip-bitscore-heuristic

##############################################################
fi
##############################################################

###### STEP-A4 - perform hmm searches vs custom databases, estimate taxonomy [optional]
if [[ -z "${step}" ]] || [[ "$step" =~ "A4" ]]; then
echo -e "${cyan}\t\tRUNNING STEP A4${nocolor}"

cd $OUTDIR

#### run hmm profile, summarise contig stats, and identify COGs
### !!! both anvi-run-hmms and anvi-run-ncbi-cogs do NOT include viruses and most eukaryotes

### custom anvio HMM profile dbs 

for f in "$bioinfdb"/HMM/*; do
    anvi-run-hmms -c ${project}.db -T $THREADS -H $f --hmmer-program hmmsearch --just-do-it
done


### run single core gene taxonomy - this is built into anvio and cant be customised currently. It provides quick taxonomy for prokaryotes.

  anvi-run-scg-taxonomy -c ${project}.db -T $THREADS
  anvi-estimate-scg-taxonomy -c ${project}.db --metagenome-mode -T $THREADS


###### export HMMs hits per gene and import them as annotations
for f in "$bioinfdb"/HMM/*; do
    rm -fr ${project}_$(basename $f)-annot; mkdir ${project}_$(basename $f)-annot
    anvi-script-get-hmm-hits-per-gene-call -c ${project}.db --hmm-source $(basename $f) -o ${project}_$(basename $f)-annot/hits_per_gene_call.tmp1
done

### get rid of the empty hmm db hit results so anvio doesnt lose its shit
find ${project}_*-annot/* -type f -empty -delete
find ${project}_*-annot -type d -empty -delete

### final reformat for the HMM profiles with hits
for f in ${project}_*-annot; do
    awk -F "\t" '!a[$2]++' ${f}/hits_per_gene_call.tmp1 | cut --complement -f 1 | sed 's/$/\t1/g' | sed 's/gene_name\tgene_hmm_id\t1/function\taccession\te_value/1' > ${f}/hits_per_gene_call.txt
done


###### finally,  import the custom hmm hits as gene annotations
for f in ${project}_*-annot; do
    anvi-import-functions -c ${project}.db -i ${f}/hits_per_gene_call.txt
done

##############################################################
fi
##############################################################

###### STEP-A5 - get taxonomy for genes with KAIJU and import into contigs database [optional]
if [[ -z "${step}" ]] || [[ "$step" =~ "A5" ]]; then
echo -e "${cyan}\t\tRUNNING STEP A5${nocolor}"

cd $OUTDIR

### get gene calls 
anvi-get-sequences-for-gene-calls -c ${project}.db -o ${project}_gene_calls.fa

### run kaiju for all kaiju databases in the bioinfdb
for f in "$bioinfdb"/KAIJU/*; do
    kaiju -t ${f}/nodes.dmp -f ${f}/*.fmi -i ${project}_gene_calls.fa -o ${project}_gene_calls.kaiju.$(basename $f) -z $THREADS -v
    sort -t $'\t' -V -k 2,2 ${project}_gene_calls.kaiju.$(basename $f) -o ${project}_gene_calls.kaiju.$(basename $f)
    kaiju-addTaxonNames -t ${f}/nodes.dmp -n ${f}/names.dmp -i ${project}_gene_calls.kaiju.$(basename $f) -o ${project}_gene_calls.kaiju.$(basename $f).names -r superkingdom,phylum,order,class,family,genus,species
done

#### cat all the kaiju hits together, sort by gene call name and then kaiju score in reverse order and keep the first line for each gene, i.e. the hit with the highest score. Kaiju seems to have similar scores even with different database sizes
awk '{print $0"\t0\tNA\tNA\tNA\tNA"}' ${project}_gene_calls.kaiju.*.names | cut -f 1-8 | sort -t $'\t' -k 2,2V -k 4,4nr | awk -F "\t" '!a[$2]++' | sed 's/0\tNA\tNA\tNA\tNA//g' > ${project}.kaiju.merge

anvi-import-taxonomy-for-genes -c ${project}.db -i ${project}.kaiju.merge -p kaiju --just-do-it

##############################################################
fi
##############################################################


###### STEP-A6 - mapping reads vs contigs to create bam files [nomap=skip]
if [[ -z "${step}" ]] || [[ "$step" =~ "A6" ]] && [[ "$nomap" == "false" ]]; then
    echo -e "${cyan}\t\tRUNNING STEP A6${nocolor}"

cd $OUTDIR
rm -fr MAPPING-${project}
mkdir MAPPING-${project}

#### make BAM file of mapping reads vs contigs and sort and index it
### the same concatenate contig fasta must be used for all read files

for SETS in $(cat $input)
do
    READS=$(echo $SETS | cut -d ";" -f 2)
    PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
    minimap2 -t $THREADS -ax map-ont ${project}.fa $READS | samtools sort -O BAM - > MAPPING-${project}/${PROFILENAME}.bam
    samtools index -@ $THREADS -b MAPPING-${project}/${PROFILENAME}.bam
done

##############################################################
fi
##############################################################

###### STEP-A7 - create profile for each sample [nomap=pseudoskip]

##############################################################
if [[ -z "${step}" ]] || [[ "$step" =~ "A7" ]]; then
##############################################################
echo -e "${cyan}\t\tRUNNING STEP A7${nocolor}"

  cd $OUTDIR

####### with mapping
for SETS in $(cat $input)
do
    PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
    if [[ "$nomap" == "false" ]]; then
        anvi-profile -i MAPPING-${project}/${PROFILENAME}.bam -c ${project}.db -o profile_${PROFILENAME} -S "${PROFILENAME}" --min-coverage-for-variability 10 -T $THREADS -M $mincontigsize -W
    fi

####### with NO mapping
    if [[ "$nomap" == "true" ]]; then
        anvi-profile -i MAPPING-${project}/${PROFILENAME}.bam -c ${project}.db -o profile_${PROFILENAME} -S "${PROFILENAME}" --min-coverage-for-variability 10 -T $THREADS -M $mincontigsize -W --blank-profile
    fi

done

##############################################################
fi
##############################################################


###### STEP-A8 - Diamond blastx whole contigs [optional]
if [[ -z "${step}" ]] || [[ "$step" =~ "A8" ]]; then
    echo -e "${cyan}\t\tRUNNING STEP A8${nocolor}"

cd $OUTDIR

anvi-export-contigs --splits-mode -c ${project}.db -o ${project}-SPLITS.fa

### Diamond blastx whole splits against custom dmnd databases
### reduce to best hit per split

####MAXTARGETS=$(seqkit stat -T ${project}-SPLITS.fa | sed '1,1d' | awk -F "\t" '{printf "%0.0f\n", $8/100}')

for f in ${bioinfdb}/DMND/*; do
    ##diamond blastx --max-target-seqs $MAXTARGETS --evalue 1e-20 --sensitive -p $THREADS -d $f -q ${project}-SPLITS.fa -f 6 -o $(basename $f .dmnd).evalue_1e-20.dmnd.blastx
    diamond blastx --range-culling --top 20 -F 15 --evalue 1e-20 --sensitive -p $THREADS -d $f -q ${project}-SPLITS.fa -f 6 -o $(basename $f .dmnd).evalue_1e-20.dmnd.blastx
done

### best hit per split subject combo
for f in *.evalue_*.dmnd.blastx; do
    awk -F "\t" '!a[$1,$2]++' $f > ${f}.tmp1
done


### reduce info for each split to split name, number of hits

for f in *.evalue_*.dmnd.blastx.tmp1; do
    cut -f 1 $f | sort -V | uniq -c | sed -E 's/ +([0-9]+) (.*)/\2\t\1/g' > $(basename $f .tmp1).tmp2
done

for f in *.evalue_*.dmnd.blastx.tmp1; do
    awk -F "\t" '!a[$1]++' $f | cut -f 1,2 | sort -V -t $'\t' -k 1,1 | cut -f 2 | paste $(basename $f .tmp1).tmp2 - > $(basename $f .tmp1).txt
        sed -i -z "s/^/split\t${f}_num_hits\t${f}_tophit\n/1" $(basename $f .tmp1).txt
        sed -i 's/\.tmp1//g' $(basename $f .tmp1).txt
done

for SETS in $(cat $input)
do
    PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
        for f in *.evalue_*.dmnd.blastx.txt; do
            anvi-import-misc-data -p profile_${PROFILENAME}/PROFILE.db --target-data-table items --just-do-it ${f}
        done
done


rm -f *.evalue_*.dmnd.blastx.tmp*

####################### BLAST #############################
### IF any blast database files are found in bioinfdb/BLAST, then also do a BLASTx search against those

for f in $(ls ${bioinfdb}/BLAST | cut -d "." -f 1 | sort -Vu); do
    blastx -evalue 1e-20 -num_threads $THREADS -db ${bioinfdb}/BLAST/$f -query ${project}-SPLITS.fa -outfmt 6 -out $(basename $f).evalue_1e-20.ncbi.blastx
done

### best hit per split subject combo
for f in *.evalue_*.ncbi.blastx; do
    awk -F "\t" '!a[$1,$2]++' $f > ${f}.tmp1
done


### reduce info for each split to split name, number of hits

for f in *.evalue_*.ncbi.blastx.tmp1; do
    cut -f 1 $f | sort -V | uniq -c | sed -E 's/ +([0-9]+) (.*)/\2\t\1/g' > $(basename $f .tmp1).tmp2
done

for f in *.evalue_*.ncbi.blastx.tmp1; do
    awk -F "\t" '!a[$1]++' $f | cut -f 1,2 | sort -V -t $'\t' -k 1,1 | cut -f 2 | paste $(basename $f .tmp1).tmp2 - > $(basename $f .tmp1).txt
        sed -i -z "s/^/split\t${f}_num_hits\t${f}_tophit\n/1" $(basename $f .tmp1).txt
        sed -i 's/\.tmp1//g' $(basename $f .tmp1).txt
done

for SETS in $(cat $input)
do
    PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
        for f in *.evalue_*.ncbi.blastx.txt; do
            anvi-import-misc-data -p profile_${PROFILENAME}/PROFILE.db --target-data-table items --just-do-it ${f}
        done
done


rm -f *.evalue_*.ncbi.blastx.tmp*

##############################################################
fi
##############################################################



##############################################################
###### STEP-A9 - merge the profiles together [nomap=pseudoskip]
##############################################################

if [[ -z "${step}" ]] || [[ "$step" =~ "A9" ]]; then
    echo -e "${cyan}\t\tRUNNING STEP A9${nocolor}"

cd $OUTDIR

### import a collection that is just all splits in a single meaningless bin (this allows us to summarize a nomap project that is too big to opened - we can then filter the anvi-summarise results and select a contig subset based on HMM hits etc to create smaller anvio projects that maybe we can actually open)
anvi-export-contigs --splits-mode -c ${project}.db -o ${project}-SPLITS.fa
grep ">" ${project}-SPLITS.fa | sed -E 's/>(.*)/\1\tBin_1/g' > ${project}-SPLITS-COLLECTION-ALL

    if [[ "$nomap" == "false" ]]; then

        PROFILELIST=$(cut -d ";" -f 3 ${input} | sed -E 's/(.+)/profile_\1\/PROFILE.db /g' | tr -d '\n' | sed 's/$/\n/g')
        rm -fr profile_${project}-MERGED

        anvi-merge ${PROFILELIST} -o profile_${project}-MERGED -c ${project}.db -W
        a="${OUTDIR}/profile_${project}-MERGED/PROFILE.db"

            mergedfile="profile_${project}-MERGED/PROFILE.db"
            if [[ -f "$mergedfile" ]]; then
                for f in *.evalue_*.blastx.txt; do
                    anvi-import-misc-data -p profile_${project}-MERGED/PROFILE.db --target-data-table items --just-do-it ${f}
                done

                anvi-import-collection -C SPLITS_COLLECTION_ALL -p profile_${project}-MERGED/PROFILE.db -c ${project}.db ${project}-SPLITS-COLLECTION-ALL
            fi

    fi


    if [[ "$nomap" == "true" ]]; then
        for SETS in $(cat $input)
            do
                PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
                a="${OUTDIR}/profile_${PROFILENAME}/PROFILE.db"
                anvi-import-collection -C SPLITS_COLLECTION_ALL -p ${OUTDIR}/profile_${PROFILENAME}/PROFILE.db -c ${project}.db ${project}-SPLITS-COLLECTION-ALL

            done
    fi


#### prepare command for interactive interface and visualise
### NOT for job submission or batch; interactive only
## this generates a file that can be copy pasted from or use bash command on to run visualisation step in terminal
echo -e "\n##### To run your interactive interface !!! \n#######\n${cyan}conda activate bioinftools\n##############${nocolor}\nIF logged in with:\nssh -L 8080:localhost:8080\ndo\n${cyan}anvi-interactive -p ${a} -c ${OUTDIR}/${project}.db --server-only --ip-address 127.0.0.1 --port-number 8080\n\t and then copy paste the server address that anvio shows into your browser (use chrome), e.g.: http://127.0.0.1:8080${nocolor}" > BINNING_GUIDE_${project}

echo -e "\n##### OR to run the interactive interface through your local system (performance can be MUCH worse this way - also specify the path to your browser if needed): ##############\n\n${cyan}anvi-interactive -p ${a} -c ${OUTDIR}/${project}.db --browser-path /usr/bin/chromium-browser${nocolor}\n" >> BINNING_GUIDE_${project}

echo -e "\n##############\nAfter MANUAL binning done, save bins as collection (named manual in this example) via interactive interface, then do:\n\n${green}anvi-summarize -p ${a} -c ${OUTDIR}/${project}.db -o ${OUTDIR}/${project}-SUMMARY -C manual${nocolor}\n" >> BINNING_GUIDE_${project}

echo -e "\n##############\nAfter anvi-summarize is run, you can optionally refine individual bins. Note which bins need refined (i.e. which bins are contaminated?), then for each bad bin, do:\n${green}anvi-refine -p ${a} -c ${OUTDIR}/${project}.db --browser-path /usr/bin/chromium-browser -C manual -b [BIN_NAME]${nocolor}" >> BINNING_GUIDE_${project}

echo -e "\n###########\nTo add additional info to the layers, visible in the layer plots: !!! \ncreate a tsv file e.g.:\nsamples density varroa_per_100_bees\nSU_TX_22_6083   0.0000000935    1.333333333\nSU_TX_22_6084   0.0000000935    6.333333333\nEQ_TX_22_5997   0.000000194     0\n\nthen do\nanvi-import-misc-data ${project}-additional_info.txt -p ${a} --target-data-table layers\n" >> BINNING_GUIDE_${project}

echo -e "\n###########\nTo reorder the layers: !!!\n create a tsv file with the data name, type, and sample names listed in the desired order e.g.:\nitem_name   data_type   data_value\nphylogeny   newick  (c2:0.0370199,(c1:0.0227268,c3:0.0227268)Int3:0.0370199)\ndensity   basic   c3,c2,c1\n\nthen do\nanvi-import-misc-data ${OUTDIR}/${project}-layer_orders.txt -p ${a} --target-data-table layer_orders\n" >> BINNING_GUIDE_${project}

##############################################################
fi
##############################################################

##############################################################
###### STEP-B1 - cluster the contigs using concoct - essentially automatic binning [nomap=skip]
##############################################################
if ([[ -z "${step}" ]] || [[ "$step" =~ "B1" ]]) && [[ "$concoct" == "true" ]]; then
        echo -e "${cyan}\t\tRUNNING STEP B1${nocolor}"
cd $OUTDIR

if [[ "$nomap" == "false" ]]; then
    mergedfile="profile_${project}-MERGED/PROFILE.db"
        if [[ -f "$mergedfile" ]]; then
            anvi-cluster-contigs -p profile_${project}-MERGED/PROFILE.db -c ${project}.db -C concoct --driver concoct --just-do-it
        fi
fi

##############################################################
fi
##############################################################

##################################################################
#################### X - prepare normalised counts of coverage (RPKM) [optional, nomap-skip]

#### prepare matrix of split coverage for normalisation using python script
#if [[ "$step" =~ "removed" ]] || [[ "$step" =~ "X" ]] && [[ ! "$step" =~ "nomap" ]]; then
 # cd ${PWD}/ANVIO_${project}

  #module purge; conda deactivate
  #conda activate ${bioinftools}

#anvi-export-splits-and-coverages --splits-mode -p profile_${project}-MERGED/PROFILE.db -c ${project}.db -o . -O ${project}

#$SEQKIT fx2tab -nl ${project}-SPLITS.fa | cut -f 2 | sed -z 's/^/length\n/g' | paste ${project}-COVs.txt - > tmp_${project} && mv tmp_${project} ${project}-COVs.txt

#$COVERAGE2RPKM -i ${project}-COVs.txt -o ${project}-RPKM.txt

#cut -f 1 ${project}-RPKM.txt | sed 's/_split.*//g' | sed 's/^contig$/__parent__/g' | paste ${project}-RPKM.txt - > tmp_${project} && mv tmp_${project} ${project}-RPKM.txt


##############################################################
###fi
##############################################################
##################################################################

