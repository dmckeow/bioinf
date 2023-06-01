#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=240GB
#SBATCH --partition ag2tb
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

echo -e "COMMAND LINE JOB SUBMISSSION:\n\tsbatch /home/dcschroe/dmckeow/projects/DWV/script/anvio2.sh $@"

Help()
{
echo -e "########### This script is designed to run on a single project per job submission, i.e. every sample in your input will be analysed together and end up in a single anvi'o figure/analyses ###########\n"
echo -e "########### REQUIRED arguments ###########\n-i --input\tfull path to a file that lists input file locations, semi-colon delimited, 3 fields:\n"
echo -e "\tFIELD1=contig fasta, separate line per assembly (full path). Can be individual assemblies or co-assemblies."
echo -e "\tFIELD2=reads fq.gz, corresponding to each assembly contig fasta (full path)"
echo -e "\tFIELD3=profile name; a short name for each profile which will appear in all outputs and figures."
echo -e "!!! profile names MUST NOT begin with a digit, and MUST ONLY include letters, digits, and underscores and MUST be UNIQUE within your samples being run!!!\n"
echo -e "! example file for --input:\n\n/home/data/19AUG22_1_DM_barcode01.contigs.fasta;/home/data/19AUG22_1_DM_barcode01-trimmedbcs.fastq.gz;location1\n/home/data/19AUG22_1_DM_barcode03.contigs.fasta;/home/data/19AUG22_1_DM_barcode03-trimmedbcs.fastq.gz;location2\n/home/data/19AUG22_1_DM_barcode06.contigs.fasta;/home/data/19AUG22_1_DM_barcode06-trimmedbcs.fastq.gz;location3\n\n"
echo -e "-s --step\tWhich script steps to run. A0 = run the whole script. \nMultiple steps can be specified, e.g.: A1_A2_A3 will only run steps A1-A3.\n\nIf running WITHOUT mapping - i.e. you only have contigs, but no reads, then add nomap to -s , e.g. A0_nomap\n\nALSO, some steps can be resumed without deleting previous progress for that step (see below), by adding the flag resume to -s but if you resume I recommend that you delete the last output file of that step, in case it is an incomplete output file\n\n IF using nomap then submit one job per sample as they cannot be merged and set -prefix to be the same as the profile name\n\n!!! STEPS AVAILABLE: !!!"
awk '/#+ STEP /' /home/dcschroe/dmckeow/projects/DWV/script/anvio2.sh | sed -E -e 's/#+ *//g' -e 's/^STEP /> /g'
echo -e "\n-m --mincontigsize\tinteger to set minimum size of contigs in bp to include in analyses (don't go below 1000)"
echo -e "\n-o --output\tabsolute path for final output directory. Will be created if it does not already exist. An output folder named after ANVIO_(prefix) will be created here\n"
echo -e "\n-p --prefix\ta meaningful name for output folder. Any pre-existing driectory at this location and identical name will be deleted. Will also be the name of the merged analyses files. It should be a UNIQUE name that cannot possibly be identical to any other folder within your output destination e.g. -p 11OCT22project will make the folder ANVIO_11OCT22project\n"
echo -e "########### OPTIONAL arguments ###########\n\n"
echo -e "\n-S --splitlength\tinteger to set minimum size of split in bp (default is 20000)\n"
echo -e "\n! example job submission (to run the whole script with mapping). In these examples, the script will make a folder called ANVIO_JAN23_PHEX within the /home/data folder:\nsbatch /home/dcschroe/dmckeow/projects/DWV/script/anvio.sh -i /home/data/anvio_JAN23_PHEX_SAMPLES1.txt -s A0 -m 2500 -o /home/data -p JAN23_PHEX\n"
echo -e "\n! example job submission (to run only steps 1-3 and 8 without mapping, and to resume any steps that can be resumed):\nsbatch /home/dcschroe/dmckeow/projects/DWV/script/anvio.sh -i /home/data/anvio_JAN23_PHEX_SAMPLES1.txt -s A1_A2_A3_A8_nomap_resume -m 2500 -o /home/data -p JAN23_PHEX\n"
echo -e "-h, --help\tshow this help message and exit\n"
}


while getopts i:s:m:o:p:S:h option
do 
    case "${option}" in
        i)input=${OPTARG};;
        s)step=${OPTARG};;
    m)mincontigsize=${OPTARG};;
    o)output=${OPTARG};;
    p)prefix=${OPTARG};;
    S)splitlength=${OPTARG};;
    h)Help; exit;;
    esac

done

if [[ -z "${input}" ]]; then echo "-i, --input REQUIRED"; Help; exit; fi
if [[ -z "${step}" ]]; then echo "-s, --step REQUIRED"; Help; exit; fi
if [[ -z "${mincontigsize}" ]]; then echo "-m, --mincontigsize REQUIRED"; Help; exit; fi
if [[ -z "${output}" ]]; then echo "-o, --output REQUIRED"; Help; exit; fi
if [[ -z "${prefix}" ]]; then echo "-p, --prefix REQUIRED"; Help; exit; fi


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

#### a directory that contains ALL of the CUSTOM HMM profiles you want to run
CUSTOM_HMMS="/panfs/jay/groups/27/dcschroe/dmckeow/custom_HMM_anvio"
CUSTOM_DMND="/panfs/jay/groups/27/dcschroe/dmckeow/custom_DMND"

### KAIJU databases
kaiju_rvdb="/panfs/jay/groups/27/dcschroe/shared/tools/kaijudb/rvdb"
kaiju_nr_euk="/panfs/jay/groups/27/dcschroe/shared/tools/kaijudb/nr_euk" ## NOTE this is a modified version of nr_euk from Kaiju that actually includes all eukaryotes

########################################################################
########################################################################

######### SOFTWARE/OTHER SCRIPTS THAT NEED SETUP BEFORE RUNNING:

COVERAGE2RPKM="/home/dcschroe/dmckeow/projects/DWV/script/coverage_to_RPKM.py" ## custom script that converts coverage data into RPKM

##### Seqkit - must be installed locally, e.g.:
SEQKIT="/panfs/jay/groups/27/dcschroe/shared/tools/seqkit"

##### taxonkit - must be installed locally, e.g.:
TAXONKIT="/home/dcschroe/dmckeow/bin/taxonkit"
TAXONKITDB="/panfs/jay/groups/27/dcschroe/dmckeow/.taxonkit"

##### ANVIO 7.1 or later - INSTALLED via CONDA
ANVIO="/panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/anvio-7.1"

##### Kaiju - INSTALLED via CONDA
KAIJU="/panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/kaiju"

##### Diamond - installed locally
DIAMOND="/panfs/jay/groups/27/dcschroe/shared/tools/diamond"


######################################################################
#################VOGDB#####################################################


####################### THE SCRIPT #####################################

#################### STEP A0 - runs STEPS A1-A9 & B1 #################################

################################################################################################
################################################################################################
########################################################################
########################################################################

#################### STEP A1 - pre-process your contig fasta(s) #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A1" ]]; then
  module purge; conda deactivate
  conda activate $ANVIO

### make output parent directories if they do not already exist
mkdir -p ${output}

dos2unix $input
sed -i 's/\t/;/g' $input

##### make the final output directory
rm -fr ${output}/ANVIO_${prefix}
mkdir ${output}/ANVIO_${prefix}
cd ${output}/ANVIO_${prefix}

rm -fr ${prefix}-reformat
mkdir ${prefix}-reformat
cd ${prefix}-reformat

### fix fasta deflines and save "deflinekey" files to show corresponding original format and new format deflines
for SETS in $(cat $input)
do
  CONTIGS=$(echo $SETS | cut -d ";" -f 1)
  PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
  anvi-script-reformat-fasta $CONTIGS -o ${PROFILENAME}.fa --simplify-names --report-file ${PROFILENAME}.deflinekey -l $mincontigsize --seq-type NT --prefix "contig_${PROFILENAME}"
done

cd ${output}/ANVIO_${prefix}
rm -f ${prefix}.fa
rm -f ${prefix}.deflinekey
touch ${prefix}.fa
touch ${prefix}.deflinekey

cat ${prefix}-reformat/*.fa >> ${prefix}.fa
cat ${prefix}-reformat/*.deflinekey >> ${prefix}.deflinekey
rm -fr ${prefix}-reformat

##############################################################
fi
##############################################################

#################### STEP A2 - generate the contigs database [nomap resume] #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A2" ]]; then
  module purge; conda deactivate
  conda activate $ANVIO

cd ${output}/ANVIO_${prefix}

### generate the contigs database

####### delete existing db if no resume flag
    if [[ ! "$step" =~ "resume" ]]; then rm -f ${prefix}.db; fi

##### uses user provided split length from command line if provided

    if [ -z "$splitlength" ]; then
        anvi-gen-contigs-database -f ${prefix}.fa -o ${prefix}.db -n "${prefix}" -T $SLURM_CPUS_PER_TASK
    fi

    if [ ! -z "$splitlength" ]; then
        anvi-gen-contigs-database -f ${prefix}.fa -o ${prefix}.db -n "${prefix}" -T $SLURM_CPUS_PER_TASK --split-length ${splitlength}
    fi

##############################################################
fi
##############################################################

#################### STEP A3 - perform hmm searches vs default databases (optional) #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A3" ]]; then
  module purge; conda deactivate
  conda activate $ANVIO
cd ${output}/ANVIO_${prefix}

#### run hmm profile, summarise contig stats, and identify COGs
### !!! both anvi-run-hmms and anvi-run-ncbi-cogs do NOT include viruses and most eukaryotes
### !!! as an alternative to the default hmm profiles of anvio, we will later run and import the results of Virsorter2

### default anvio HMM profile dbs
    anvi-run-hmms -c ${prefix}.db -T $SLURM_CPUS_PER_TASK --just-do-it
    anvi-run-ncbi-cogs -c ${prefix}.db -T $SLURM_CPUS_PER_TASK
    anvi-run-kegg-kofams -c ${prefix}.db -T $SLURM_CPUS_PER_TASK --just-do-it --skip-bitscore-heuristic

##############################################################
fi
##############################################################

#################### STEP A4 - perform hmm searches vs custom databases, estimate taxonomy (optional)) #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A4" ]]; then
  module purge; conda deactivate
  conda activate $ANVIO
cd ${output}/ANVIO_${prefix}

#### run hmm profile, summarise contig stats, and identify COGs
### !!! both anvi-run-hmms and anvi-run-ncbi-cogs do NOT include viruses and most eukaryotes

### custom anvio HMM profile dbs - best with HIGH (E-20), MED(E-60), LOW(E-180) evalue copies of each db

for f in $(realpath $CUSTOM_HMMS/*); do
    anvi-run-hmms -c ${prefix}.db -T $SLURM_CPUS_PER_TASK -H $f --hmmer-program hmmsearch --just-do-it
done


### run single core gene taxonomy - this is built into anvio and cant be customised currently. It provides quick taxonomy for prokaryotes.
### !!! the alternative is to use Kraken2 and then import the results, which will provide gene level taxonomy

  anvi-run-scg-taxonomy -c ${prefix}.db -T $SLURM_CPUS_PER_TASK
  anvi-estimate-scg-taxonomy -c ${prefix}.db --metagenome-mode -T $SLURM_CPUS_PER_TASK


###### export HMMs hits per gene and import them as annotations
for f in $(realpath $CUSTOM_HMMS/*); do
    rm -fr ${prefix}_$(basename $f)-annot; mkdir ${prefix}_$(basename $f)-annot
    anvi-script-get-hmm-hits-per-gene-call -c ${prefix}.db --hmm-source $(basename $f) -o ${prefix}_$(basename $f)-annot/hits_per_gene_call.tmp1
    awk -F "\t" '!a[$2]++' ${prefix}_$(basename $f)-annot/hits_per_gene_call.tmp1 | cut --complement -f 1 | sed 's/$/\t1/g' | sed 's/gene_name\tgene_hmm_id\t1/function\taccession\te_value/1' > ${prefix}_$(basename $f)-annot/hits_per_gene_call.txt
done

### get rid of the empty hmm db hit results so anvio doesnt lose its shit
find ${prefix}_*-annot -type f -empty -delete
find ${prefix}_*-annot -type d -empty -delete

###### finally,  import the custom hmm hits as gene annotations
for f in $(realpath $CUSTOM_HMMS/*); do
    anvi-import-functions -c ${prefix}.db -i ${prefix}_$(basename $f)-annot/hits_per_gene_call.txt
done

##############################################################
fi
##############################################################



#################### STEP A5 - get taxonomy for genes with KAIJU and import into contigs database (optional) #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A5" ]]; then
  module purge; conda deactivate
  conda activate $ANVIO
cd ${output}/ANVIO_${prefix}

### get taxonomy for genecalls using Kraken2 and reformat:
anvi-get-sequences-for-gene-calls -c ${prefix}.db -o ${prefix}_gene_calls.fa

### run kaiju for RVDB
conda deactivate
conda activate $KAIJU

kaiju -t ${kaiju_rvdb}/nodes.dmp -f ${kaiju_rvdb}/*.fmi -i ${prefix}_gene_calls.fa -o ${prefix}_gene_calls.rvdb.kaiju -z $SLURM_CPUS_PER_TASK -v
sort -t $'\t' -V -k 2,2 ${prefix}_gene_calls.rvdb.kaiju -o ${prefix}_gene_calls.rvdb.kaiju

kaiju -t ${kaiju_nr_euk}/nodes.dmp -f ${kaiju_nr_euk}/*.fmi -i ${prefix}_gene_calls.fa -o ${prefix}_gene_calls.nr_euk.kaiju -z $SLURM_CPUS_PER_TASK -v
sort -t $'\t' -V -k 2,2 ${prefix}_gene_calls.nr_euk.kaiju -o ${prefix}_gene_calls.nr_euk.kaiju

kaiju-addTaxonNames -t ${kaiju_rvdb}/nodes.dmp -n ${kaiju_rvdb}/names.dmp -i ${prefix}_gene_calls.rvdb.kaiju -o ${prefix}_gene_calls.rvdb.kaiju.names -r superkingdom,phylum,order,class,family,genus,species

kaiju-addTaxonNames -t ${kaiju_nr_euk}/nodes.dmp -n ${kaiju_nr_euk}/names.dmp -i ${prefix}_gene_calls.nr_euk.kaiju -o ${prefix}_gene_calls.nr_euk.kaiju.names -r superkingdom,phylum,order,class,family,genus,species

##### prepare both kaiju runs. We remove any viral hits from the nr database because we trust RVDB more for this
awk '{print $0"\t0\tNA\tNA\tNA\tNA"}' ${prefix}_gene_calls.rvdb.kaiju.names | cut -f 1-8 > ${prefix}.kaijumerge.tmp1
awk '{print $0"\t0\tNA\tNA\tNA\tNA"}' ${prefix}_gene_calls.nr_euk.kaiju.names | cut -f 1-8 | awk -F "\t" 'BEGIN{OFS="\t"};{if($8 ~ "Viruses;") print "U",$2,"0","0","NA","NA","NA","NA"; else print $0}' > ${prefix}.kaijumerge.tmp2

#### align the two kaiju outputs, keeping whichever hit has the higher score. Kaiju seems to have similar scores even with different database sizes
paste ${prefix}.kaijumerge.tmp1 ${prefix}.kaijumerge.tmp2 | awk -F "\t" '{if($4 >= $12) print $0}' | cut -f 1-8 > ${prefix}.kaijumerge.tmp3
paste ${prefix}.kaijumerge.tmp1 ${prefix}.kaijumerge.tmp2 | awk -F "\t" '{if($4 < $12) print $0}' | cut -f 9-16 > ${prefix}.kaijumerge.tmp4

cat ${prefix}.kaijumerge.tmp3 ${prefix}.kaijumerge.tmp4 | sort -t $'\t' -V -k 2,2 -o ${prefix}_gene_calls.rvdb_nr_euk.merged.kaiju

sed -i 's/0\tNA\tNA\tNA\tNA//g' ${prefix}_gene_calls.rvdb_nr_euk.merged.kaiju

rm -f *.kaijumerge.tmp*

conda deactivate
conda activate $ANVIO

anvi-import-taxonomy-for-genes -c ${prefix}.db -i ${prefix}_gene_calls.rvdb_nr_euk.merged.kaiju -p kaiju --just-do-it

##############################################################
fi
##############################################################


#################### STEP A6 - mapping reads vs contigs to create bam files [skipped if nomap] #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A6" ]] && [[ ! "$step" =~ "nomap" ]]; then
  module purge; conda deactivate
  module load minimap2/2.17
  module load samtools
cd ${output}/ANVIO_${prefix}

#### make BAM file of mapping reads vs contigs and sort and index it
### the same concatenate contig fasta must be used for all read files

for SETS in $(cat $input)
do
    READS=$(echo $SETS | cut -d ";" -f 2)
    PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
    minimap2 -t $SLURM_CPUS_PER_TASK -ax map-ont ${prefix}.fa $READS | samtools sort -O BAM - > ${PROFILENAME}.bam
    samtools index -@ $SLURM_CPUS_PER_TASK -b ${PROFILENAME}.bam
done

##############################################################
fi
##############################################################

#################### STEP A7 - create profile for each sample [resume nomap] #################################

##############################################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A7" ]]; then
##############################################################
  module purge; conda deactivate
  conda activate $ANVIO
  cd ${output}/ANVIO_${prefix}

####### with mapping, first attempt
for SETS in $(cat $input)
do
    PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
    if [[ ! "$step" =~ "nomap" ]] && [[ ! "$step" =~ "resume" ]]; then
        anvi-profile -i ${PROFILENAME}.bam -c ${prefix}.db -o profile_${PROFILENAME} -S "${PROFILENAME}" --min-coverage-for-variability 10 -T $SLURM_CPUS_PER_TASK -M $mincontigsize -W
    fi

####### with NO mapping, first attempt
    if [[ "$step" =~ "nomap" ]] && [[ ! "$step" =~ "resume" ]]; then
        anvi-profile -i ${PROFILENAME}.bam -c ${prefix}.db -o profile_${PROFILENAME} -S "${PROFILENAME}" --min-coverage-for-variability 10 -T $SLURM_CPUS_PER_TASK -M $mincontigsize -W --blank-profile
    fi
####### with mapping, resume

    if [[ ! "$step" =~ "nomap" ]] && [[ "$step" =~ "resume" ]]; then
        anvi-profile -i ${PROFILENAME}.bam -c ${prefix}.db -o profile_${PROFILENAME} -S "${PROFILENAME}" --min-coverage-for-variability 10 -T $SLURM_CPUS_PER_TASK -M $mincontigsize 
    fi

####### with NO mapping, resume
    if [[ "$step" =~ "nomap" ]] && [[ "$step" =~ "resume" ]]; then
        anvi-profile -i ${PROFILENAME}.bam -c ${prefix}.db -o profile_${PROFILENAME} -S "${PROFILENAME}" --min-coverage-for-variability 10 -T $SLURM_CPUS_PER_TASK -M $mincontigsize --blank-profile
    fi

done

##############################################################
fi
##############################################################


#################### STEP A8 - Diamond blastx whole contigs against VOGDB (optional) #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A8" ]]; then
cd ${output}/ANVIO_${prefix}

module purge; conda deactivate
conda activate $ANVIO

anvi-export-contigs --splits-mode -c ${prefix}.db -o ${prefix}-SPLITS.fa

### Diamond blastx whole splits against custom dmnd databases
### reduce to best hit per split and VOG hit only (by evalue)

MAXTARGETS=$($SEQKIT stat -T ${prefix}-SPLITS.fa | sed '1,1d' | awk -F "\t" '{printf "%0.0f\n", $8/100}')

for f in $(realpath $CUSTOM_DMND/*); do
    $DIAMOND blastx --max-target-seqs $MAXTARGETS --evalue 1e-180 --sensitive -p $SLURM_CPUS_PER_TASK -d $f -q ${prefix}-SPLITS.fa -f 6 -o $(basename $f).evalue_min.dmnd.blastx
    $DIAMOND blastx --max-target-seqs $MAXTARGETS --evalue 1e-60 --sensitive -p $SLURM_CPUS_PER_TASK -d $f -q ${prefix}-SPLITS.fa -f 6 -o $(basename $f).evalue_med.dmnd.blastx
    $DIAMOND blastx --max-target-seqs $MAXTARGETS --evalue 1e-20 --sensitive -p $SLURM_CPUS_PER_TASK -d $f -q ${prefix}-SPLITS.fa -f 6 -o $(basename $f).evalue_max.dmnd.blastx
done

for f in *.evalue_m*.dmnd.blastx; do
    awk -F "\t" '!a[$1,$2]++' $f > ${f}.tmp1
done


### reduce info for each split to split name, number of vog hits, evalues comma separated, VOG codes comma separated

for f in *.evalue_m*.dmnd.blastx.tmp1; do
    cut -f 1 $f | sort -V | uniq -c | sed -E 's/ +([0-9]+) (.*)/\2\t\1/g' > $(basename $f .tmp1).txt
    sed -i -z "s/^/split\t$f\n/1" $(basename $f .tmp1).txt
    sed -i -E 's/\t(.+evalue_...)\.dmnd\.blastx\.tmp1/\t\1/g' $(basename $f .tmp1).txt
done

for SETS in $(cat $input)
do
    PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
        for f in *.evalue_m*.dmnd.blastx.txt; do
            anvi-import-misc-data -p profile_${PROFILENAME}/PROFILE.db --target-data-table items --just-do-it ${f}
        done
done


rm -f *.evalue_m*.dmnd.blastx.tmp1


##############################################################
fi
##############################################################



##############################################################
#################### STEP A9 - merge the profiles together [if nomap, each db is finalised separately] #################################
##############################################################

if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A9" ]]; then
  module purge; conda deactivate
  conda activate $ANVIO
cd ${output}/ANVIO_${prefix}

    if [[ ! "$step" =~ "nomap" ]]; then

        PROFILELIST=$(cut -d ";" -f 3 ${input} | sed -E 's/(.+)/profile_\1\/PROFILE.db /g' | tr -d '\n' | sed 's/$/\n/g')
        rm -fr profile_${prefix}-MERGED
        ##### use --enforce-hierarchical-clustering if your profile has >20000 splits
        anvi-merge ${PROFILELIST} -o profile_${prefix}-MERGED -c ${prefix}.db -W
        a="${output}/ANVIO_${prefix}/profile_${prefix}-MERGED/PROFILE.db"
        b=" -V ${output}/ANVIO_${prefix}/${prefix}-RPKM.txt"


            if [[ ! "$step" =~ "A8" ]]; then

                for f in *.evalue_m*.dmnd.blastx.txt; do
                    anvi-import-misc-data -p profile_${prefix}-MERGED/PROFILE.db --target-data-table items --just-do-it ${f}
                done
            fi

    fi


    if [[ "$step" =~ "nomap" ]]; then
        a="${output}/ANVIO_${prefix}/profile_${PROFILENAME}/PROFILE.db"
        b=""
    fi

#### prepare command for interactive interface and visualise
### NOT for job submission or batch; interactive only
## this generates a file that can be copy pasted from or use bash command on to run visualisation step in terminal
echo -e "\n##### \n!!! To run your interactive interface !!! \n#######\nconda activate $ANVIO\n##############\nanvi-interactive -p ${a} -c ${output}/ANVIO_${prefix}/${prefix}.db --browser-path /usr/bin/chromium-browser${b}\n\n !!! OR (for much better performance) !!!\n\nIF logged in with:\nssh -L 8080:localhost:8080\ndo\nanvi-interactive -p ${a} -c ${output}/ANVIO_${prefix}/${prefix}.db --server-only" > ${prefix}_ANVIO_INTERACTIVE_GUIDE

echo -e "\n##############\nAfter MANUAL binning done, save bins as collection (named manual) via interactive interface, then do:\nconda activate $ANVIO; anvi-summarize -p ${a} -c ${output}/ANVIO_${prefix}/${prefix}.db -o ${output}/ANVIO_${prefix}/${prefix}-SUMMARY -C manual\nIt would also be good to save an .svg of this figure with all bins visible\n\nAny images saved from interactive interface will be sent to ~/Downloads" >> ${prefix}_ANVIO_INTERACTIVE_GUIDE

echo -e "\n##############\nFor an interactive summary of results:\n/usr/bin/chromium-browser ${output}/ANVIO_${prefix}/${prefix}-SUMMARY/index.html" >> ${prefix}_ANVIO_INTERACTIVE_GUIDE

echo -e "\n##############\nAfter SUMMARY done, note which bins need refined (i.e. which bins are contaminated?), then for each bad bin, do:\nanvi-refine -p ${a} -c ${output}/ANVIO_${prefix}/${prefix}.db --browser-path /usr/bin/chromium-browser -C manual -b [BIN_NAME]" >> ${prefix}_ANVIO_INTERACTIVE_GUIDE


echo -e "\n###########\nTo add additional info to the layers, visible in the layer plots: !!! \ncreate a tsv file e.g.:\nsamples density varroa_per_100_bees\nSU_TX_22_6083   0.0000000935    1.333333333\nSU_TX_22_6084   0.0000000935    6.333333333\nEQ_TX_22_5997   0.000000194     0\n\nthen do\nanvi-import-misc-data ${prefix}-additional_info.txt -p ${a} --target-data-table layers\n" >> ${prefix}_ANVIO_INTERACTIVE_GUIDE

echo -e "\n###########\nTo reorder the layers: !!!\n create a tsv file with the data name, type, and sample names listed in the desired order e.g.:\nitem_name   data_type   data_value\nphylogeny   newick  (c2:0.0370199,(c1:0.0227268,c3:0.0227268)Int3:0.0370199)\ndensity   basic   c3,c2,c1\n\nthen do\nanvi-import-misc-data ${output}/ANVIO_${prefix}/${prefix}-layer_orders.txt -p ${a} --target-data-table layer_orders\n" >> ${prefix}_ANVIO_INTERACTIVE_GUIDE

##############################################################
fi
##############################################################

##################################################################
#################### STEP B1 - prepare normalised counts of coverage (RPKM) (optional) [skipped if nomap] #################################

#### prepare matrix of split coverage for normalisation using python script
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "B1" ]] && [[ ! "$step" =~ "nomap" ]]; then
  cd ${output}/ANVIO_${prefix}

  module purge; conda deactivate
  conda activate $ANVIO

anvi-export-splits-and-coverages --splits-mode -p profile_${prefix}-MERGED/PROFILE.db -c ${prefix}.db -o . -O ${prefix}

$SEQKIT fx2tab -nl ${prefix}-SPLITS.fa | cut -f 2 | sed -z 's/^/length\n/g' | paste ${prefix}-COVs.txt - > tmp_${prefix} && mv tmp_${prefix} ${prefix}-COVs.txt

$COVERAGE2RPKM -i ${prefix}-COVs.txt -o ${prefix}-RPKM.txt

cut -f 1 ${prefix}-RPKM.txt | sed 's/_split.*//g' | sed 's/^contig$/__parent__/g' | paste ${prefix}-RPKM.txt - > tmp_${prefix} && mv tmp_${prefix} ${prefix}-RPKM.txt


##############################################################
fi
##############################################################
##################################################################

############# enable permission for all in group ##################
chmod 777 -R ${output}/ANVIO_${prefix}
