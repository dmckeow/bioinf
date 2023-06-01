#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=600GB
#SBATCH --partition ag2tb
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

echo -e "COMMAND LINE JOB SUBMISSSION:\n\tsbatch /home/dcschroe/dmckeow/projects/DWV/script/anvio.sh $@"

Help()
{
echo -e "########### REQUIRED arguments ###########\n-i --input\tfull path to a file that lists input file locations, semi-colon delimited, 4 fields:\n"
echo -e "\tFIELD1=contig fasta, separate line per assembly (full path). Can be individual assemblies or co-assemblies."
echo -e "\tFIELD2=reads fq.gz, corresponding to each assembly contig fasta (full path)"
echo -e "\tFIELD3=group name; samples sharing the same group name will be included in the same analyses and anvio figures"
echo -e "\tFIELD4=profile name; a short name for each profile which will appear in all outputs and figures."
echo -e "!!! group and profile names MUST NOT begin with a digit, and MUST ONLY include letters, digits, and underscores and MUST be UNIQUE within your samples being run!!!\n"
echo -e "! example file for --input - in this example, there will be separate analyses and interactive figures for projectA and projectB:\n\n/home/data/19AUG22_1_DM_barcode01.contigs.fasta;/home/data/19AUG22_1_DM_barcode01-trimmedbcs.fastq.gz;projectA;hive1\n/home/data/19AUG22_1_DM_barcode03.contigs.fasta;/home/data/19AUG22_1_DM_barcode03-trimmedbcs.fastq.gz;projectA;hive2\n/home/data/19AUG22_1_DM_barcode05.contigs.fasta;/home/data/19AUG22_1_DM_barcode05-trimmedbcs.fastq.gz;projectB;hive1\n/home/data/19AUG22_1_DM_barcode06.contigs.fasta;/home/data/19AUG22_1_DM_barcode06-trimmedbcs.fastq.gz;projectB;hive2\n\n"
echo -e "-s --step\tWhich script steps to run. A0 = run the whole script. \nMultiple steps can be specified, e.g.: A1_A2_A3 will only run steps A1-A3.\n\nIf running WITHOUT mapping - i.e. you only have contigs, but no reads, then add nomap to -s , e.g. A0_nomap\n\nALSO, some steps can be resumed without deleting previous progress for that step (see below), by adding the flag resume to -s but if you resume I recommend that you delete the last output file of that step, in case it is an incomplete output file\n\n!!! STEPS AVAILABLE: !!!"
awk '/#+ STEP /' /home/dcschroe/dmckeow/projects/DWV/script/anvio.sh | sed -E -e 's/#+ *//g' -e 's/^STEP /> /g'
echo -e "\n-m --mincontigsize\tinteger to set minimum size of contigs in bp to include in analyses (don't go below 1000)"
echo -e "\n-o --output\tabsolute path for final output directory that ALREADY EXISTS. An output folder named after ANVIO_(prefix) will be created here\n"
echo -e "\n-p --prefix\ta meaningful name for output folder that must NOT ALREADY EXIST. It should be a UNIQUE name that cannot possibly be identical to any other folder within your output destination e.g. -p sequencingrun11OCT22 will make the folder ANVIO_sequencingrun11OCT22\n"
echo -e "\n! example job submission (to run the whole script with mapping):\nsbatch /home/dcschroe/dmckeow/projects/DWV/script/anvio.sh -i /home/data/anvio_JAN23_PHEX_SAMPLES1.txt -s A0 -m 2500 -o /home/data -p JAN23_PHEX\n"
echo -e "\n! example job submission (to run only steps 1-3 and 8 without mapping, and to resume any steps that can be resumed):\nsbatch /home/dcschroe/dmckeow/projects/DWV/script/anvio.sh -i /home/data/anvio_JAN23_PHEX_SAMPLES1.txt -s A1_A2_A3_A8_nomap_resume -m 2500 -o /home/data -p JAN23_PHEX\n"
echo -e "-h, --help\tshow this help message and exit\n"
}


while getopts i:s:m:o:p:h option
do 
    case "${option}" in
        i)input=${OPTARG};;
        s)step=${OPTARG};;
    m)mincontigsize=${OPTARG};;
    o)output=${OPTARG};;
    p)prefix=${OPTARG};;
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
###### KRAKEN2
touch /scratch.global/dcschroe/kraken2db_nt/*
KDBNT="/scratch.global/dcschroe/kraken2db_nt"; ### a DB comprised of NCBI's nt: see Kraken2 database bulding for instructions. NOTE - this will be very large (~> 1 TB), so do it in scratch
KDBNT_s3="s3://dmckeow/kraken2db_nt/" ### location for krakendb to be stored on 2nd tier storage

#### make a file somewhere that lists full paths to custom HMMs you want anvi'o run on
CUSTOM_HMMS="/panfs/jay/groups/27/dcschroe/dmckeow/customhmmslist" ### a hmmer profile made by the anvio team for NCLDVs - available online

VS2DB="/home/dcschroe/dmckeow/vs2/db"

VOGDB="/home/dcschroe/dmckeow/VOGDB"
DRAMDB="/home/dcschroe/dmckeow/DRAMDB/"
DIAMONDDB="/home/dcschroe/dmckeow/VOGDB/DBVOG.dmnd"

########################################################################
########################################################################

######### SOFTWARE/OTHER SCRIPTS THAT NEED SETUP BEFORE RUNNING:

VS2_ANVIO="/home/dcschroe/dmckeow/VirSorter2_to_Anvio/virsorter_to_anvio.py" ### a custom script made by the VS2 developers to reformat VS2 results for import into anvio
COVERAGE2RPKM="/home/dcschroe/dmckeow/projects/DWV/script/coverage_to_RPKM.py" ## custom script that converts coverage data into RPKM

##### Seqkit - must be installed locally, e.g.:
SEQKIT="/panfs/jay/groups/27/dcschroe/shared/tools/seqkit"

##### taxonkit - must be installed locally, e.g.:
TAXONKIT="/home/dcschroe/dmckeow/bin/taxonkit"
TAXONKITDB="/panfs/jay/groups/27/dcschroe/dmckeow/.taxonkit"
##### KRAKEN2 - must be installed locally
KRAKEN2="/home/dcschroe/dmckeow/kraken2local/kraken2"

##### ANVIO 7.1 or later - INSTALLED via CONDA

##### Virsorter2 - INSTALLED via CONDA

##### Diamond - installed locally
DIAMOND="/panfs/jay/groups/27/dcschroe/shared/tools/diamond"


######################################################################
######################################################################


### SEE script /home/dcschroe/dmckeow/projects/DWV/script/database_setup.sh to setup the DBs needed for this script !!!!!!!



####################### THE SCRIPT #####################################

#################### STEP A0 - runs STEPS A1-A9, B1-B2 #################################

################################################################################################
################################################################################################
########################################################################
########################################################################

#################### STEP A1 - pre-process your contig fasta(s) #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A1" ]]; then
  module purge; conda deactivate
  conda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/anvio-7.1

##### make the final output directory
rm -fr ${output}/ANVIO_${prefix}
mkdir ${output}/ANVIO_${prefix}
cd ${output}/ANVIO_${prefix}

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
echo -e "\n!!! STEP A1 DONE !!!\n"
##############################################################
fi
##############################################################

#################### STEP A2 - generate the contigs database [nomap resume] #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A2" ]]; then
  module purge; conda deactivate
  conda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/anvio-7.1

cd ${output}/ANVIO_${prefix}


####### with mapping, first attempt
    if [[ ! "$step" =~ "nomap" ]] && [[ ! "$step" =~ "resume" ]]; then
### generate the contigs database
### anvi-gen-contigs-database does not have overwrite option (-W), so rm -f db before is necessary
for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    rm -f ${GROUP1}.db
    anvi-gen-contigs-database -f ${GROUP1}.fa -o ${GROUP1}.db -n "${GROUP1}" -T $SLURM_CPUS_PER_TASK
done

echo -e "\n!!! STEP A2 (with mapping) DONE !!!\n"
    fi

####### with NO mapping, first attempt

    if [[ "$step" =~ "nomap" ]] && [[ ! "$step" =~ "resume" ]]; then

### generate the contigs database
### anvi-gen-contigs-database does not have overwrite option (-W), so rm -f db before is necessary
for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    rm -f ${GROUP1}.db
    anvi-gen-contigs-database --skip-mindful-splitting -f ${GROUP1}.fa -o ${GROUP1}.db -n "${GROUP1}" -T $SLURM_CPUS_PER_TASK
done

echo -e "\n!!! STEP A2 (without mapping) DONE !!!\n"
    fi

####### with mapping, resume
    if [[ ! "$step" =~ "nomap" ]] && [[ "$step" =~ "resume" ]]; then
### generate the contigs database
for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    anvi-gen-contigs-database -f ${GROUP1}.fa -o ${GROUP1}.db -n "${GROUP1}" -T $SLURM_CPUS_PER_TASK
done

echo -e "\n!!! STEP A2 (with mapping and resumed) DONE !!!\n"
    fi

####### with NO mapping, resume

    if [[ "$step" =~ "nomap" ]] && [[ "$step" =~ "resume" ]]; then

### generate the contigs database
for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    anvi-gen-contigs-database --skip-mindful-splitting -f ${GROUP1}.fa -o ${GROUP1}.db -n "${GROUP1}" -T $SLURM_CPUS_PER_TASK
done

echo -e "\n!!! STEP A2 (without mapping and resumed) DONE !!!\n"
    fi


##############################################################
fi
##############################################################

#################### STEP A3 - perform hmm searches vs default databases (optional) #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A3" ]]; then
  module purge; conda deactivate
  conda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/anvio-7.1
cd ${output}/ANVIO_${prefix}

#### run hmm profile, summarise contig stats, and identify COGs
### !!! both anvi-run-hmms and anvi-run-ncbi-cogs do NOT include viruses and most eukaryotes
### !!! as an alternative to the default hmm profiles of anvio, we will later run and import the results of Virsorter2

### default anvio HMM profile dbs
for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    anvi-run-hmms -c ${GROUP1}.db -T $SLURM_CPUS_PER_TASK --just-do-it
    anvi-run-ncbi-cogs -c ${GROUP1}.db -T $SLURM_CPUS_PER_TASK
    anvi-run-kegg-kofams -c ${GROUP1}.db -T $SLURM_CPUS_PER_TASK --just-do-it --skip-bitscore-heuristic
    ########anvi-display-contigs-stats ${GROUP1}.db -o ${GROUP1}.contigstats ## this causes loop to hang
done

echo -e "\n!!! STEP A3 DONE !!!\n"
##############################################################
fi
##############################################################

#################### STEP A4 - perform hmm searches vs custom databases, estimate taxonomy (optional)) #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A4" ]]; then
  module purge; conda deactivate
  conda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/anvio-7.1
cd ${output}/ANVIO_${prefix}

#### run hmm profile, summarise contig stats, and identify COGs
### !!! both anvi-run-hmms and anvi-run-ncbi-cogs do NOT include viruses and most eukaryotes
### !!! as an alternative to the default hmm profiles of anvio, we will later run and import the results of Virsorter2

### custom anvio HMM profile dbs

for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
    for f in $(cat $CUSTOM_HMMS)
    do
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    
    anvi-run-hmms -c ${GROUP1}.db -T $SLURM_CPUS_PER_TASK -H $f --hmmer-program hmmsearch --just-do-it
    
  done
done


### run single core gene taxonomy - this is built into anvio and cant be customised currently. It provides quick taxonomy for prokaryotes.
### !!! the alternative is to use Kraken2 and then import the results, which will provide gene level taxonomy
for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
  GROUP1=$(echo $SETS | cut -d ";" -f 3)
  anvi-run-scg-taxonomy -c ${GROUP1}.db -T $SLURM_CPUS_PER_TASK
  anvi-estimate-scg-taxonomy -c ${GROUP1}.db --metagenome-mode -T $SLURM_CPUS_PER_TASK
done

##############################################################
fi
##############################################################



#################### STEP A5 - get taxonomy for genes with Kraken2 and import into contigs database (optional) #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A5" ]]; then
  module purge; conda deactivate
  conda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/anvio-7.1
cd ${output}/ANVIO_${prefix}

### get taxonomy for genecalls using Kraken2 and reformat:
for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
  GROUP1=$(echo $SETS | cut -d ";" -f 3)
  anvi-get-sequences-for-gene-calls -c ${GROUP1}.db -o ${GROUP1}_gene_calls.fa
##### first run Kraken2 for microbial and some eukaryote databases
  $KRAKEN2 --db $KDBNT ${GROUP1}_gene_calls.fa --use-mpa-style --threads $SLURM_CPUS_PER_TASK --use-names --output ${GROUP1}.results.kraken2 --report ${GROUP1}.report.kraken2 --unclassified-out ${GROUP1}.unclassified.kraken2 --classified-out ${GROUP1}.classified.kraken2

done

### and reformat:
for SETS in $(cat $input | sort -Vu -t ";" -k 3,3)
do
  GROUP1=$(echo $SETS | cut -d ";" -f 3)
  
##### reformat kraken results for anvio import
cut -f 2,3 ${GROUP1}.results.kraken2 | sed -E 's/(.*)\(taxid ([0-9]+)\)/\1\t\2/g' > tmp1_${GROUP1}

$TAXONKIT reformat --data-dir $TAXONKITDB --threads $SLURM_CPUS_PER_TASK -F -I 3 tmp1_${GROUP1} | awk -F "\t" 'BEGIN{OFS="\t"};{if($2 ~ "unclassified" || $4 =="") print $1,$2,$3,";;;;;;"; else print $0}' | cut -f 1,4 | sed 's/;/\t/g' | sed -z 's/^/gene_callers_id\tt_domain\tt_phylum\tt_class\tt_order\tt_family\tt_genus\tt_species\n/1' > ${GROUP1}_input_matrix.txt

anvi-import-taxonomy-for-genes -c ${GROUP1}.db -i ${GROUP1}_input_matrix.txt -p default_matrix

done

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
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    PROFILENAME=$(echo $SETS | cut -d ";" -f 4)
    minimap2 -t $SLURM_CPUS_PER_TASK -ax map-ont ${GROUP1}.fa $READS | samtools sort -O BAM - > ${PROFILENAME}.bam
    samtools index -@ $SLURM_CPUS_PER_TASK -b ${PROFILENAME}.bam
done

##############################################################
fi
##############################################################

#################### STEP A7 - create profile for each contig database [resume nomap] #################################
#### profile the contig databases WITH CONCATENATED FASTA OR CO_ASSEMBLIES ONLY - MULTIPLE SEQUENCING SAMPLES (SET OF READS) AND ONE CONTIG FASTA
##############################################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A7" ]]; then
##############################################################
  module purge; conda deactivate
  conda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/anvio-7.1
  cd ${output}/ANVIO_${prefix}

####### with mapping, first attempt
    if [[ ! "$step" =~ "nomap" ]] && [[ ! "$step" =~ "resume" ]]; then
for SETS in $(cat $input)
do
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    PROFILENAME=$(echo $SETS | cut -d ";" -f 4)
    anvi-profile -i ${PROFILENAME}.bam -c ${GROUP1}.db -o profile_${PROFILENAME} -S "${PROFILENAME}" --min-coverage-for-variability 10 -T $SLURM_CPUS_PER_TASK -W
done

    fi

####### with NO mapping, first attempt

    if [[ "$step" =~ "nomap" ]] && [[ ! "$step" =~ "resume" ]]; then
for SETS in $(cat $input)
do
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    anvi-profile --blank-profile -i ${GROUP1}.bam -c ${GROUP1}.db -o profile_${GROUP1} -S "${GROUP1}" -T $SLURM_CPUS_PER_TASK -W
done

fi

####### with mapping, to resume
    if [[ ! "$step" =~ "nomap" ]] && [[ "$step" =~ "resume" ]]; then
for SETS in $(cat $input)
do
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    PROFILENAME=$(echo $SETS | cut -d ";" -f 4)
    anvi-profile -i ${PROFILENAME}.bam -c ${GROUP1}.db -o profile_${PROFILENAME} -S "${PROFILENAME}" --min-coverage-for-variability 10 -T $SLURM_CPUS_PER_TASK
done

    fi

####### with NO mapping, to resume

    if [[ "$step" =~ "nomap" ]] && [[ "$step" =~ "resume" ]]; then
for SETS in $(cat $input)
do
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    anvi-profile --blank-profile -i ${GROUP1}.bam -c ${GROUP1}.db -o profile_${GROUP1} -S "${GROUP1}" -T $SLURM_CPUS_PER_TASK
done

fi

##############################################################
fi
##############################################################




#################### STEP A8 - run Virsorter2 on gene calls (optional) #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A8" ]]; then
module purge; conda deactivate
cd ${output}/ANVIO_${prefix}

### run Virsorter2 and import results into Anvio2 (uses custom scripts made by VS2 developer - virsorter_to_anvio.py - anvi'o does not support this officially for VS2 yet)
for SETS in $(cat $input)
do
        if [[ ! "$step" =~ "nomap" ]]; then
  GROUP1=$(echo $SETS | cut -d ";" -f 3)
    PROFILENAME=$(echo $SETS | cut -d ";" -f 4)
        fi
########## GROUP name is used for both db and profile is mapping is not used
         if [[ "$step" =~ "nomap" ]]; then
  GROUP1=$(echo $SETS | cut -d ";" -f 3)
    PROFILENAME=$(echo $SETS | cut -d ";" -f 3)
        fi

  conda deactivate
  conda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/anvio-7.1
  anvi-export-table --table splits_basic_info ${GROUP1}.db -o ${GROUP1}_splits_basic_info.txt
  anvi-export-gene-calls --gene-caller prodigal -c ${GROUP1}.db -o ${GROUP1}_all_gene_calls.txt
  conda deactivate
  conda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/vs2
  rm -fr ${GROUP1}_VS2_results
  virsorter run -i ${GROUP1}.fa -w ${GROUP1}_VS2_results/ --min-length $mincontigsize --max-orf-per-seq 25 --provirus-off --keep-original-seq --prep-for-dramv --hallmark-required-on-short --viral-gene-required --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae -j $SLURM_CPUS_PER_TASK
  $VS2_ANVIO -i ${GROUP1}_VS2_results/ -s ${GROUP1}_splits_basic_info.txt -n ${GROUP1}_all_gene_calls.txt -d ${VS2DB} -A ${GROUP1}_virsorter_additional_info.txt -C ${GROUP1}_virsorter_collection.txt -F ${GROUP1}_virsorter_annotations.txt
  conda deactivate
  conda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/anvio-7.1
  anvi-import-functions -c ${GROUP1}.db -i ${GROUP1}_virsorter_annotations.txt
    anvi-import-collection ${GROUP1}_virsorter_collection.txt -c ${GROUP1}.db -p profile_${PROFILENAME}/PROFILE.db -C VS2_${GROUP1}
    anvi-import-misc-data ${GROUP1}_virsorter_additional_info.txt -p profile_${PROFILENAME}/PROFILE.db --target-data-table items

done

##############################################################
fi
##############################################################


#################### STEP A9 - Diamond blastx whole contigs against VOGDB (optional) #################################
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "A9" ]]; then
cd ${output}/ANVIO_${prefix}

module purge; conda deactivate
  conda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/anvio-7.1

  for SETS in $(cat $input)
do
  GROUP2=$(echo $SETS | cut -d ";" -f 3)

anvi-export-contigs --splits-mode -c ${GROUP2}.db -o ${GROUP2}-SPLITS.fa

### Diamond blastx whole splits against VOGDB
  $DIAMOND blastx -p $SLURM_CPUS_PER_TASK -d $DIAMONDDB -q ${GROUP2}-SPLITS.fa -f 6 -o ${GROUP2}.dmnd.blastx


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

anvi-import-misc-data ${GROUP2}_dmndblastx_additional_info.txt -p profile_${GROUP2}/PROFILE.db --target-data-table items

rm -f tmp1_${GROUP2}_*

done

##############################################################
fi
##############################################################



##############################################################
#################### STEP B1 - merge the profiles together [if nomap, each db is finalised separately] #################################
##############################################################

if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "B1" ]]; then
  module purge; conda deactivate
  conda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/anvio-7.1
cd ${output}/ANVIO_${prefix}

        if [[ ! "$step" =~ "nomap" ]]; then
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

rm -fr profile_${GROUP2}-MERGED

##### use --enforce-hierarchical-clustering if your profile has >20000 splits
anvi-merge ${PROFILELIST} -o profile_${GROUP2}-MERGED -c ${GROUP2}.db -W

#### post merge Virsorter2 import
###anvi-import-collection ${GROUP2}_virsorter_collection.txt -c ${GROUP2}.db -p profile_${GROUP2}-MERGED/PROFILE.db -C VS2_${GROUP2}
###anvi-import-misc-data ${GROUP2}_virsorter_additional_info.txt -p profile_${GROUP2}-MERGED/PROFILE.db --target-data-table items

#### prepare command for interactive interface and visualise
### NOT for job submission or batch; interactive only
## this generates a file that can be copy pasted from or use bash command on to run visualisation step in terminal
echo -e "\n##### \n!!! To run your interactive interface !!! \n#######\nconda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/anvio-7.1\n##############\nanvi-interactive -p ${output}/ANVIO_${prefix}/profile_${GROUP2}-MERGED/PROFILE.db -c ${output}/ANVIO_${prefix}/${GROUP2}.db --browser-path /usr/bin/chromium-browser -V ${output}/ANVIO_${prefix}/${GROUP2}-RPKM.txt" > ${GROUP2}_ANVIO_INTERACTIVE_GUIDE

echo -e "\n##############\nAfter MANUAL binning done, save bins as collection (named manual) via interactive interface, then do:\nconda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/anvio-7.1; anvi-summarize -p ${output}/ANVIO_${prefix}/profile_${GROUP2}-MERGED/PROFILE.db -c ${output}/ANVIO_${prefix}/${GROUP2}.db -o ${output}/ANVIO_${prefix}/${GROUP2}-SUMMARY -C manual\nIt would also be good to save an .svg of this figure with all bins visible\n\nAny images saved from interactive interface will be sent to ~/Downloads" >> ${GROUP2}_ANVIO_INTERACTIVE_GUIDE

echo -e "\n##############\nFor an interactive summary of results:\n/usr/bin/chromium-browser ${output}/ANVIO_${prefix}/${GROUP2}-SUMMARY/index.html" >> ${GROUP2}_ANVIO_INTERACTIVE_GUIDE

echo -e "\n##############\nAfter SUMMARY done, note which bins need refined (i.e. which bins are contaminated?), then for each bad bin, do:\nanvi-refine -p ${output}/ANVIO_${prefix}/profile_${GROUP2}-MERGED/PROFILE.db -c ${output}/ANVIO_${prefix}/${GROUP2}.db --browser-path /usr/bin/chromium-browser -C manual -b [BIN_NAME]" >> ${GROUP2}_ANVIO_INTERACTIVE_GUIDE


echo -e "\n###########\nTo add additional info to the layers, visible in the layer plots: !!! \ncreate a tsv file e.g.:\nsamples density varroa_per_100_bees\nSU_TX_22_6083   0.0000000935    1.333333333\nSU_TX_22_6084   0.0000000935    6.333333333\nEQ_TX_22_5997   0.000000194     0\n\nthen do\nanvi-import-misc-data ${GROUP2}-additional_info.txt -p ${output}/ANVIO_${prefix}/profile_${GROUP2}-MERGED/PROFILE.db --target-data-table layers\n" >> ${GROUP2}_ANVIO_INTERACTIVE_GUIDE

echo -e "\n###########\nTo reorder the layers: !!!\n create a tsv file with the data name, type, and sample names listed in the desired order e.g.:\nitem_name   data_type   data_value\nphylogeny   newick  (c2:0.0370199,(c1:0.0227268,c3:0.0227268)Int3:0.0370199)\ndensity   basic   c3,c2,c1\n\nthen do\nanvi-import-misc-data ${output}/ANVIO_${prefix}/${GROUP2}-layer_orders.txt -p ${output}/ANVIO_${prefix}/profile_${GROUP2}-MERGED/PROFILE.db --target-data-table layer_orders\n" >> ${GROUP2}_ANVIO_INTERACTIVE_GUIDE


done

        fi

#######################################
        if [[ "$step" =~ "nomap" ]]; then

for SETS in $(cat $input)
do
    GROUP1=$(echo $SETS | cut -d ";" -f 3)
    
#### post merge Virsorter2 import
###anvi-import-collection ${GROUP2}_virsorter_collection.txt -c ${GROUP2}.db -p profile_${GROUP2}-MERGED/PROFILE.db -C VS2_${GROUP2}
###anvi-import-misc-data ${GROUP2}_virsorter_additional_info.txt -p profile_${GROUP2}-MERGED/PROFILE.db --target-data-table items

#### prepare command for interactive interface and visualise
### NOT for job submission or batch; interactive only
## this generates a file that can be copy pasted from or use bash command on to run visualisation step in terminal
echo -e "\n##### \n!!! To run your interactive interface !!! \n#######\nconda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/anvio-7.1\nanvi-interactive -p ${output}/ANVIO_${prefix}/profile_${GROUP1}/PROFILE.db -c ${output}/ANVIO_${prefix}/${GROUP1}.db --browser-path /usr/bin/chromium-browser" > ${GROUP1}_ANVIO_INTERACTIVE_GUIDE

echo -e "\n##############\n!!! After MANUAL binning done, save bins as collection (named manual) via interactive interface, then do: !!!\nconda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/anvio-7.1; anvi-summarize -p ${output}/ANVIO_${prefix}/profile_${GROUP1}/PROFILE.db -c ${output}/ANVIO_${prefix}/${GROUP1}.db -o ${output}/ANVIO_${prefix}/${GROUP1}-SUMMARY -C manual\n\nAny images saved from interactive interface will be sent to ~/Downloads" >> ${GROUP1}_ANVIO_INTERACTIVE_GUIDE

echo -e "\n##############\n!!! For an interactive summary of results:!!! \n/usr/bin/chromium-browser ${output}/ANVIO_${prefix}/${GROUP1}-SUMMARY/index.html" >> ${GROUP1}_ANVIO_INTERACTIVE_GUIDE

echo -e "\n##############\n!!! After SUMMARY done, note which bins need refined (i.e. which bins are contaminated?), then for each bad bin, do: !!! \nanvi-refine -p ${output}/ANVIO_${prefix}/profile_${GROUP1}/PROFILE.db -c ${output}/ANVIO_${prefix}/${GROUP1}.db --browser-path /usr/bin/chromium-browser -C manual -b [BIN_NAME]" >> ${GROUP1}_ANVIO_INTERACTIVE_GUIDE


echo -e "\n###########\n !!! To add additional info to the layers, visible in the layer plots: !!!\ncreate a tsv file e.g.: \nsamples density varroa_per_100_bees\nSU_TX_22_6083   0.0000000935    1.333333333\nSU_TX_22_6084   0.0000000935    6.333333333\nEQ_TX_22_5997   0.000000194     0\n\nthen do\nanvi-import-misc-data YOUR_EXTRA_INFO_FILE -p ${output}/ANVIO_${prefix}/profile_${GROUP1}/PROFILE.db --target-data-table layers\n" >> ${GROUP1}_ANVIO_INTERACTIVE_GUIDE

echo -e "\n###########\n !!! To reorder the layers: !!! \n create a tsv file with the data name, type, and sample names listed in the desired order e.g.: \nitem_name   data_type   data_value\nphylogeny   newick  (c2:0.0370199,(c1:0.0227268,c3:0.0227268)Int3:0.0370199)\ndensity   basic   c3,c2,c1\n\nthen do\nanvi-import-misc-data YOUR_EXTRA_INFO_FILE -p ${output}/ANVIO_${prefix}/profile_${GROUP1}/PROFILE.db --target-data-table layer_orders\n" >> ${GROUP1}_ANVIO_INTERACTIVE_GUIDE


done

        fi


##############################################################
fi
##############################################################

##################################################################
#################### STEP B2 - prepare normalised counts of coverage (RPKM) (optional) [skipped if nomap] #################################

#### prepare matrix of split coverage for normalisation using python script
if [[ "$step" =~ "A0" ]] || [[ "$step" =~ "B2" ]] && [[ ! "$step" =~ "nomap" ]]; then
  cd ${output}/ANVIO_${prefix}

  module purge; conda deactivate
  conda activate /panfs/jay/groups/27/dcschroe/dmckeow/.conda/envs/anvio-7.1

  for SETS in $(cat SAMPLES2.txt)
  do
  GROUP2=$(echo $SETS | cut -d ";" -f 1)

anvi-export-splits-and-coverages --splits-mode -p profile_${GROUP2}-MERGED/PROFILE.db -c ${GROUP2}.db -o $output -O ${GROUP2}

$SEQKIT fx2tab -nl ${GROUP2}-SPLITS.fa | cut -f 2 | sed -z 's/^/length\n/g' | paste ${GROUP2}-COVs.txt - > tmp_${GROUP2} && mv tmp_${GROUP2} ${GROUP2}-COVs.txt

$COVERAGE2RPKM -i ${GROUP2}-COVs.txt -o ${GROUP2}-RPKM.txt

cut -f 1 ${GROUP2}-RPKM.txt | sed 's/_split.*//g' | sed 's/^contig$/__parent__/g' | paste ${GROUP2}-RPKM.txt - > tmp_${GROUP2} && mv tmp_${GROUP2} ${GROUP2}-RPKM.txt

done



##############################################################
fi
##############################################################
##################################################################

############# enable permission for all in group ##################
chmod 777 -R $output