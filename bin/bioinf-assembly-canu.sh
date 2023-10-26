#!/bin/bash

umask 002
green='\033[32m';red='\033[31m';cyan='\033[36m';nocolor='\033[m'
echo -e "${green}RUNNING bioinf-assembly-canu with the following parameters: $@ ${nocolor}"

Help()
{
echo -e "bioinf-assembly-canu is for de novo assembly of ONT sequencing data of metagenomes, using CANU. \nAll reads provided as input will be assembled together, so provide only the reads for a single sequencing barcode per script job submission (unless you intend to co-assemble multiple samples/runs).\nThis script will NOT perform any changes to your original input data.\nIt will use Porechop to remove Oxford Nanopore adapters.\n${red}This script is intended for job submission to a SLURM cluster using sbatch - running it another way will NOT work.${nocolor}\n"
echo -e "${green}########### REQUIRED arguments ###########${nocolor}\n"
echo -e "-i -input\t Provide the path to a single fastq.gz file OR a single directory containing multiple fastq.gz files. OR provide a plain text file that lists the absolute paths to multiple read FILES in a single column (this is useful for assembling together reads that are within multiple directories.\n${red}Read files MUST be in the fastq.gz format AND have the .fastq.gz file extension${nocolor}\n"
echo -e "-p -project\tA meaningful and unique name for your project. Output files will contain this name\n"
echo -e "${green}########### OPTIONAL arguments ###########${nocolor}\n"
echo -e "-s, -step\tWhich script steps to run. Steps available:\n"
awk '/^###### STEP-/ {print "\t\t"$0}' $(which bioinf-assembly-canu.sh)
echo -e "\n-G, -genomesize [default = 10000]\texpected genome size in bp - use the smallest genome size expected (going too small is fine, but a genome size that is too big may cause assembly failure)\n"
echo -e "-R, -minread [default = 1000]\tminimum read length to use in bp - use a read length that is short enough so that CANU doesn't exclude too much of your data\n"
echo -e "-V, -minoverlap [default = 500]\tminimum read overlaps to allow in bp - CANU's default is 500 bp overlap with a minimum read length of 1000 bp - you might need to try proportionately shorter overlaps than this, especially if your reads are shorter than 1000 bp. For 300-500 bp long reads, I have found that an overlap of 20-50 succeeds at generating good contigs, but going higher results in few or no contigs\n"
echo -e "-r, -resume [default = false]\trestart a previous run from where it last ended. Useful for resuming interrupted jobs without losing progress\n"
echo -e "-I, -isoncorrect [default = false]\tPerform additional correction of reads using isONcorrect. For cDNA data\n"
echo -e "!!! Checking your ASSEMBLY run !!!\n check your slurm.out and slurm.err logs and the canu.out file in your final output folder. OR check what files you have in your temporary and final output folders - a successful run will have corrected reads, trimmed reads, AND contigs.fasta\n this script will also generate a series of plots to summarise the taxa present in your dataset and the fate of the reads"
echo -e "\n!!! IMPORTANT NOTES !!!\n There is some error that causes canu jobs to end early without generating contigs. This usually happens during the final contig generation stage. It is possibly a memory overuse issue, BUT the job might complete if you resume it. If your job has trimmed reads file and corrected reads file, but no contigs, then run it again with --resume. In fact it is probably best to just try and resume all failed canu runs anyway"
echo -e "\nCANU has a complex set of options, and this script in its original state has certain set parameters for canu, including: -nanopore, 32g minimum memory, use grid. If you wish to change canus' parameters, the only way to do so is to edit your copy of the script itself. Feel free to do so"
echo -e "\nEXAMPLE OF RUNNING SCRIPT (SLURM sbatch):\n${cyan}\tsbatch --time=12:00:00 --cpus-per-task=12 --mem=64GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err bioinf-assembly-canu.sh -i /panfs/jay/groups/27/dcschroe/shared/data/gridion_dmckeown/02FEB23DM1/02FEB23DM1/20230202_1617_X1_FAU43675_8c3f2c30/fastq_pass/barcode07 -p 02FEB23DM1_barcode07 -G 5000 -R 200 -V 25${nocolor}\n"
echo -e "-h, --help\tshow this help message and exit\n"
}

### set euclidean options to default setting here
resume="false"
isoncorrect="false"
###

while getopts i:p:s:G:R:V:Irh option
do 
    case "${option}" in
        i)input=${OPTARG};;
        p)project=${OPTARG};;

        s)step=${OPTARG};;
    G)genomesize=${OPTARG};;
    R)minread=${OPTARG};;
    V)minoverlap=${OPTARG};;
    r)resume="true";;
    I)isoncorrect="true";;
    h)Help; exit;;
    esac
done


####################### SOFTWARE #####################################
### LOAD software available via shared environment on server:
module purge
eval "$(conda shell.bash hook)"
conda activate bioinftools

### scripts
bioinf_read_check_py=$(which bioinf-read-check.py)


####################### SET AND CHECK VARIABLES/ARGUMENTS #####################################
mkdir -p ${bioinftmp} ## make tmp data storage in scratch
OUTDIR="${PWD}/ASSEMBLY-CANU_${project}"
TMPDIR="${bioinftmp}/ASSEMBLY-CANU_${project}"
echo -e "${cyan}Your temporary outputs for this run will be at: ${TMPDIR}\n Your final outputs for this run will be at:${OUTDIR}${nocolor}\n"

### Stupidly, SLURM does not have env variables for these, so we must set them here
SLURM_TIME=$(squeue -j $SLURM_JOB_ID -h --Format TimeLimit)
SLURM_PARTITION=$(squeue -j $SLURM_JOB_ID -h --Format Partition)

######### from commandline arugments ######
if [[ -z "${input}" ]]; then echo -e "${red}-i, --input missing"; Help; exit; fi
if [[ -z "${project}" ]]; then echo -e "${red}-p, --project missing"; Help; exit; fi

### set defaults for optional arguments
if [[ -z "${genomesize}" ]]; then genomesize="10000"; echo -e "${green}Bioinf-assembly-canu is using default genome length 10000 bp ${nocolor}"; fi
if [[ -z "${minread}" ]]; then minread="1000"; echo -e "${green}Bioinf-assembly-canu is using default read length 1000 bp ${nocolor}"; fi
if [[ -z "${minoverlap}" ]]; then minoverlap="500"; echo -e "${green}Bioinf-assembly-canu is using default read overlap length 500 bp ${nocolor}"; fi

### THREADS

if [[ ! -z "${SLURM_CPUS_PER_TASK}" ]]; then THREADS="${SLURM_CPUS_PER_TASK}"; echo -e "${green}Using threads set by SLURM_CPUS_PER_TASK:${nocolor}"; fi
if [[ -z "${SLURM_CPUS_PER_TASK}" ]]; then THREADS="4"; echo -e "${green}No SLURM_CPUS_PER_TASK or --threads set, using default threads:${nocolor}"; fi
echo -e "\t${THREADS} threads"

#######
if [[ -z $(conda env list | grep "bioinftools") ]]; then echo "NO conda environment for bioinftools found - see README"; else echo "conda environment for bioinftools FOUND"; fi

if [ -z ${bioinftmp+x} ]; then echo -e "${red}bioinf temporary working space is unset - see README and bioinf-setup.sh in bioinf directory ${nocolor}"; exit; fi
echo -e "${green}bioinf temporary working space is set to '$bioinftmp' ${nocolor}"

####################### DATABASES #####################################

################################################################################################
####################### RUNNING THE SCRIPT #####################################################
################################################################################################

###### STEP-A1 locate read files, trim ONT barcodes with porechop, summarise read stats

################################################################################################
if [[ -z "${step}" ]] || [[ "$step" == "A1" ]] && [[ "$resume" == "false" ]]; then
################################################################################################

### single/multiple fastq.gz file(s) OR single/multiple directories containing fastq.gz file(s)
rm -fr ${TMPDIR}
mkdir ${TMPDIR}
rm -fr ${OUTDIR}
cd ${TMPDIR}

### if input is a directory, then find fastq.gz within it
### if input is a file, and its name includes "fastq.gz", then use it as input
### if the input is a file, but its name does not include "fastq.gz", and it contains "fastq.gz", then use the list of reads within it as input
if
        [[ -d ${input} ]]; then find $input -name "*.fastq.gz" > ${TMPDIR}/${project}.bctrimmedreads.list
    elif
        [[ -f ${input} ]] && [[ "$input" =~ ".fastq.gz" ]]; then echo "${input}" > ${TMPDIR}/${project}.bctrimmedreads.list
    elif
        [[ -f ${input} ]] && [[ ! "$input" =~ ".fastq.gz" ]] && [[ $(cat "$input") =~ ".fastq.gz" ]]; then cat "${input}" > ${TMPDIR}/${project}.bctrimmedreads.list
fi

### quit if no read fast.gz found
if [[ ! $(cat "${TMPDIR}/${project}.bctrimmedreads.list") =~ ".fastq.gz" ]]; then echo -e "${red}\t\tFAILED to find read files. Check your --input. Are your files fastq.gz, both in name and format?${nocolor}"; exit; fi

### concatenate all of the input read files together
rm -f ${TMPDIR}/${project}.bctrimmedreads.fastq.tmp; touch ${TMPDIR}/${project}.bctrimmedreads.fastq.tmp
for f in $(cat ${TMPDIR}/${project}.bctrimmedreads.list)
do
    zcat $f >>  ${TMPDIR}/${project}.bctrimmedreads.fastq.tmp
done


###### STORE whole deflines from fastq
awk '/^@/' ${TMPDIR}/${project}.bctrimmedreads.fastq.tmp > ${TMPDIR}/${project}.originalreads.deflines

################################################################################################
### Perform isONcorrect
if [[ "$isoncorrect" == "true" ]]; then

echo -e "${cyan}Performing isONcorrect correction on reads before assembly${nocolor}"
# Pipeline to get high-quality full-length reads from ONT cDNA sequencing

#cdna_classifier.py ${TMPDIR}/${project}.bctrimmedreads.fastq.tmp ${TMPDIR}/IOC-full_length.fq -t ${THREADS}
cp ${TMPDIR}/${project}.bctrimmedreads.fastq.tmp ${TMPDIR}/IOC-full_length.fq ### lazy way of skipping the cdna_classifier
isONclust --t ${THREADS} --ont --fastq ${TMPDIR}/IOC-full_length.fq --outfolder ${TMPDIR}/IOC-clustering

isONclust write_fastq --N 1 --clusters ${TMPDIR}/IOC-clustering/final_clusters.tsv --fastq ${TMPDIR}/IOC-full_length.fq --outfolder ${TMPDIR}/IOC-clustering/fastq_files

run_isoncorrect --t ${THREADS} --fastq_folder ${TMPDIR}/IOC-clustering/fastq_files --outfolder ${TMPDIR}/IOC-correction/

# MERGE ALL CORRECTED READS INTO ONE FILE
rm -f ${TMPDIR}/${project}.bctrimmedreads.fastq.tmp
touch ${TMPDIR}/${project}.bctrimmedreads.fastq.tmp
OUTFILES=${TMPDIR}"/IOC-correction/"*"/corrected_reads.fastq"
for f in $OUTFILES
do 
  cat $f >> ${TMPDIR}/${project}.bctrimmedreads.fastq.tmp
done

fi
################################################################################################

### TRIM the barcodes
porechop -t $THREADS -i ${TMPDIR}/${project}.bctrimmedreads.fastq.tmp -o ${TMPDIR}/${project}.bctrimmedreads.fastq
pigz -p $THREADS ${TMPDIR}/${project}.bctrimmedreads.fastq
rm -f ${TMPDIR}/${project}.bctrimmedreads.fastq.tmp

##### print read stats to slurm.out
echo -e "\n\t${green}READ STATS\n"
seqkit stat -T -j $THREADS ${TMPDIR}/${project}.bctrimmedreads.fastq.gz | sed '1,1d' | awk -F "\t" '{print "\nfile\t"$1"\nformat\t"$2"\ntype\t"$3"\nnum_seqs\t"$4"\nsum_len\t"$5"\nmin_len\t"$6"\navg_len\t"$7"\nmax_len\t"$8"\n\n"}'
echo -e "\tREAD STATS${nocolor}"

################################################################################################
fi
################################################################################################

###### STEP-A2 run the CANU assembly via submission of job array to SLURM. The script will wait for all subprocesses to finish before moving the final outputs from temporary storage to the output directory
################################################################################################
if [[ -z "${step}" ]] || [[ "$step" == "A2" ]]; then
################################################################################################
cd ${TMPDIR}

rm -f "${TMPDIR}"/"${project}"
echo "RUNNING" > "${TMPDIR}"/"${project}"

if [[ "$isoncorrect" == "true" ]]; then
    trimassemble="-trim-assemble"
elif [[ "$isoncorrect" == "false" ]]; then
    trimassemble=""
fi

canu -p "${project}" -d "${TMPDIR}" \
genomeSize="${genomesize}" \
maxInputCoverage=10000 corOutCoverage=all corMinCoverage=0 corMhapSensitivity=high \
minReadLength="${minread}" minOverlapLength="${minoverlap}" $trimassemble \
useGrid=true -nanopore "${TMPDIR}/${project}.bctrimmedreads.fastq.gz" \
gridOptionsJobName="canu.${SLURM_JOB_ID}" \
gridOptions="--time=${SLURM_TIME} --partition ${SLURM_PARTITION}" \
-maxMemory=$(($SLURM_MEM_PER_NODE -1)) \
-maxThreads=$THREADS \
-minMemory=32g \
saveOverlaps=false purgeOverlaps=aggressive stageDirectory="${bioinftmp}/canu.${SLURM_JOB_ID}" \
onSuccess="sed -i 's/RUNNING/SUCCESS/g'" \
onFailure="sed -i 's/RUNNING/FAILURE/g'"

################################################################################################
fi
################################################################################################

################################################################################################
if [[ -z "${step}" ]] || [[ "$step" == "A2" ]]; then
################################################################################################
cd ${TMPDIR}

########### the following steps check that the job is running every minute. Once the job is found to be no longer running, the outputs are moved to the final OUTDIR folder from the temporary one and the tmp files are deleted (if run is successful or if it fails due to insufficient read coverage)
#### this is done to avoid filling up the temporary storage with assembly temporary files, which can be very large
#### runs that fail due to running out of time are left as they are and can be restarted by running only step A2 again

ISJOBRUNNING=$(grep -s "RUNNING" "${TMPDIR}"/"${project}" | wc -l)

##### while the slurm jobs running list does contain your jobname (not 0), the script waits. Once the job name is not found (0), then the loops breaks and continues with the following steps

while [[ ! $ISJOBRUNNING -le 0 ]]
do
  sleep 60s
  ISJOBRUNNING=$(grep -s "RUNNING" "${TMPDIR}"/"${project}" | wc -l)

if [[ $ISJOBRUNNING -le 0 ]]; then
    break
fi

done

##### make the final output directory
rm -fr ${OUTDIR}
mkdir ${OUTDIR}

#### read stats to canu.out
echo -e "\nREAD STATS:\n\n" >> ${TMPDIR}/canu.out
seqkit stat -T -j $THREADS ${TMPDIR}/${project}.bctrimmedreads.fastq.gz | sed '1,1d' | awk -F "\t" '{print "\nfile\t"$1"\nformat\t"$2"\ntype\t"$3"\nnum_seqs\t"$4"\nsum_len\t"$5"\nmin_len\t"$6"\navg_len\t"$7"\nmax_len\t"$8"\n\n"}' >> ${TMPDIR}/canu.out

###### success?
JOBSUCCESS=$(grep -s "SUCCESS" "${TMPDIR}"/"${project}" | wc -l)
if [[ ! $JOBSUCCESS -le 0 ]]; then

echo -e "\nYOUR JOB WAS SUCCESSFUL\n" >> ${TMPDIR}/canu.out

##### stages of completion check
REA=$(ls ${TMPDIR} | grep -s "${project}.bctrimmedreads.fastq.gz" - | wc -l | sed -e 's/0/FAILED/g' -e 's/1/SUCCEEDED/g')
COR=$(ls ${TMPDIR} | grep -s "${project}.correctedReads.fasta.gz" - | wc -l | sed -e 's/0/FAILED/g' -e 's/1/SUCCEEDED/g')
TRI=$(ls ${TMPDIR} | grep -s "${project}.trimmedReads.fasta.gz" - | wc -l | sed -e 's/0/FAILED/g' -e 's/1/SUCCEEDED/g')
CON=$(ls ${TMPDIR} | grep -s "${project}.contigs.fasta" - | wc -l | sed -e 's/0/FAILED/g' -e 's/1/SUCCEEDED/g')

echo -e "PROCESS SUMMARY:\n1) Pre-processing the reads? $REA\n2) Correcting the processed reads? $COR\n3) Trimming the corrected reads? $TRI\n4) Assembling the trimmed and corrected reads into contigs? $CON\n" >> ${TMPDIR}/canu.out

cp ${TMPDIR}/canu.out ${OUTDIR}/${project}.canu.out
cp ${TMPDIR}/${project}.correctedReads.fasta.gz ${OUTDIR}/
cp ${TMPDIR}/${project}.trimmedReads.fasta.gz ${OUTDIR}/
cp ${TMPDIR}/${project}.report ${OUTDIR}/
cp ${TMPDIR}/${project}.contigs.layout.tigInfo ${OUTDIR}/
cp ${TMPDIR}/${project}.contigs.layout.readToTig ${OUTDIR}/
cp ${TMPDIR}/${project}.unassembled.fasta ${OUTDIR}/
cp ${TMPDIR}/${project}.contigs.fasta ${OUTDIR}/
cp ${TMPDIR}/${project}.bctrimmedreads.fastq.gz ${OUTDIR}/
cp ${TMPDIR}/${project}.bctrimmedreads.list ${OUTDIR}/
cp ${TMPDIR}/${project}.bctrimmedreads.fastq.wholedeflines ${OUTDIR}/

seqkit seq --min-len 200 --min-qual 9 ${TMPDIR}/${project}.bctrimmedreads.fastq.gz > ${OUTDIR}/${project}.MinLen200-minQ9.fastq
seqkit stat -a -b -T ${OUTDIR}/${project}*.fast* | sed '1,1d' | sed -E -e 's/\.fastq|\.fasta|\.gz//g' -e 's/\t/;/g' | sed "s/${project}./${project};/g" > ${OUTDIR}/${project}.stats.csv
seqkit stat -a -b -T ${OUTDIR}/${project}.stats.csv | head -1 | sed 's/\t/;/g' | cat - ${OUTDIR}/${project}.stats.csv > ${OUTDIR}/${project}.stats.csv.tmp && mv ${OUTDIR}/${project}.stats.csv.tmp ${OUTDIR}/${project}.stats.csv
rm -fr ${OUTDIR}/${project}.MinLen200-minQ9.fastq

cd ${OUTDIR}

rm -fr ${TMPDIR}

#### create file that shows which reads got assembled
# first get a key to convert canu's post-correction names for the reads (readID) to the original read names
zcat ${OUTDIR}/${project}.correctedReads.fasta.gz | grep ">" | sed -E 's/>(.+) id=(.+)/\1\t\2/g' | awk '{print "s/^"$2"$/"$1"/g"}' > ${OUTDIR}/canuid-origid.tmp

# now convert the readIDs back to the original IDs so that we can identify them
cut -f 1 ${OUTDIR}/${project}.contigs.layout.readToTig | sed -f ${OUTDIR}/canuid-origid.tmp - | paste - ${OUTDIR}/${project}.contigs.layout.readToTig | cut -f 1,3 | sed '1,1d' > ${OUTDIR}/origid-canutigid.tmp

# reduce the previous key down to only the lines that are relevant to the contigs assembled
grep ">" ${OUTDIR}/${project}.contigs.fasta | sed -E 's/>tig([0-9]+).*/\1/g' | sed -E 's/^0+//g' | grep -w -f - ${OUTDIR}/origid-canutigid.tmp > ${OUTDIR}/origid-canutigid-assembled.tmp

grep ">" ${OUTDIR}/${project}.contigs.fasta | sed -E 's/>tig([0-9]+).*/\1/g' | sed -E 's/^(0+)([1-9]+[0-9]*)/\2\ttig\1\2/g' | awk '{print "s/\\t"$1"$/\\t"$2"/g"}' | sed -f - ${OUTDIR}/origid-canutigid-assembled.tmp > ${OUTDIR}/${project}.contigs.layout.readToTig.assembledonly.origNames

rm -f *canu*id*.tmp

### put together other info for R scripts
cut -f 2 ${OUTDIR}/${project}.contigs.layout.readToTig.assembledonly.origNames | sed -E 's/tig0+//g' > tmp.tig.1

sed '1,1d' ${OUTDIR}/${project}.contigs.layout.tigInfo | awk -F "\t" -v N=$(basename ${OUTDIR}/${project}.contigs.layout.tigInfo) '{print "/^"$1"$/ s/^"$1"$/"$2"\\t"$3"\\t"N"/g"}' > tmp.tig.2

sed -f tmp.tig.2 tmp.tig.1 > tmp.tig.3; paste ${OUTDIR}/${project}.contigs.layout.readToTig.assembledonly.origNames tmp.tig.3 | sed -z 's/^/read\tcontig\tcontig_length\tcontig_coverage\tsample\n/g' > ${OUTDIR}/${project}.contigs.layout.readToTig.assembledonly.origNames.lengthcov

rm -f tmp.tig.1 tmp.tig.2 tmp.tig.3
#### generate plots

python $bioinf_read_check_py -c ${OUTDIR}

exit
fi

####### failed?
JOBFAIL=$(grep -s "FAILURE" "${TMPDIR}"/"${project}" | wc -l)
if [[ ! $JOBFAIL -le 0 ]]; then

echo -e "\nYOUR JOB HAS FAILED. \nTO RESTART A FAILED JOB, RESTART WITH -s A2\nIF STEPS 1 AND/OR 2 (SEE BELOW) FAILED, THEN THERE IS PROBABLY SOMETHING WRONG WITH YOUR DATA OR INPUT COMMAND OR PARAMETERS. IF STEPS 3 AND/OR 4 FAILED, THEN TRY RESTARTING THE JOB WITH -s A2" >> ${TMPDIR}/canu.out

##### stages of completion check
REA=$(ls ${TMPDIR} | grep -s "${project}.bctrimmedreads.fastq.gz" - | wc -l | sed -e 's/0/FAILED/g' -e 's/1/SUCCEEDED/g')
COR=$(ls ${TMPDIR} | grep -s "${project}.correctedReads.fasta.gz" - | wc -l | sed -e 's/0/FAILED/g' -e 's/1/SUCCEEDED/g')
TRI=$(ls ${TMPDIR} | grep -s "${project}.trimmedReads.fasta.gz" - | wc -l | sed -e 's/0/FAILED/g' -e 's/1/SUCCEEDED/g')
CON=$(ls ${TMPDIR} | grep -s "${project}.contigs.fasta" - | wc -l | sed -e 's/0/FAILED/g' -e 's/1/SUCCEEDED/g')

echo -e "PROCESS SUMMARY:\n1) Pre-processing the reads? $REA\n2) Correcting the processed reads? $COR\n3) Trimming the corrected reads? $TRI\n4) Assembling the trimmed and corrected reads into contigs? $CON\n" >> ${TMPDIR}/canu.out


cp ${TMPDIR}/canu.out ${OUTDIR}/${project}.canu.out
cp ${TMPDIR}/${project}.correctedReads.fasta.gz ${OUTDIR}/
cp ${TMPDIR}/${project}.trimmedReads.fasta.gz ${OUTDIR}/
cp ${TMPDIR}/${project}.report ${OUTDIR}/
cp ${TMPDIR}/${project}.contigs.layout.tigInfo ${OUTDIR}/
cp ${TMPDIR}/${project}.contigs.layout.readToTig ${OUTDIR}/
cp ${TMPDIR}/${project}.unassembled.fasta ${OUTDIR}/
cp ${TMPDIR}/${project}.contigs.fasta ${OUTDIR}/
cp ${TMPDIR}/${project}.bctrimmedreads.fastq.gz ${OUTDIR}/
cp ${TMPDIR}/${project}.bctrimmedreads.list ${OUTDIR}/
cp ${TMPDIR}/${project}.bctrimmedreads.fastq.wholedeflines ${OUTDIR}/


seqkit seq --min-len 200 --min-qual 9 ${TMPDIR}/${project}.bctrimmedreads.fastq.gz > ${OUTDIR}/${project}.MinLen200-minQ9.fastq
seqkit stat -a -b -T ${OUTDIR}/${project}*.fast* | sed '1,1d' | sed -E -e 's/\.fastq|\.fasta|\.gz//g' -e 's/\t/;/g' | sed "s/${project}./${project};/g" > ${OUTDIR}/${project}.stats.csv
seqkit stat -a -b -T ${OUTDIR}/${project}*.fast* | head -1 | sed 's/\t/;/g' | cat - ${OUTDIR}/${project}.stats.csv > ${OUTDIR}/${project}.stats.csv.tmp && mv ${OUTDIR}/${project}.stats.csv.tmp ${OUTDIR}/${project}.stats.csv
rm -fr ${OUTDIR}/${project}.MinLen200-minQ9.fastq

exit
fi

################################################################################################
fi
################################################################################################
