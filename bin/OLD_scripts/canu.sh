#!/bin/bash -l
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=150GB
#SBATCH --partition agsmall
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

####################### PREREQUISITES #####################################
## nanopore seq reads data
## rm -f slurm* tmp* ## cleanup between tests
####################### SCRIPT PURPOSE #####################################

####################### HOW TO RUN SCRIPT #####################################

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################
## $1 = input file directory (this will also be location of results)
## $2 = genome size for de novo metagenome assembly
## $3 = sequencing run ID/name
## $4 = sequencing run barcode
## $5 = output directory for all results
## $6 = final destination for storing results (without intermediate files)
## $7 = flag to control which steps are run:
## 	A0 = run all steps. If runnning on an assembly for the first time, use this option
##	A1 = run only STEP 001, which prepares raw reads for assembly AND will DELETE all previous assembly data for that run (if in the same output directories)
##	A2 = run only STEP 002, which runs the assembly. To restart a previous run that, for example, timed out on the server, use this option
## e.g.:
## sbatch /home/dcschroe/dmckeow/projects/DWV/script/viral_pipeline0.sh /home/dcschroe/dmckeow/data/05_31_2022_Honeybee/fastq_pass/barcode49 10000 05_31_2022_Honeybee barcode49 /scratch.global/dcschroe /home/dcschroe/dmckeow/data/assembly_Albert

####################### OUTPUT FILES #####################################

####### some useful code to check output and troubleshoot
##find . -name "info.txt" | awk '/DENOVO/' | awk -F ".seqStore" '!a[$1]++' | sort -V > tmp; for f in $(cat tmp); do realpath $f; done > tmp2; for f in $(cat tmp2); do echo "$f" ; cat $f; done > CANU_read_info ## check read inputs and processing

####################### LOAD SOFTWARE #####################################
module load minimap2/2.17
module load samtools
module load porechop/0.2.4
module load pigz
CANU="/home/dcschroe/dmckeow/canu-2.2/bin/canu"
KRAKEN="/panfs/roc/msisoft/kraken/2.0.7beta/kraken2"
KDB="/panfs/roc/msisoft/kraken/kraken_db"
SEQKIT="/home/dcschroe/dmckeow/seqkit"
SCRATCH="/scratch.global/dcschroe"

####################### SET VARIABLES/ARGUMENTS #####################################
HOUSE="$HOME/data/house_ref_genomes/*.fna.gz" ## house reference genomes location
ONAME="${3}_${4}"
mkdir /scratch.global/dcschroe/CANU
################################################################################################
####################### RUNNING THE SCRIPT #####################################################
################################################################################################

################################ STEP 000 ###########################################
#### use this step to generate a convenient file to submit jobs with
## bash /home/dcschroe/dmckeow/projects/DWV/script/viral_pipeline0.sh RUNFILE 11OCT22_1_DM 5000
if [[ "$1" =~ "RUNFILE" ]]; then
DATASOURCE="/scratch.global/dcschroe/data/gridion_dmckeown" ## /home/dcschroe/shared/data
date=`date +%d_%b_%Y`
RUN="$2"
RUNNAME=`echo $RUN | sed 's/[^A-Za-z0-9_]/_/g'`
GENOMESIZE="$3"
find ${DATASOURCE}/${RUN} -path "*/fastq_pass/barcode*" -name barcode* | sed -E 's/.*\/(barcode[0-9]*)$/\1/g' > tmp_${RUN}_BARCODES
SCRATCH="/scratch.global/dcschroe/CANU"
FINALOUT="/home/dcschroe/dmckeow/data/CANU_DENOVO__${RUN}__GS${3}__${date}"
find ${DATASOURCE}/${RUN} -path "*/fastq_pass/barcode*" -name barcode* | paste -d ";" - tmp_${RUN}_BARCODES | awk -F ";" -v RN=$RUNNAME -v G=$GENOMESIZE -v S=$SCRATCH -v F=$FINALOUT '{print "sbatch /home/dcschroe/dmckeow/projects/DWV/script/viral_pipeline0.sh "$1" "G" "RN" "$2" "S" "F" A0"}' > RUN_CANU_DENOVO__${RUN}__GS${3}__${date}
rm -f tmp_*BARCODES
exit
######## FOR CO ASSEMBLY, change $1 and $4 in the RUNFILE to target a higher folder containing reads of multiples barcodes/runs, AND a name indicating co-assembly

fi

################################ STEP 001 ###########################################

############ TO RUN ON A SINGLE ONT RUN
##### read correction and filtering

######### DE NOVO ASSEMBLY - CANU (AND READ TRIMMING/CORRECTION)
###### you may need to lower genome size if Canu complains that read coverage allowance is too low
##### initial steps
#### SKIP these steps IF restarting a run, otherwise all progress from previous run will be deleted
### this is useful if your previous run timed out before completion

if [[ "$7" =~ "A0" ]] || [[ "$7" =~ "A1" ]]; then

rm -fr $5/$ONAME
mkdir $5/$ONAME
touch $5/$ONAME/tmp_allreads_${ONAME}.fastq
find $1 -name "*.fastq.gz" > $5/$ONAME/tmp_allreads_${ONAME}.list
sed -i '/\/unclassified\//d' $5/$ONAME/tmp_allreads_${ONAME}.list
for f in $(cat $5/$ONAME/tmp_allreads_${ONAME}.list)
do
zcat $f >> $5/$ONAME/tmp_allreads_${ONAME}.fastq
done

porechop -t $SLURM_CPUS_PER_TASK -i $5/$ONAME/tmp_allreads_${ONAME}.fastq -o $5/$ONAME/tmp
$SEQKIT seq -j $SLURM_CPUS_PER_TASK -m 50 $5/$ONAME/tmp > $5/$ONAME/tmp_allreads_${ONAME}.fastq
rm -f $5/$ONAME/tmp
pigz -p $SLURM_CPUS_PER_TASK $5/$ONAME/tmp_allreads_${ONAME}.fastq

fi

################################ STEP 002 ###########################################

####### Assembly Step
##### with flag A2, you will rnu only this step, which means CANU will restart previous run attempt, without losing any progress
if [[ "$7" =~ "A0" ]] || [[ "$7" =~ "A2" ]]; then

READMIN=$($SEQKIT stat -T -j $SLURM_CPUS_PER_TASK $5/$ONAME/tmp_allreads_${ONAME}.fastq.gz | cut -f 7 | sed "1,1d" | sed -E "s/(.+)/\1-49/g" | bc | sed -E "s/([0-9]+)\.[0-9]+/\1/g")

#### try to reduce overlap intermediate size (corMax), porechopped reads
$CANU -p $ONAME -d $5/$ONAME genomeSize="$2" maxInputCoverage=10000 corOutCoverage=all corMhapSensitivity=high corMinCoverage=0 minReadLength=$READMIN minOverlapLength=50 useGrid=true -nanopore "$5/$ONAME/tmp_allreads_${ONAME}.fastq.gz" gridOptionsJobName=$ONAME gridOptions="--time=96:00:00 --partition agsmall" -maxMemory=149g -maxThreads=$SLURM_CPUS_PER_TASK -minMemory=16g saveOverlaps=false purgeOverlaps=aggressive stageDirectory=$SCRATCH/$SLURM_JOBID

fi

################################ STEP 003 ###########################################
###### wait until canu is DONE before moving outputs to the final destination

if [[ "$7" =~ "B0" ]] || [[ "$7" =~ "A2" ]]; then

i1=$(grep -s "\-\- Bye\.\|\-\- ERROR:  Read coverage .* lower than allowed\|ABORT:" $5/$ONAME/canu.out | wc -l)

while [ $i1 -le 0 ]
do
  sleep 60
  i1=$(grep -s "\-\- Bye\.\|\-\- ERROR:  Read coverage .* lower than allowed\|ABORT:" $5/$ONAME/canu.out | wc -l)
if [[ "$i1" == '1' ]]; then
    break

  fi

done

i2=$(grep -s "\-\- Bye\." $5/$ONAME/canu.out | wc -l)
if [[ "$i2" == '1' ]]; then
mkdir $6
rm -f $6/${ONAME}*
cp $5/$ONAME/canu.out $6/${ONAME}.FINISHED.canu.out
cp $5/$ONAME/${ONAME}.correctedReads.fasta.gz $6/
cp $5/$ONAME/${ONAME}.trimmedReads.fasta.gz $6/
cp $5/$ONAME/${ONAME}.report $6/
cp $5/$ONAME/${ONAME}.contigs.layout.tigInfo $6/
cp $5/$ONAME/${ONAME}.contigs.layout.readToTig $6/
cp $5/$ONAME/${ONAME}.unassembled.fasta $6/
cp $5/$ONAME/${ONAME}.contigs.fasta $6/
cp $5/$ONAME/allreads_${ONAME}.stats $6/
cp $5/$ONAME/tmp_allreads_${ONAME}.fastq.gz $6/
rm -fr $5/$ONAME
fi

i3=$(grep -s "\-\- ERROR:  Read coverage .* lower than allowed" $5/$ONAME/canu.out | wc -l)
if [[ "$i3" == '1' ]]; then
mkdir $6
rm -f $6/${ONAME}*
cp $5/$ONAME/canu.out $6/${ONAME}.LOWCOVERAGE.canu.out
cp $5/$ONAME/${ONAME}.correctedReads.fasta.gz $6/
cp $5/$ONAME/${ONAME}.trimmedReads.fasta.gz $6/
cp $5/$ONAME/${ONAME}.report $6/
cp $5/$ONAME/${ONAME}.contigs.layout.tigInfo $6/
cp $5/$ONAME/${ONAME}.contigs.layout.readToTig $6/
cp $5/$ONAME/${ONAME}.unassembled.fasta $6/
cp $5/$ONAME/${ONAME}.contigs.fasta $6/
cp $5/$ONAME/allreads_${ONAME}.stats $6/
cp $5/$ONAME/tmp_allreads_${ONAME}.fastq.gz $6/
rm -fr $5/$ONAME
fi

##### runs that ABORTED should be run again with flag A2
i4=$(grep -s "ABORT:" $5/$ONAME/canu.out | wc -l)
if [[ "$i4" == '1' ]]; then
mkdir $6
cp $5/$ONAME/canu.out $6/${ONAME}.ABORTED.canu.out
cp $5/$ONAME/allreads_${ONAME}.stats $6/
fi

fi
