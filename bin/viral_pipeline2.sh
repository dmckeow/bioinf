#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=8
#SBATCH --mem=200GB
#SBATCH --tmp=200GBB
#SBATCH --partition amdlarge
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err


####################### PREREQUISITES #####################################
## nanopore seq reads data
## rm -f slurm* tmp* ## cleanup between tests
####################### SCRIPT PURPOSE #####################################

####################### HOW TO RUN SCRIPT #####################################

####################### INPUTS (FILES, COMMAND LINE ARGUMENTS, ETC) #####################################
## $1 = input file directory (this will also be location of results)

####################### OUTPUT FILES #####################################
## in placed in directory flye at location of input file

####################### LOAD SOFTWARE #####################################
module load minimap2/2.17
module load samtools
CANU="/home/dcschroe/dcschroe/dmckeow/canu-2.2/bin/canu"
KRAKEN="/panfs/roc/msisoft/kraken/2.0.7beta/kraken2"
KDB="/panfs/roc/msisoft/kraken/kraken_db"
SEQKIT="/home/dcschroe/dcschroe/dmckeow/seqkit"
####################### SET VARIABLES/ARGUMENTS #####################################

################################################################################################
####################### RUNNING THE SCRIPT #####################################################
################################################################################################
cd $1
cd REF_"$(basename $1)"

#################### STEP 002 #################################
rm -f tmp_check; ls tmp_*.bam > tmp_check; ls tmp_*.ref1.fna.gz >> tmp_check; touch tmp_check; if [ -s tmp_check ];
then
  echo "! manually remove the tmp_ prefix from the tmp_*.bam files AND corresponding tmp_*.ref1.fna.gz files you want to keep, then do rm -f tmp* - use tmp_manual_rename for convenience!"
else
### extract reads that mapped to reference genomes

### reads previously corrected and trimmed by canu, assembly only
for f in "$(basename $1)"*.bam; do samtools fastq $f | sed -E '/^\"+$/d' > "$(basename $f .bam)".mapped.fastq; done

for f in "$(basename $1)"*.mapped.fastq; do

  GSIZE=$($SEQKIT stat -T "$(basename $f .mapped.fastq).ref1.fna.gz" | sed '1,1d' | cut -f 5) ## get reference genome size

  $CANU -p "$(basename $f .mapped.fastq)" -d "$(basename $f .mapped.fastq)" genomeSize=$GSIZE maxInputCoverage=all corOutCoverage=all corMhapSensitivity=high corMinCoverage=0 redMemory=32 oeaMemory=32 batMemory=64 -nanopore $f -assemble -trimmed -corrected
done

fi
