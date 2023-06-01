#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition long
#SBATCH --cpus-per-task 6
#SBATCH --mem 48GB

############# MANUAL STEPS AND PREREQUISITES ################

##### 1. list all genomes to be analysed, then assign them to groups of 5 or less, and assign a number to each group - these are the "run groups"
#### this assumes you have a lot of genomes (e.g. 50); you may wish to use smaller or larger run groups (you can use single genome per run)
## for your own reference, record these rungroups and the genomes assigned to them, e.g.:
## [genome_assembly_name][run_group]
## Chordaria-linearis.fa 1
## Chrysoparadoxa-australica 1
## Desmarestia-dudresnayi 1
## Desmarestia-herbacea_FEMALE 2
## Desmarestia-herbacea_MALE 2
## Dictyota-dichotoma_MALE 2
## CAUTION: do not place genomes sharing identical gene or contig names in same run group
## 2. this script will create a subdirectory per genome entry (required) and move the final .gff and .faa outputs into the parent directory
## 3. ensure inputs files for formatted correctly (see IVEX000_fix_genomes.sh)
## 4. ensure all genomes to be analysed are placed in input/IVEX000_genomes/ (see $gf and $fa) - OR any other folder in which you have read write permissions (ensure you redirect $gf and $fa) - this script will use assembly fasta even if not gzipped - it will gzip input assembly at end

### run from project folder (which contains archives input script tmp finalresult)

####################### INPUTS #####################################
###### load software ######
module load bedtools/2.29.2; module load seqkit/0.14.0
GeneMarkS="/projet/fr2424/sib/dmckeown/genemark_suite_linux_64/gmsuite/gmsn.pl" ## GeneMarkS not on cluster - can be installed locally

##### command line arguments required #####
## $1 = genome assembly fasta (nucleotide .fa.gz) ABSOLUTE PATH
## $2 = rungroup number e.g. 1


##### input files for analysis (named by run groups) #####
n1=$(echo "$1" | sed 's/.*\///g' | cut -d "." -f 1) ## set name
n2=$(echo "$1" | cut -d "." -f 1,2) ## set name2
pwd=$(pwd) ## get path to project folder
FA1="original_gms"$2"_"$n1".faa"
GF1="original_gms"$2"_"$n1".gff"
gf=$pwd"/input/IVEX000_genomes/"$n1".gff"
fa=$pwd"/input/IVEX000_genomes/"$n1"_proteins.fa"

###################### SCRIPT ##############################

### create and move into genome subdirectory
cd input/IVEX000/
mkdir "$n1"_gms;
cd "$n1"_gms;

#### run GeneMarkS with virus option
gunzip $1
$GeneMarkS --clean --virus --faa --format GFF3 --name $n1 --output $n1 --species $n1 $n2

######### fix formats of GeneMarkS outputs ######
#### rename .gff and .faa
mv $n1 $GF1
mv "$n1".faa $FA1

##### fix protein/gene names
## remove hashtag lines, empty lines, and gene/CDS entries (they are identical to mRNA entries)
sed -i -e '/^#/d' -e '/^$/d' $GF1;
sed -i -Ee '/.+\t.+\tgene\t/d' -Ee '/.+\t.+\tCDS\t/d' -Ee 's/;Parent=gene_[0-9]+//g' $GF1;
## fix mRNA names
awk -F"\t" '{gsub("=tran_","="$1".")}1' $GF1 > tmp && mv tmp $GF1;
## remove empty lines and fix protein names to match mRNA names
sed -i -e '/^$/d' -Ee 's/^>gene_([0-9]+)\|.+>(.+)/>\2.\1/g' $FA1;

###### consolidate GeneMarkS predictions with Eukaryotic ones #####
## Only report those entries in GeneMarkS prediction (-a) that have no overlap in Gmove prediction (-b) (mRNA features only)
awk -F"\t" '{if($3 =="mRNA") print $0}' $gf | bedtools intersect -v -a $GF1 -b - > gms"$2"_"$n1".gff
## merge non-overlapping GeneMarkS predictions with Gmove predictions (this is the final gff)
cat $gf gms"$2"_"$n1".gff | sort -t $'\t' -V -k 1,1 -k 4,5 > tmp && mv tmp gms"$2"_"$n1".gff
## extract only GeneMarkS proteins that do not overlap with Gmove (and all Gmove proteins) (this is final protein fasta)
awk -F"\t|;" '{if($3 =="mRNA") print $9}' gms"$2"_"$n1".gff | sed 's/^ID=//g' | sort -Vu > tmp
cat $FA1 $fa | seqkit grep -n -f tmp - > gms"$2"_"$n1".faa;

### move final outputs to IVEX000 directory
mv gms"$2"_"$n1".gff ..
mv gms"$2"_"$n1".faa ..
rm -f tmp*
gzip -f $n2
