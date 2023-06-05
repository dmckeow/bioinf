#!/bin/bash

umask 002
green='\033[32m';red='\033[31m';cyan='\033[36m';nocolor='\033[m'
echo -e "${green}RUNNING bioinf-setup with the following parameters: $@ ${nocolor}"

Help()
{
echo -e "${green}This script will setup the databases needed for the all of the Schroeder lab scripts and pipelines. Run this script from the directory in which it is found${nocolor}"
echo -e "${red}NOTE 1 - some of the databases generated are quite large - currently around 200 GB for the kaiju databases. If someone in your group has already ran this full script, then you could save space and time by simply running this script with the --preexisting option${nocolor}"
echo -e "${green}NOTE 2 - Similarly, if you want to update your system's paths to your own scripts and databases, without re-downloading or setting up anything, simply run Step A1 of this script on its own - e.g. -s A1 ${nocolor}"
echo -e "-t --tmp\tThe absolute path to a temporary large storage folder - e.g. /scratch.global - a folder named after your username will be made there"
echo -e "\n-s --step\t[optional] Which script steps to run. If blank, whole script is run. You can run multiple specific steps - e.g. A1_A2_A6 will only run steps A1, A2 and A6"
echo -e "\t${green}STEPS AVAILABLE: ${nocolor}"
awk '/^###### STEP-/ {print "\t\t"$0}' bioinf-setup.sh
echo -e "-P --preexisting\t[optional] The absolute path to a previously setup bioinf db directory e.g. /home/bioinf/db - with this option provided, the script will make soft links to all the database files found within the pre-existing bioinf setup, and link your setup to them. Doing this avoids having to re-download and rebuild databases that might already exist within your system or computing group"
echo -e "-d --dmnd\t[optional] A protein fasta file to generate an additional custom dmnd database. Use -s B1 to run this in isolation. Can be run as many times as you want, and all databases will used in the binning pipeline"
echo -e "-k --kaiju\t[optional] Specify one of kaiju's databses to build, in addition to the defaults made by this script. Do kaiju-makedb --help to see the options (with the conda environment bioinftools active) - e.g. -k fungi . Use -s B2 to run this in isolation. Can be run as many times as you want, and all databases will used in the binning pipeline"
echo -e "-T, --threads\tnumber of threads for job (SLURM --cpus-per-task does the same thing); default 1"
echo -e "-h, --help\tshow this help message and exit\n"
echo -e "\nEXAMPLES FOR RUNNING SCRIPT:\n\tSUBMIT TO SLURM (see README for more info):\nEXAMPLES FOR RUNNING SCRIPT:${cyan}sbatch --cpus-per-task=12 --time=96:00:00 --mem=240GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err bioinf-setup.sh -t /scratch.global${nocolor}"
}

while getopts t:s:P:d:k:T:h option
do 
    case "${option}" in
        s)step=${OPTARG};;
        t)tmp=${OPTARG};;
        P)preexisting=${OPTARG};;
        d)dmnd=${OPTARG};;
        k)kaiju=${OPTARG};;
        T)threads=${OPTARG};;
        h)Help; exit;;
    esac

done


####################### SET AND CHECK VARIABLES/ARGUMENTS #####################################
if [[ -z "${tmp}" ]]; then echo -e "${red}-t, --tmp REQUIRED. You must provide the script a location to place the large temporary files from various processes${nocolor}"; exit; fi

if [[ ! -z "${preexisting}" ]]; then step="A1_A2"; echo -e "${red}Running with --preexisting option: setup will create soft links to all files/folders within: ${preexisting}${nocolor}"; fi


### THREADS

if [[ ! -z "${SLURM_CPUS_PER_TASK}" ]]; then THREADS="${SLURM_CPUS_PER_TASK}"; echo -e "${green}Using threads set by SLURM_CPUS_PER_TASK:${nocolor}"; fi

if [[ ! -z "${threads}" ]]; then THREADS="${threads}"; echo -e "${green}Using threads set by --threads flag:${nocolor}"; fi

if [[ -z "${threads}" ]] && [[ -z "${SLURM_CPUS_PER_TASK}" ]]; then THREADS="1"; echo -e "${green}No SLURM_CPUS_PER_TASK or --threads set, using default threads:${nocolor}"; fi
echo -e "\t${THREADS} threads"

####################### CHECK FOR DEPENDENCIES #####################################

if [[ ! -f ./bioinf-setup.sh ]]; then echo "NOT in the bioinf directory where this script is located - cd there and try again"; exit; fi
if [[ -z $(conda env list | grep "bioinftools") ]]; then echo "NO conda environment for bioinftools found - see README"; else echo "conda environment for bioinftools FOUND"; fi

####################### DATABASES #####################################


####################### SOFTWARE #####################################
module purge
eval "$(conda shell.bash hook)"

conda activate bioinftools


####################################### SCRIPT #########################################################

###### STEP-A1: setup environmental variables and local folders
######################################################################
if [[ -z "${step}" ]] || [[ "$step" =~ "A1" ]]; then
######################################################################

sed -i '/.*BIOINF - SCHROEDER.*/d' $HOME/.bashrc ## remove any previous environment variables for bioinf schroeder

sed -i -z "s|$|\n###### ENV VARS FOR BIOINF - SCHROEDER GROUP UMN ######\nexport PATH=\"$PWD/bin:\$PATH\" ###### BIOINF - SCHROEDER\n|g" $HOME/.bashrc ## export all scripts in bin to the path

sed -i -z "s|$|export bioinfdb=$PWD/db ###### BIOINF - SCHROEDER\n|g" $HOME/.bashrc ## export path for database directory 
sed -i -z "s|$|export bioinftmp=$tmp/$USER ###### BIOINF - SCHROEDER\n|g" $HOME/.bashrc ## export path for tmp directory 

source $HOME/.bashrc

mkdir -p "$bioinftmp"
mkdir -p "$bioinfdb"/{CUSTOM_DMND,CUSTOM_HMMS,kaiju_nr_euk,kaiju_rvdb,VOGDB} ### make folders for the databases

######################################################################
fi
######################################################################

###### STEP-A2: create links to a pre-existing database setup for bioinf
######################################################################
if ([[ -z "${step}" ]] || [[ "$step" =~ "A2" ]]) && [[ ! -z "${preexisting}" ]]; then
######################################################################

### make the empty directories if the pre-exsiting bioinfdb has extra folders
cd "$bioinfdb"
for f in $preexisting/db/*; do
    mkdir -p $(basename $f)
done
### create symlinks for everything in the main directories
for f in $preexisting/db/*; do
    cd $(basename $f)
        for g in ${f}/*; do
          ln -s $g
        done
    cd "$bioinfdb"
done

exit
######################################################################
fi
######################################################################

###### STEP-A3: setup VOGDB
######################################################################
if [[ -z "${step}" ]] || [[ "$step" =~ "A3" ]]; then
######################################################################
T_O="tmp.A3"

cd "$bioinfdb"/VOGDB
rm -fr "$bioinfdb"/VOGDB/*
##### get latest VOGDB - a database of curated core genes for a wide range of virus groups
## also creates a file that will add more all useful info to VOGs
### get latest VOGDB
wget --cut-dirs 2 -r -np -nH -R "index.html*" https://fileshare.csb.univie.ac.at/vog/latest/

### make a file to add more useful info the VOG codes
zcat vog.annotations.tsv.gz | sed '1,1d' | sort -t $'\t' -k 1,1V > "$T_O".1
zcat vog.lca.tsv.gz | sed '1,1d' | sort -t $'\t' -k 1,1V > "$T_O".2
paste "$T_O".1 "$T_O".2 | awk -F "\t" '{print $1"\t"$5"\t"$9}' | awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/.*\|/,"",$2)}1' | awk -F"\t" 'BEGIN{OFS="\t"};{gsub(/.*;/,"",$3)}1' | sed -E 's/[^a-zA-Z0-9\t]/_/g' | sed -E -e 's/_+/_/g' -e 's/\t_|_\t/\t/g' | awk -F "\t" '{print "/_____"$1"_____/s/_____"$1"_____/"$1"__"$2"__"$3"/g"}' > VOGDB_add_info

rm -f *"$T_O"*
######################################################################
fi
######################################################################

###### STEP-A4: setup HMMs VOGDB for ANVI'O to use
######################################################################
if [[ -z "${step}" ]] || [[ "$step" =~ "A4" ]]; then
######################################################################
T_O="tmp.A4"

cd "$bioinfdb"/CUSTOM_HMMS
rm -fr VOGDB; mkdir VOGDB; cd VOGDB

rm -fr "$T_O"_VOGDB; mkdir "$T_O"_VOGDB

pigz -f -p $THREADS -dc "$bioinfdb"/VOGDB/vog.hmm.tar.gz | tar --directory ./"$T_O"_VOGDB -xf -

rm -f "$T_O".1; touch "$T_O".1
for f in "$T_O"_VOGDB/VOG*.hmm; do
    cat $f >> "$T_O".1
done

sed -i -E 's/^(NAME +)(VOG[0-9]+)$/\1_____\2_____/g' "$T_O".1
sed 's/^\/_____/\/^NAME  _____/g' "$bioinfdb"/VOGDB/VOGDB_add_info > "$T_O".2

rm -f *"$T_O".1.part.*

N=$(wc -l "$T_O".1 | awk -F " " -v S="$THREADS" '{printf "%0.0f\n" ,$1/S}')
split -d -l $N "$T_O".1 "$T_O".1.part.
for f in "$T_O".1.part.*; do
    sed -f "$T_O".2 $f > "$T_O".2.${f} &
done
wait

cat "$T_O".2."$T_O".1.part.* > genes.hmm

pigz -f -p $THREADS genes.hmm

zcat genes.hmm.gz | awk '/^NAME *VOG[0-9]+/' | sed -E 's/^NAME +//g' | sort -Vu | awk '{print $0"\tnone""\tVOGDB"}' | sed -z 's/^/gene\taccession\thmmsource\n/1' > genes.txt

echo -e "-E 1e-20" > noise_cutoff_terms.txt
echo -e "VOGDB" > kind.txt
echo "https://vogdb.org/" > reference.txt
echo -e "AA:GENE" > target.txt

rm -fr *"$T_O"*

######################################################################
fi
######################################################################

###### STEP-A5: setup VOGDB for DMND
######################################################################
if [[ -z "${step}" ]] || [[ "$step" =~ "A5" ]]; then
######################################################################
T_O="tmp.A5"

cd "$bioinfdb"/CUSTOM_DMND
rm -fr "$bioinfdb"/CUSTOM_DMND/*
rm -fr "$bioinfdb"/VOGDB/tmp_VOG_faa; mkdir "$bioinfdb"/VOGDB/tmp_VOG_faa;
pigz -f -p $THREADS -dc "$bioinfdb"/VOGDB/vog.faa.tar.gz | tar --directory "$bioinfdb"/VOGDB/tmp_VOG_faa -xf -

for f in "$bioinfdb"/VOGDB/tmp_VOG_faa/VOG*.faa; do
  V=$(basename $f | sed 's/\.faa//g');
  awk '{if($0 ~ ">") gsub(/ .*/,"",$0)}1' $f | sed "s/>/>$V./g" > ${f}.2
done

rm -f "$bioinfdb"/VOGDB/tmp_VOG_faa/tmp4; touch "$bioinfdb"/VOGDB/tmp_VOG_faa/tmp4
for f in "$bioinfdb"/VOGDB/tmp_VOG_faa/VOG*.faa.2; do
    cat $f >> "$bioinfdb"/VOGDB/tmp_VOG_faa/tmp4
done

sed -i -E 's/>(VOG[0-9]+)(\..+)/>_____\1_____\2/g' "$bioinfdb"/VOGDB/tmp_VOG_faa/tmp4

N=$(wc -l "$bioinfdb"/VOGDB/tmp_VOG_faa/tmp4 | awk -F " " -v S="$THREADS" '{printf "%0.0f\n" ,$1/S}')
split -d -l $N "$bioinfdb"/VOGDB/tmp_VOG_faa/tmp4 "$bioinfdb"/VOGDB/tmp_VOG_faa/tmp4.part.

for f in "$bioinfdb"/VOGDB/tmp_VOG_faa/tmp4.part.*; do
    sed -f "$bioinfdb"/VOGDB/VOGDB_add_info $f > ${f}.tmp5 &
done
wait

rm -f "$bioinfdb"/VOGDB/vog.modified.faa; touch "$bioinfdb"/VOGDB/vog.modified.faa
for f in "$bioinfdb"/VOGDB/tmp_VOG_faa/tmp4.part.*.tmp5; do
    cat $f >> "$bioinfdb"/VOGDB/vog.modified.faa
done

rm -f vog.dmnd
diamond makedb --threads $THREADS --in "$bioinfdb"/VOGDB/vog.modified.faa -d vog.dmnd

pigz -f -p $THREADS "$bioinfdb"/VOGDB/vog.modified.faa

rm -fr "$bioinfdb"/VOGDB/tmp_VOG_faa

######################################################################
fi
######################################################################

###### STEP-A6: build the databases for Kaiju - one for cells (Archaea,Bacteria,Eukaryota - nr_euk) and one for viruses (rvdb)
######################################################################
if [[ -z "${step}" ]] || [[ "$step" =~ "A6" ]]; then
######################################################################
T_O="tmp.A6"

cd "$bioinftmp"

### replace all the eukaryotic taxa to include all Eukaryote, plus Archaea and Bacteria in the nr_euk db built by kaiju
rm -f $(which kaiju)-*taxonlistEuk.tsv
echo -e "Archaea\nBacteria\nEukaryota" | taxonkit name2taxid - | awk '{print $2"\t"$1}' > $(which kaiju)-*taxonlistEuk.tsv

rm -fr "$bioinfdb"/kaiju_rvdb/* rvdb
kaiju-makedb -s rvdb -t $THREADS
cp rvdb/kaiju_db_rvdb.fmi "$bioinfdb"/kaiju_rvdb/kaiju_db_rvdb.fmi
cp names.dmp "$bioinfdb"/kaiju_rvdb/names.dmp
cp nodes.dmp "$bioinfdb"/kaiju_rvdb/nodes.dmp
rm -fr rvdb names.dmp nodes.dmp

rm -fr "$bioinfdb"/kaiju_nr_euk/* nr_euk
kaiju-makedb -s nr_euk -t $THREADS
cp nr_euk/kaiju_db_nr_euk.fmi "$bioinfdb"/kaiju_nr_euk/kaiju_db_nr_euk.fmi
cp names.dmp "$bioinfdb"/kaiju_nr_euk/names.dmp
cp nodes.dmp "$bioinfdb"/kaiju_nr_euk/nodes.dmp
rm -fr nr_euk names.dmp nodes.dmp

######################################################################
fi
######################################################################


###### STEP-B1: setup your own specific DMND database
######################################################################
if ([[ -z "${step}" ]] || [[ "$step" =~ "B1" ]]) && [[ ! -z "${dmnd}" ]]; then
######################################################################
cd "$bioinfdb"/CUSTOM_DMND

dbname=$(basename ${dmnd} | sed 's/\./_/g')
diamond makedb --threads $THREADS --in ${dmnd} -d "$bioinfdb"/VOGDB/${dname}.dmnd

######################################################################
fi
######################################################################


###### STEP-B2: setup your own specific kaiju database
######################################################################
if ([[ -z "${step}" ]] || [[ "$step" =~ "B2" ]]) && [[ ! -z "${kaiju}" ]]; then
######################################################################

cd "$bioinftmp"


rm -fr "$bioinfdb"/kaiju_${kaiju}/* ${kaiju}
kaiju-makedb -s ${kaiju} -t $THREADS
mv ${kaiju}/kaiju_db_*.fmi "$bioinfdb"/kaiju_${kaiju}
mv names.dmp "$bioinfdb"/kaiju_${kaiju}
mv nodes.dmp "$bioinfdb"/kaiju_${kaiju}
rm -fr ${kaiju} names.dmp nodes.dmp


######################################################################
fi
######################################################################
