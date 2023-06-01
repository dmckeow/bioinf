#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 16
#SBATCH --mem 5GB

###### load software ######
module load diamond/0.9.36; module load seqtk/1.3;

###### input file variables ######
#ofa="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/04__FINAL_PROTEOMES/phex_1-20/*.fa" ## original fastas amino acid
#ofa="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/04__FINAL_PROTEOMES/phex_21-33/*.fa" ## original fastas amino acid
#ofa="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/00__PUBLIC_GENOMES_PROTEOMES/*.fa" ## original fastas amino acid

###### get list of all proteins with any viral hit
#srun awk -F ";" '{if($23 ~"VC_category=VR_" || $23 ~"VC_category=Vr_" || $23 ~"VC_category=Vu_") print $0}' UBE03_D001 | cut -d ";" -f 1-7 | sort -V -k 1,1 -k 4,4 -k 5,5 | awk -F"\t|;" '{print $9}' | sed 's/ID=//g' | sed -E 's/.$/& assembled CDS/g' > UBE04_A001;

###### get fasta files for diamond search
#for file in $ofa; do srun seqtk subseq $file UBE04_A001 | sed 's/\tassembled/ assembled/g' > UBE04_A001_ncvog_"$(basename $file)"; done;

###### rename some fastas (only needed for public genomes)
#for file in UBE04_A001_ncvog_UBE00_*_proteins.fa; do srun mv "$file" "${file//UBE04_A001_ncvog_UBE00_/UBE04_A001_ncvog_}"; done;

###### list input fastas for diamond search, then split into sublists of 5 genomes each
#srun find . -name "UBE04_A001_ncvog_*.fa" | split -l 5 - ; for file in x*; do srun mv "$file" "${file//x/UBE04_B001_x}"; done;

###### run this step once per UBE04_B001_x file (as $1 variable from terminal job submission), run all simultaneously if possible
for file in `(cat $1)`; do srun diamond blastp -d /projet/fr2424/sib/dmckeown/db/virus/ncvog.dmnd -q $file --more-sensitive -o $(basename "$file" .fa).blast -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles; done;
