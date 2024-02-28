#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition long
#SBATCH --cpus-per-task 16
#SBATCH --mem 5GB

###### load software ######

###### input file variables ######

### Public genomes: 
#og="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/00__PUBLIC_GENOMES_ANNOTATIONS/*.gff"
#ofn="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/00__PUBLIC_GENOMES_ASSEMBLIES/*.fa"
#ofa="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/00__PUBLIC_GENOMES_PROTEOMES/*.fa"

### Phaeoexplorer genomes:
#og="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/02__FINAL_ANNOTATIONS/phex_1-20/*.gff" ## original gffs
#ofn="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/01__FINAL_GENOME_ASSEMBLIES/phex_1-20/*.fa" ## original fastas nucleotide
#ofa="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/04__FINAL_PROTEOMES/phex_1-20/*.fa" ## original fastas amino acid

og="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/02__FINAL_ANNOTATIONS/phex_21-33/*.gff" ## original gffs
ofn="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/01__FINAL_GENOME_ASSEMBLIES/phex_21-33/*.fa" ## original fastas nucleotide
ofa="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/04__FINAL_PROTEOMES/phex_21-33/*.fa" ## original fastas amino acid

##### UBE03_P007_nt/aa fields: (gff) $1-$8 -delimiter "\t", $9 -delimiter "\t;", $10-$15 -delimiter ";", (blast result - best viral hit) $16-$30 -delimiter ";", (blast result - best cellular hit) $31-$45 -delimiter ";", $46 -delimiter ";" best viral hit relative bitscore, $47 -delimiter ";" best cellular relative bitscore, $48 -delimiter ";" virus vs cellular category
##### now paste UBE03_P007_nt and UBE03_P007_aa together, and select either the best virus or cellular hit based on the virus vs cellular category
##### UBE03_P010 has the gff info with a single blast hit to each gene, either nt or aa. Whether nt or aa or viral or cellular hit was kept was based on the following priority list (the first being kept over the next and so on): VR_C[r,u,o] aa, V[r,u,o]_CR aa, Vr_C[r,u,o] aa, V[u,o]_Cr aa, Vu_C[u,o] aa, V[o]_Cu aa, THEN REPEATS SAME FOR NT, any combination of Vo_Co and NA/none

#srun paste -d ";" UBE03_P007_nt UBE03_P007_aa | cut -d ";" -f 40,80 | sed 's/$/;-/g' > UBE03_P008;

#srun awk -F ";" '{if($3 =="-" && $2 ~ "VR_") print $1";"$2";Vaa"; else print $0}' UBE03_P008 | awk -F ";" '{if($3 =="-" && $2 ~ "Vr_") print $1";"$2";Vaa"; else print $0}' | awk -F ";" '{if($3 =="-" && $2 ~ "Vu_Cu") print $1";"$2";Vaa"; else print $0}' | awk -F ";" '{if($3 =="-" && $2 ~ "Vu_Co") print $1";"$2";Vaa"; else print $0}' | awk -F ";" '{if($3 =="-" && $2 ~ "Vo_Cu") print $1";"$2";Caa"; else print $0}' | awk -F ";" '{if($3 =="-" && $2 ~ "_Cr") print $1";"$2";Caa"; else print $0}' | awk -F ";" '{if($3 =="-" && $1 ~ "Vr_") print $1";"$2";Vnt"; else print $0}' | awk -F ";" '{if($3 =="-" && $1 ~ "Vu_Cu") print $1";"$2";Vnt"; else print $0}' | awk -F ";" '{if($3 =="-" && $1 ~ "Vu_Co") print $1";"$2";Vnt"; else print $0}' | awk -F ";" '{if($3 =="-" && $1 ~ "Vo_Cu") print $1";"$2";Cnt"; else print $0}' | awk -F ";" '{if($3 =="-" && $1 ~ "_Cr") print $1";"$2";Cnt"; else print $0}' | awk -F ";" '{if($3 =="-" && $1 =="NA" && $2 =="NA") print $1";"$2";NA"; else print $0}' | awk -F ";" '{if($3 =="-" && $1 =="NA" && $2 ~ "Vo_Co") print $1";"$2";none"; else print $0}' | awk -F ";" '{if($3 =="-" && $1 ~ "Vo_Co" && $2 =="NA") print $1";"$2";none"; else print $0}' | awk -F ";" '{if($3 =="-" && $1 ~ "Vo_Co" && $2 ~ "Vo_Co") print $1";"$2";none"; else print $0}' > UBE03_P009;

#srun cut -d ";" -f 3 UBE03_P009 | paste -d ";" UBE03_P007_nt UBE03_P007_aa - | awk -F ";" '{print > "UBE03_P009_"$81}';
#srun cut -d ";" -f 1-7,8-22,40 UBE03_P009_Vnt > UBE03_P010; 
#srun cut -d ";" -f 1-7,23-37,40 UBE03_P009_Cnt >> UBE03_P010; cut -d ";" -f 41-47,48-62,80 UBE03_P009_Vaa >> UBE03_P010; 
#srun cut -d ";" -f 41-47,63-77,80 UBE03_P009_Caa >> UBE03_P010; 
#srun cut -d ";" -f 1-7,8-22,40 UBE03_P009_none | sed 's/;NA/;none/g' >> UBE03_P010; 
#srun cut -d ";" -f 1-7,8-22,40 UBE03_P009_NA >> UBE03_P010;
#srun sort -V -k 1,1 -k 4,4 -k 5,5 UBE03_P010 -o UBE03_P010;

#srun cut -d ";" -f 9 UBE03_P010 | sed 's/=/E_E_E/1' | sed -E 's/!|#|>|<|=|;|:|\[|\]|\|| |\/|\(|\)|,/_/g' | sed -E 's/___/_/g' | sed 's/__/_/g' | sed 's/E_E_E/=/g' > UBE03_P011;
#srun cut -d ";" -f 20 UBE03_P010 | sed 's/=/E_E_E/1' | sed -E 's/!|#|>|<|=|;|:|\[|\]|\|| |\/|\(|\)|,/_/g' | sed -E 's/___/_/g' | sed 's/__/_/g' | sed 's/E_E_E/=/g' > UBE03_P012;
#srun paste -d ";" UBE03_P010 UBE03_P011 UBE03_P012 | awk -F ";" 'BEGIN{OFS=";"} ; {print $1,$2,$3,$4,$5,$6,$7,$8,$24,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$25,$21,$22,$23}' > UBE03_D001; ## fix problematic characters in blast hit names

#### make UBE03_D002 which has both blast hits per genes and virus/cell category - "full data file", wheres D001 is the virus ID file
srun sed 's/=/E_E_E/1' UBE03_P007_nt | awk -F ";" 'BEGIN{OFS=";";a[9]=a[20]=a[24]=a[35]}; {for(x in a)gsub(/!|#|>|<|=|;|:|\[|\]|\|| |\/|\(|\)|,/,"_",$x);print}' | sed -E 's/___/_/g' | sed 's/__/_/g' | sed 's/E_E_E/=/g' | sed 's/_sseqid_/_sseqid=/g' | sed 's/_ssalltitles_/_ssalltitles=/g' > UBE03_Q001_nt;
srun sed 's/=/E_E_E/1' UBE03_P007_aa | awk -F ";" 'BEGIN{OFS=";";a[9]=a[20]=a[24]=a[35]}; {for(x in a)gsub(/!|#|>|<|=|;|:|\[|\]|\|| |\/|\(|\)|,/,"_",$x);print}' | sed -E 's/___/_/g' | sed 's/__/_/g' | sed 's/E_E_E/=/g' | cut --complement -d ";" -f 1-7 | sed 's/_sseqid_/_sseqid=/g' | sed 's/_ssalltitles_/_ssalltitles=/g' > UBE03_Q001_aa;
srun paste -d ";" UBE03_Q001_nt UBE03_Q001_aa > UBE03_D002;


