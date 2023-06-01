#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 8
#SBATCH --mem 24GB

###### load software ######
module load bedtools/2.29.2; module load seqkit/0.14.0;

###### input file variables ######

### NCVOG defintion files
NCVOGFUN="/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/NCVOG_functional_categories"
NCVOGDEF="/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/NCVOG_NCVOGdef"

### NCBI taxonomy files
TA1="/projet/fr2424/sib/dmckeown/db/NCBI_taxdump/rankedlineage.dmp"
TA2="/projet/fr2424/sib/dmckeown/db/NCBI_taxdump/fullnamelineage.dmp"

###### input file variables ######
PN1="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/00__PUBLIC_GENOMES_ASSEMBLIES/*.fa";
PN2="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/01__FINAL_GENOME_ASSEMBLIES/phex_1-20/*.fa";
PN3="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/01__FINAL_GENOME_ASSEMBLIES/phex_21-33/*.fa";
N1="/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/UBE_public_1/UBE05_B005_*_SUMMARY";
N2="/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/UBE_phex_1-20/UBE05_B005_*_SUMMARY";
N3="/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/UBE_phex_21-33/UBE05_B005_*_SUMMARY";
Q1="/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/UBE_public_1/UBE05_B005_*.gff";
Q2="/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/UBE_phex_1-20/UBE05_B005_*.gff";
Q3="/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/UBE_phex_21-33/UBE05_B005_*.gff";


###### prepare NCVOG keys ######
### prep gi,NCVOG code, funct name key for figs etc
#srun awk -F '\t' 'BEGIN{OFS= "\t"} ; {if($0 ~ /NCVOG/ && $10 != "") print $1,$7,"ncvog_prot_name="$10}' $NCVOGDEF | sort -Vu -k 1,1 | sort -V -k 2,2 | sed 's/NCVOG/ncvog=/g' > UBE05_K001;

## UBE05_K002 = $1 gi, $2 ncvog code, $3 protein name, $4 functional group
#srun cut -f 2 UBE05_K001 | sed 's/ncvog=/NCVOG/g' | awk 'FNR==NR { array[$1]=$2; next } { for (i in array) gsub("\\<"i"\\>", array[i]) }1'  $NCVOGFUN - | sed 's/^/ncvfg=/g' | sed 's/NCVOG[0-9]*/NA/g' | paste UBE05_K001 - | sed 's/#REF!/uncharacterized_protein/g' | sed -E 's/!|#|>|<|=|;|:|\[|\]|\|| |\/|\(|\)|,/_/g' | sed -E 's/___/_/g' | sed 's/__/_/g' | sed -E 's/\tncvog_([0-9]+)/\tncvog=\1/g' | sed -E 's/\tncvfg_/\tncvfg=/g' | sed -E 's/\tncvog_prot_name_/\tncvog_prot_name=/g' > UBE05_K002;

##### prepare ncvog labels to genes ######
### get best hit per gene vs NCVOG db + fix problematic characters in blast hit names
#srun awk -F "\t" '!a[$1]++' UBE04_A001_ncvog*blast | sed -E 's/!|#|>|<|=|;|:|\[|\]|\|| |\/|\(|\)|,/_/g' | sed -E 's/_+/_/g' > UBE05_A001;
#srun cut -f 1,2 UBE05_A001 | sed -E 's/(\tgi_[0-9]+).*/\1/g' | sed 's/$/_/g' > UBE05_A002;

### make gi - ncvog info key
#srun awk '{print "gi_"$1"_\t"$2";"$3";"$4}' UBE05_K002 | awk 'FNR==NR { array[$1]=$2; next } { for (i in array) gsub("\\<"i"\\>", array[i]) }1' - UBE05_A002 | sed 's/^/UBE05_ID=/g' | sed 's/\t/_\t/1' | sed -E '/gi_[0-9]+_/d' > UBE05_A003;

### get all gene names with NCVOG hits and split into smaller parts for array
#srun awk -F "\t|;" '{if($3 =="mRNA") print "UBE05_"$9"_" ; else print "NA;NA;NA"}' UBE03_D001 > UBE05_A005;
#srun split -l 100000 UBE05_A005 UBE05_B001_;

##### add ncvog labels to gff/blast info- slow ~10 minutes per ~200,000 lines - run this (simultaneously ideally) step for each UBE05_B001 file at once in fast.q
### target is duplicated gene name at end of D001 which is replaced by NCVOG if matched

#srun awk 'FNR==NR { array[$1]=$2; next } { for (i in array) gsub(i, array[i]) }1'  UBE05_A003 $1 | sed -E 's/UBE05.*/ncvog=NA;ncvog_prot_name=NA;ncvfg=NA/g' > UBE05_B002_"$1";

##### merge ncvog labels withg rest of data and make fnial .gffs

#for file in UBE05_B002_UBE05_B001_*; do srun mv "$file" "${file//UBE05_B002_UBE05_B001_/UBE05_B002_}"; done;
#srun cat UBE05_B002_* | paste -d ";" UBE03_D001 - > UBE05_B003; ## merge the output of previous step
#srun cat UBE05_B002_* | paste -d ";" UBE03_D002 - > UBE05_B004; ## merge the output of previous step

#srun sort -V -k 1,1 -k 4,4 -k 5,5 UBE05_B004 | cut -f 1 | sed -E 's/^(.*)_.ontig[0-9]*/\1/g' | sed 's/_sdr_f//g' | paste - UBE05_B004 | awk -F "\t" '{print > "UBE05_B006_"$1}'; for file in  UBE05_B006_*; do srun cut --complement -f 1 $file | sed 's/;NA//g' | sed -E 's/(;exons=[0-9]+);ncvog=NA;ncvog_prot_name=NA;ncvfg=NA/\1/g' > "$file".gff; done; ## split by genome, both cellular virus best hit and relative reciprocal blast scores per gene

#srun sort -V -k 1,1 -k 4,4 -k 5,5 UBE05_B004 | awk -F ";" 'BEGIN{OFS=";"}; {print $53,$71,$68,$72}' - > UBE05_B004.1;

#srun sort -V -k 1,1 -k 4,4 -k 5,5 UBE05_B003 | cut -f 1 | sed -E 's/^(.*)_.ontig[0-9]*/\1/g' | sed 's/_sdr_f//g' | paste - UBE05_B003 | paste -d";" - UBE05_B004.1 | awk -F "\t" '{print > "UBE05_B005_"$1}'; for file in  UBE05_B005_*; do srun cut --complement -f 1 $file | sed 's/;NA//g' | sed -E 's/(;exons=[0-9]+);ncvog=NA;ncvog_prot_name=NA;ncvfg=NA/\1/g' > "$file".gff; done; ## split by genome, one blast annotation with virus priority

##### simplified gffs to assess in excel, with help from count files- make notes in count file copy to ID contigs/clusters to be visualised
#echo "contig;start;end;gene;exons;pident;length;evalue;bitscore;ssalltitles;hit_type;sskingdoms;VC_category;NCVOG;NCVOG_name;NCVFG;vaa_ssalltitles;vaa_relative_bitscore;caa_ssalltitles;caa_relative_bitscore" | sed 's/;/\t/g' > UBE05_S001_header;
#for file in UBE05_B006_*.gff; do awk '$3 =="mRNA"' $file | awk -F ";" 'BEGIN{OFS="\t"}; {print $53,$71,$68,$72}' | sed 's/..._ssalltitles=\|..._relative_bitscore=//g' > UBE05_S002_"$(basename $file .gff)"; done; for file in UBE05_S002_UBE05_B006_*; do srun mv "$file" "${file//UBE05_S002_UBE05_B006_/UBE05_S002_UBE05_B005_}"; done;
#for file in UBE05_B005_*.gff; do awk '$3 =="mRNA"' $file | cut --complement -f 2,3,6-8 | cut -d ";" -f 1,7,10-11,18-20,22-27 | sed 's/;/\t/g' | sed -E 's/[a-z]([a-z][a-z])_sskingdoms=/\1\t/g' | sed 's/ncvog=\|ncvog_prot_name=\|ncvfg=\|ID=\|exons=\|VC_category=\|[a-z][a-z][a-z]_pident=\|[a-z][a-z][a-z]_ssalltitles=\|[a-z][a-z][a-z]_length=\|[a-z][a-z][a-z]_bitscore=\|[a-z][a-z][a-z]_evalue=//g' | sed 's/;/\t/g' | paste - UBE05_S002_"$(basename $file .gff)" | cat UBE05_S001_header - > $(basename $file .gff)_SUMMARY; done;

###### run on all genomes - and in separate subfolder as C and D make many tmp files #####

###### D get contexts of ALL contigs ######
### n/N=no result (because monoexonic with no viral hit or non-ME), c/C=cellular hit, o/O=orfan, v/V=virus hit; lowercase=non-ME, uppercase=ME;
#mkdir D_UBE05;
#cd D_UBE05; for file in $N1 $N2 $N3; do sed '1d' $file | cut -f 1,5,12 | awk '{if($2 =="1" && $3 =="NA") print $1"\tN"; else print $0}' | awk '{if($2 =="1" && $3 =="Cellular") print $1"\tC"; else print $0}' | awk '{if($2 =="1" && $3 =="none") print $1"\tO"; else print $0}' | awk '{if($2 =="1" && $3 =="Viruses") print $1"\tV"; else print $0}' | sed -E -e 's/\t[0-9]+\tNA/\tn/g' -e 's/\t[0-9]+\tCellular/\tc/g' -e 's/\t[0-9]+\tnone/\to/g' -e 's/\t[0-9]+\tViruses/\tv/g' | uniq -c | sed -E 's/ +([0-9]+) (.+)/\2\1-/g' | awk '{print $2 > "UBE05_D001_"$1}'; done;
#for file in UBE05_D001_*; do awk '{if(NR==1) print FILENAME"\t"$0; else print $0}' $file | tr -d '\n' | sed 's/-$/\n/g' | sed 's/UBE05_D001_//g'; done > UBE05_D002;
#for f in UBE05_D001_*; do rm -f $f; done;

#for file in $N1 $N2 $N3; do sed '1d' $file | cut -f 1,10 | awk '{print $2 > "UBE05_D003_"$1}'; done;
#for file in UBE05_D003_*; do awk '{if(NR==1) print FILENAME"\t"$0; else print $0}' $file | sed 's/$/:/g' | tr -d '\n' | sed 's/:$/\n/g' | sed 's/UBE05_D003_//g' | sed -E 's/(:)EsV-1-[0-9]+_Ectocarpus_siliculosus_virus_1_|_Ectocarpus_siliculosus_virus_1_|_Ectocarpus_siliculosus_/\1/g'; done > UBE05_D004;
#for f in UBE05_D003_*; do rm -f $f; done;
#cd ..;

##### N core gene count summary ######
### NCVOG count per contig
#mkdir N_UBE05; cd N_UBE05;
#srun awk -F "\t" '{if($14 ~ /[0-9]+/) print $1}' $N1 $N2 $N3 | sed '/^contig$/d' | sort -Vu > UBE05_N001; ## list contigs
#srun sed 's/ \|;/_/g' /projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/NCVOG_NCVOGdef | awk -F "\t" '{if($9 ~ /NCVOG[0-9]+/) print $9}'  | sed 's/NCVOG//g' | #sort -Vu > UBE05_N002; ## list ALL NCVOGs
#srun awk -F "\t" '{if($14 ~ /[0-9]+/) print $1"\t"$14}' $N1 $N2 $N3 | sort -V | awk -F "\t" '{print $2 > "UBE05_N003_"$1}'; 
#for f in UBE05_N003_*; do cat UBE05_N002 $f | sort -V | uniq -c | sed -E 's/ +([0-9]+) /\1\t/g' | awk -F "\t" -v f=$f '{print f"\t"$1-1"\t"$2}'; done | sed 's/UBE05_N003_//g' > UBE05_N004;
#awk '{print $2 > "UBE05_N005_"$3}' UBE05_N004;
#srun sed 's/ \|;/_/g' /projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/NCVOG_NCVOGdef | awk -F "\t" '{if($9 ~ /NCVOG[0-9]+/) print $9";"$10}'  | sed 's/NCVOG//g' | sort -Vu | sort -t ';' -V -k 1,1 | sed 's/$/\t/g' | tr -d '\n' | sed 's/\t$/\n/g' | sed 's/^/contig\tgenomic_context;n=no_result-because_monoexonic_with_no_viral_hit_or_non-ME-;c=cellular hit;o=orfan;v=virus hit;lower_or_uppercase=non-ME_or_ME\t/g' > UBE05_N006;

### get reference virus genome NCVOG presence/absence
#awk -F "\t" '{if($7 ~ /[0-9]+/ && $2 ~ "Ectocarpus") print $7}' /projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/NCVOG_NCVOGdef | sed 's/NCVOG//g' | cat - UBE05_N002 | sort -V | uniq -c | sed -E 's/ +([0-9]+) /\1\t/g' | awk -F "\t" '{print $1-1}' | sed 's/$/\t/g' | tr -d '\n' | sed -e 's/\t$/\n/g' -e 's/^/EsV-1\tNA\t/g' > UBE05_N007.1;
#awk -F "\t" '{if($7 ~ /[0-9]+/ && $2 ~ "Feldmannia") print $7}' /projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/NCVOG_NCVOGdef | sed 's/NCVOG//g' | cat - UBE05_N002 | sort -V | uniq -c | sed -E 's/ +([0-9]+) /\1\t/g' | awk -F "\t" '{print $1-1}' | sed 's/$/\t/g' | tr -d '\n' | sed -e 's/\t$/\n/g' -e 's/^/FsV-158\tNA\t/g' > UBE05_N007.2;

#paste UBE05_N005_* > UBE05_N009; ## merge all NCVOGs count (columns) and all contings (rows)
#cut -f 1 UBE05_N004 | awk '!a[$0]++' | sed 's/$/\t/g' | grep -f - ../D_UBE05/UBE05_D002 > UBE05_N010; ## get genome contexts
#paste UBE05_N010 UBE05_N009 | cat UBE05_N006 UBE05_N007.1 UBE05_N007.2 - > UBE05_N011;
### split UBE05_N011 to view in spreadsheets
#cut -f 1,2,3-800 UBE05_N011 > UBE05_N011.1; cut -f 1,2,801-1600 UBE05_N011 > UBE05_N011.2; cut -f 1,2,1601-2400 UBE05_N011 > UBE05_N011.3; cut -f 1,2,2401-3104 UBE05_N011 > UBE05_N011.4;

#echo ";NCVOG_group;1;1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;3;3;3;3;3;3;3;3;3;3;3;3;3;3;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;4;5;5;5" | sed 's/;/\t/g' > UBE05_N012;
#srun awk -F"\t" 'BEGIN{OFS="\t"}; {print $1,$2,$24,$25,$33,$40,$53,$77,$234,$280,$1106,$199,$222,$226,$255,$259,$295,$1012,$1288,$11,$12,$34,$36,$39,$42,$247,$254,$257,$261,$292,$1059,$1064,$1303,$14,$37,$63,$211,$215,$230,$231,$250,$256,$294,$911,$993,$1003,$1029,$1062,$1078,$1133,$1134,$1261,$1289,$1295,$78,$1019,$1117}' UBE05_N011 | cat UBE05_N012 - > UBE05_N013; ## reduce to group 1-5 NCVOGs only

##### get core gene counts per genome
#srun sed '1,3d' UBE05_N011 | cut --complement -f 2 | sed -E 's/_.ontig[0-9]+\t/;/g' | awk -F ";" '{print $2 > "UBE05_N014_"$1}';
#head -3 UBE05_N011 | cut --complement -f 2 > UBE05_N015.1;
#for f in UBE05_N014_*; do srun awk 'BEGIN{OFS="\t"}; {for(i=1;i<=NF;i++)$i=(a[i]+=$i)}END{print FILENAME,$0}' $f | sed 's/UBE05_N014_//g'; done | cat UBE05_N015.1 - > UBE05_N015; ## all NCVOGs

#srun sed 's/^/\t/g' UBE05_N015 | awk -F"\t" 'BEGIN{OFS="\t"}; {print $1,$2,$24,$25,$33,$40,$53,$77,$234,$280,$1106,$199,$222,$226,$255,$259,$295,$1012,$1288,$11,$12,$34,$36,$39,$42,$247,$254,$257,$261,$292,$1059,$1064,$1303,$14,$37,$63,$211,$215,$230,$231,$250,$256,$294,$911,$993,$1003,$1029,$1062,$1078,$1133,$1134,$1261,$1289,$1295,$78,$1019,$1117}' | cat UBE05_N012 - | cut --complement -f 1 > UBE05_N016; ## reduce to group 1-5 NCVOGs only

#cd ..;

###### get viral clusters #####
#mkdir C_UBE05;
#cd C_UBE05; for file in $N1 $N2 $N3; do srun sed '1d' $file | awk -F "\t" 'BEGIN{OFS="\t"}; {if($12 =="Viruses") print $0}' | awk -F "\t" 'BEGIN{OFS="\t"}; {if(f != $1) print $0"\t__SPLIT__"; else print $0} {f=$1}' | awk -F "\t" 'BEGIN{OFS="\t"}; {if($0 !~ "__SPLIT__" && $2 - i > 20000) print $0"\t__SPLIT__"; else print $0} {i=$3}' | csplit -q --prefix=UBE05_C001_"$(basename $file)"_ - '/__SPLIT__/' '{*}'; done;
#for file in UBE05_C001_UBE05_B005_*_SUMMARY_*; do srun cut -f 1,12 $file | uniq -c | sed -E 's/ +([0-9]+) (.+)/\1\t\2/g' > UBE05_C002_"$(basename $file)"; done;
#for file in UBE05_C002_UBE05_C001_*_SUMMARY_*; do srun awk -F"\t" '{if($1 > 2) print FILENAME}' $file | sort -Vu | sed 's/UBE05_C002_//g'; done > UBE05_C003;
### nucleotide
#for file in $(cat UBE05_C003); do head -1 $file | cut -f 1,2; done > UBE05_C004;
#for file in $(cat UBE05_C003); do tail -1 $file | cut -f 3; done | paste UBE05_C004 - | sed -E 's/(.+.ontig[0-9]+)\t([0-9]+)\t([0-9]+)/\1\t\2\t\3\t\1_\2-\3/g' > UBE05_C005;
#for file in $PN1 $PN2 $PN3; do bedtools getfasta -fi $file -bed UBE05_C005 > UBE05_C005_"$(basename $file)"; done;
#find UBE05_C005_*.fa -size 0 -delete;
#for file in UBE05_C005_*.fa; do sed -E 's/_.ontig([0-9]+:[0-9]+)-([0-9]+)/\t\1\t\2/g' $file | sed -E 's/_|-//g' | sed -E 's/(>[Aa-Zz][Aa-Zz])[Aa-Zz]+([Aa-Zz][Aa-Zz])/\1\2/g' | sed 's/\t/_/g' > MAUVE_"$(basename $file)"; done; cd ..;

##### get cluster size category counts, per genome
#for f in $Q1 $Q2 $Q3; do head -1 $f | cut -f 1 | sed -E 's/_.ontig[0-9]+//g' ;done | sort -V > UBE05_C000; ## list all genomes
#for file in UBE05_C002_UBE05_C001_*_SUMMARY_*; do awk -F"\t" '{if($1 >= 2 && $1 <= 5) print $2}' $file | sed -E 's/(.+)_.ontig.+/\1/g'; done | cat UBE05_C000 - | sort -V | uniq -c | sed -E 's/ +([0-9]+) (.+)/\1\t\2/g' | sort -V -k 2,2 | awk '{print $1-1}' > UBE05_C003.1; ## count clusters 2-5 genes
#for file in UBE05_C002_UBE05_C001_*_SUMMARY_*; do awk -F"\t" '{if($1 >= 6 && $1 <= 20) print $2}' $file | sed -E 's/(.+)_.ontig.+/\1/g'; done | cat UBE05_C000 - | sort -V | uniq -c | sed -E 's/ +([0-9]+) (.+)/\1\t\2/g' | sort -V -k 2,2 | awk '{print $1-1}' > UBE05_C003.2; ## count clusters 6-20 genes
#for file in UBE05_C002_UBE05_C001_*_SUMMARY_*; do awk -F"\t" '{if($1 >= 21 && $1 <= 50) print $2}' $file | sed -E 's/(.+)_.ontig.+/\1/g'; done | cat UBE05_C000 - | sort -V | uniq -c | sed -E 's/ +([0-9]+) (.+)/\1\t\2/g' | sort -V -k 2,2 | awk '{print $1-1}' > UBE05_C003.3; ## count clusters 21-50 genes
#for file in UBE05_C002_UBE05_C001_*_SUMMARY_*; do awk -F"\t" '{if($1 >= 51 && $1 <= 100) print $2}' $file | sed -E 's/(.+)_.ontig.+/\1/g'; done | cat UBE05_C000 - | sort -V | uniq -c | sed -E 's/ +([0-9]+) (.+)/\1\t\2/g' | sort -V -k 2,2 | awk '{print $1-1}' > UBE05_C003.4; ## count clusters 51-100 genes
#for file in UBE05_C002_UBE05_C001_*_SUMMARY_*; do awk -F"\t" '{if($1 >= 101 && $1 <= 200) print $2}' $file | sed -E 's/(.+)_.ontig.+/\1/g'; done | cat UBE05_C000 - | sort -V | uniq -c | sed -E 's/ +([0-9]+) (.+)/\1\t\2/g' | sort -V -k 2,2 | awk '{print $1-1}' > UBE05_C003.5; ## count clusters 101-200 genes
#for file in UBE05_C002_UBE05_C001_*_SUMMARY_*; do awk -F"\t" '{if($1 >= 201 && $1 <= 300) print $2}' $file | sed -E 's/(.+)_.ontig.+/\1/g'; done | cat UBE05_C000 - | sort -V | uniq -c | sed -E 's/ +([0-9]+) (.+)/\1\t\2/g' | sort -V -k 2,2 | awk '{print $1-1}' > UBE05_C003.6; ## count clusters 201-300 genes
#for file in UBE05_C002_UBE05_C001_*_SUMMARY_*; do awk -F"\t" '{if($1 >= 301) print $2}' $file | sed -E 's/(.+)_.ontig.+/\1/g'; done | cat UBE05_C000 - | sort -V | uniq -c | sed -E 's/ +([0-9]+) (.+)/\1\t\2/g' | sort -V -k 2,2 | awk '{print $1-1}' > UBE05_C003.7; ## count clusters 301+ genes
#echo "## counts of genes in virus gene clusters of range x to y, separated by <20,000 bp between virus genes:genome;2-5 genes;6-20 genes;21-50 genes;51-100 genes;101-200 genes;201-300 genes;301+ genes" | sed -e 's/:/\n/g' -e 's/;/\t/g' > UBE05_C003.8;
#paste UBE05_C000 UBE05_C003.1 UBE05_C003.2 UBE05_C003.3 UBE05_C003.4 UBE05_C003.5 UBE05_C003.6 UBE05_C003.7 | cat UBE05_C003.8 - > UBE05_C003.9; ## FINAL FILE, counts per genome

###### get general virus count info #####
#mkdir R_UBE05;
#cd R_UBE05;
### for all virus amino acid hits, get: $1 genome, $2 contig, $3 vaa_staxids, $4 ncvog, $5 vaa_relative_bitscore=, $6 caa_relative_bitscore=
#for f in $Q1 $Q2 $Q3; do head -1 $f | cut -f 1 | sed -E 's/_.ontig[0-9]+//g' ;done | sort -V > UBE05_R000; ## list all genomes
#for f in $Q1 $Q2 $Q3; do awk '{if($3 =="mRNA") print FILENAME"\t"$0}' $f | sed -E 's/.+UBE05_B005_(.+).gff/\1/g' | awk -F"\t|;" 'BEGIN{OFS="\t"}; {if($0 ~ "vaa_qseqid=") print $1,$2,$30,$33,$37,$39}' | sed -E 's/vaa_staxids=|ncvog=|.aa_relative_bitscore=//g' | sed -E 's/\tnone|\tUR/\t0.00/g' | sed -E 's/(.+\t.+\t[0-9]+)_.+(\t.+\t.+\t.+)/\1\2/g' ; done > UBE05_R001;
### round relative bitscores to 1 decimal place
#awk '{printf "%.1f\n", $5}' UBE05_R001 | sed -E 's/1\.[0-9]+/1.0/g' > UBE05_R001.1; awk '{printf "%.1f\n", $6}' UBE05_R001 | sed -E 's/1\.[0-9]+/1.0/g' | paste UBE05_R001 UBE05_R001.1 - | cut --complement -f 5,6 | awk -F "\t" '{print $0"\t"$5 - $6}' | sort -n -k 7,7 > UBE05_R002;
#awk -F "\t" '{print $1 > "UBE05_R002.1_"$7}' UBE05_R002;

##### per genome:
##### count all reciprocated virus aa hits, which have: virus vs cellular (absent cell hit = 0.0 relative bitscore) and virus taxa/or NCVOG
#for f in `(cat UBE05_R000)`; do grep -c $f UBE05_R002.1_* | sed 's/UBE05_R002.1_//g' | sed 's/:/\t/g' | sort -n -k 1,1 > UBE05_R002.2_"$f" ; done;
#for f in UBE05_R002.2_*; do cut -f 2 $f | sed 's/$/\t/g' | tr -d '\n' | sed 's/\t$/\n/g' > UBE05_R002.3_"$f"; done;
#echo "## These are all virus hits counted, per brown algal genome, from amino acid BLAST. Numerical categories (columns) are the virus relative BLAST bitscore (rbitscore) minus the cellular rbitscore (set to 0 if cellular hit absent). This indicates whether a gene has higher similarity to viruses or cells. Lower scores (down to the minimum of -1) are more cellular and higher scores are more viral (up to the max of 1). For example, a score of zero means the virual and cellular rbitscores were identical - this is an ambiguous virus/hit. A score of 1 is a gene with a virus hit and no cellular hit, a clear viral gene:## ;cell only hit;;;;;<- increasing cell hit;;;;;virus or host;;;;;increasing virus hit ->;;;;;virus only hit:genome;-1;-0.9;-0.8;-0.7;-0.6;-0.5;-0.4;-0.3;-0.2;-0.1;0;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;1" | sed -e 's/:/\n/g' -e 's/;/\t/g' > UBE05_R002.4;
#cat UBE05_R002.3_* | paste UBE05_R000 - | cat UBE05_R002.4 - > UBE05_R002.5;

##### HEATMAP: count all virus taxa (lowest taxa above species), per genome, aa hits, which have: virus > cellular rbitscore by 0.3
### remove known false ambiguous hits
#for f in $Q1 $Q2 $Q3; do awk '{if($3 =="mRNA") print FILENAME"\t"$0}' $f | sed -E 's/.+UBE05_B005_(.+).gff/\1/g' | awk -F"\t|;" 'BEGIN{OFS="\t"}; {if($0 ~ "vaa_qseqid=") print $1,$29,$30,$33,$37,$39}' | sed -E 's/vaa_staxids=|ncvog=|.aa_relative_bitscore=|..._ssalltitles=//g' | sed -E 's/\tnone|\tUR/\t0.00/g' | sed -E 's/(.+\t.+\t[0-9]+)_.+(\t.+\t.+\t.+)/\1\2/g'; done > UBE05_R001_taxa;
### round relative bitscores to 1 decimal place
#awk '{printf "%.1f\n", $5}' UBE05_R001_taxa | sed -E 's/1\.[0-9]+/1.0/g' > UBE05_R001.1_taxa; awk '{printf "%.1f\n", $6}' UBE05_R001_taxa | sed -E 's/1\.[0-9]+/1.0/g' | paste UBE05_R001_taxa UBE05_R001.1_taxa - | cut --complement -f 5,6 | awk -F "\t" '{print $0"\t"$5 - $6}' | sort -n -k 7,7 > UBE05_R002_taxa;
#sed -E -e 's/\t//g' -e 's/\|/\t/g' $TA2 | awk -F"\t" '{print $2"\t~"$1"~"}' | sed -E 's/!|#|>|<|=|;|:|\[|\]|\|| |\/|\(|\)|,/_/g' | sed -E 's/___/_/g' | sed 's/__/_/g' > UBE05_R002.1_taxa; ## key of virus name in NCBI and taxid
#awk '{if($3 =="NA") print $2}' UBE05_R002_taxa | sed 's/_$//g'> UBE05_R002.2_taxa; ## list result hits with no taxid
#cut -f 1 UBE05_R002.1_taxa | grep -oF -f - UBE05_R002.2_taxa | sort -Vu > UBE05_R002.3_taxa; ## reduce to relevant taxa only
#grep -wF -f UBE05_R002.3_taxa UBE05_R002.1_taxa > UBE05_R002.4_taxa; ## reduce to relevant taxa only
#awk -F "\t" 'FNR==NR { array[$1]=$2; next } { for (i in array) gsub(i, array[i]) }1' UBE05_R002.4_taxa UBE05_R002.2_taxa | sed -E 's/.+\~([0-9]+)\~.*/\1/g' > UBE05_R002.5_taxa; ## replace subject names with taxids
#awk '{if($3 =="NA") print $0}' UBE05_R002_taxa | paste - UBE05_R002.5_taxa | awk -F "\t" 'BEGIN{OFS="\t"}; {print $1,$2,$8,$4,$5,$6,$7}' | awk -F "\t" 'BEGIN{OFS="\t"}; {if($3 ~ /[a-z]+|[A-Z]+/) print $1,$2,"NA",$4,$5,$6,$7; else print $0}' > UBE05_R002.6_taxa
#awk '{if($3 !="NA") print $0}' UBE05_R002_taxa | cat - UBE05_R002.6_taxa > UBE05_R002.7_taxa; ## all missing taxids (RVDB) restored
#for f in `(cat UBE05_R000)`; do awk 'BEGIN{OFS="\t"}; {if($7 >= 0.3) print $1,$3}' UBE05_R002.7_taxa | sed '/\tNA/d' | sed -E 's/(.+\t)([0-9]+)/\1__\2__/g' > UBE05_R003; done; ## taxa, virus > cell, any quality
#cut -f 2 UBE05_R003 | sort -Vu > UBE05_R003.1;
#sed 's/\t|\t/;/g' $TA1 | sed 's/;/__\t/1' | sed 's/\t|//g' | sed 's/^/__/1' | grep -f UBE05_R003.1 - | awk -F "\t" 'FNR==NR { array[$1]=$2; next } { for (i in array) gsub("\\<"i"\\>", array[i]) }1' - UBE05_R003 > UBE05_R003.2; ## replace taxids with taxonomy, taxa count, virus > cell, any quality
#awk -F"\t" '{if($2 ~ /virus|viridae|virales|viricota/) print $2}' UBE05_R003.2 | awk -F";" '{print $4";"$3}' | sed -E -e 's/;+$|^;+//g' | sort -Vu > UBE05_R003.3;
#awk -F"\t" '{if($2 ~ /virus|viridae|virales|viricota/) print $2 > "UBE05_R003.4_"$1}' UBE05_R003.2;
#for f in UBE05_R003.4_*; do awk -F";" '{print $4";"$3}' $f | sed -E -e 's/;+$|^;+//g' | cat UBE05_R003.3 - | sort -V | uniq -c | sed -E 's/ +([0-9]+) /\1\t/g' | sort -V -k 2,2 | awk -F "\t" '{print $1-1}' > UBE05_R003.5_"$f"; done ## count taxa
#for f in UBE05_R003.5_*; do awk 'FNR==NR{s+=$1;f=FILENAME;next}; {print f"\t"$1/s*100}' $f $f | sed 's/UBE05_R003.5_UBE05_R003.4_//g' | paste - UBE05_R003.3; done | awk -F"\t" '{if($3 !="") print $0}' | sed -ze 's/^/genome\tpercentage\ttaxa\n/g' -e 's/\t0\t/\t\t/g' > UBE05_R003_taxa_heatmap; ## convert to % of v genes per genome


###### make file to plot virus vs cell rbitscores blastp ggenes, all genomes
### each UBE05_R004.3_UBE05_R000_[xx] makes a wrapped figure of 6 scattergraphs
#awk -F "\t" 'BEGIN{OFS="\t"}; {if($5 != 0.0 && $3 !~ "NA") print $1,$5,$6,$4}' UBE05_R002.7_taxa | sort -V -k 1,1 | sed -z 's/^/genome\tviral_blastp_relative_bitscore\tcellular_blastp_relative_bitscore\tNCVOG\n/' | sed 's/\tNA/\tnone/g' > UBE05_R004.1;
#echo "0022;Group_1;;0023;Group_1;;0031;Group_1;;0038;Group_1;;0052;Group_1;;0076;Group_1;;0249;Group_1;;0304;Group_1;;1164;Group_1;;0211;Group_2;;0237;Group_2;;0241;Group_2;;0272;Group_2;;0276;Group_2;;0320;Group_2;;1068;Group_2;;1353;Group_2;;0009;Group_3;;0010;Group_3;;0032;Group_3;;0034;Group_3;;0037;Group_3;;0040;Group_3;;0262;Group_3;;0271;Group_3;;0274;Group_3;;0278;Group_3;;0317;Group_3;;1117;Group_3;;1122;Group_3;;1368;Group_3;;0012;Group_4;;0035;Group_4;;0062;Group_4;;0225;Group_4;;0230;Group_4;;0245;Group_4;;0246;Group_4;;0267;Group_4;;0273;Group_4;;0319;Group_4;;0965;Group_4;;1049;Group_4;;1059;Group_4;;1087;Group_4;;1120;Group_4;;1136;Group_4;;1191;Group_4;;1192;Group_4;;1326;Group_4;;1354;Group_4;;1360;Group_4;;0077;Group_5;;1076;Group_5;;1175;Group_5:" | sed -e 's/;;/\n/g' -e 's/;/\t/g' -e 's/://g' > UBE05_R004.2;
#awk 'FNR==NR { array[$1]=$2; next } { for (i in array) gsub(i, array[i]) }1' UBE05_R004.2 UBE05_R004.1 | sed -Ee 's/[0-9][0-9][0-9][0-9]/none/g' | sed '1d' | sed -E 's/Group_2|Group_3|Group_4|Group_5/Group_2+/g' | sort -V -k 1,4 | uniq -c | sed -E 's/ +([0-9]+) (.+)/\2\t\1/g' > UBE05_R004;
#split -l 6 UBE05_R000 UBE05_R000_; ## list genomes in groups of 6
#for f in UBE05_R000_*; do grep -f $f UBE05_R004 | sort -V -k 1,1 -k 4,4r | sed -z 's/^/genome\tviral_blastp_relative_bitscore\tcellular_blastp_relative_bitscore\tNCVOG\tcount\n/' > UBE05_R005_"$f"; done

###### make file to plot virus vs cell rbitscores blastp ggenes, by contig
### each UBE05_R004.3_UBE05_R000_[xx] makes a wrapped figure of 6 scattergraphs

#for f in $Q1 $Q2 $Q3; do awk '{if($3 =="mRNA") print FILENAME"\t"$0}' $f | sed -E 's/.+UBE05_B005_(.+).gff/\1/g' | awk -F"\t|;" 'BEGIN{OFS="\t"}; {if($0 ~ "vaa_qseqid=") print $1"-----"$2,$29,$30,$33,$37,$39}' | sed -E 's/vaa_staxids=|ncvog=|.aa_relative_bitscore=|..._ssalltitles=//g' | sed -E 's/\tnone|\tUR/\t0.00/g' | sed -E 's/(.+\t.+\t[0-9]+)_.+(\t.+\t.+\t.+)/\1\2/g'; done > UBE05_R001_contig;
### round relative bitscores to 1 decimal place
#awk '{printf "%.1f\n", $5}' UBE05_R001_contig | sed -E 's/1\.[0-9]+/1.0/g' > UBE05_R001.1_contig; awk '{printf "%.1f\n", $6}' UBE05_R001_contig | sed -E 's/1\.[0-9]+/1.0/g' | paste UBE05_R001_contig UBE05_R001.1_contig - | cut --complement -f 5,6 | awk -F "\t" '{print $0"\t"$5 - $6}' | sort -n -k 7,7 > UBE05_R002_contig;
#sed -E -e 's/\t//g' -e 's/\|/\t/g' $TA2 | awk -F"\t" '{print $2"\t~"$1"~"}' | sed -E 's/!|#|>|<|=|;|:|\[|\]|\|| |\/|\(|\)|,/_/g' | sed -E 's/___/_/g' | sed 's/__/_/g' > UBE05_R002.1_contig; ## key of virus name in NCBI and taxid
#awk '{if($3 =="NA") print $2}' UBE05_R002_contig | sed 's/_$//g'> UBE05_R002.2_contig; ## list result hits with no taxid
#cut -f 1 UBE05_R002.1_contig | grep -oF -f - UBE05_R002.2_contig | sort -Vu > UBE05_R002.3_contig; ## reduce to relevant taxa only
#grep -wF -f UBE05_R002.3_contig UBE05_R002.1_contig > UBE05_R002.4_contig; ## reduce to relevant taxa only
#awk -F "\t" 'FNR==NR { array[$1]=$2; next } { for (i in array) gsub(i, array[i]) }1' UBE05_R002.4_contig UBE05_R002.2_contig | sed -E 's/.+\~([0-9]+)\~.*/\1/g' > UBE05_R002.5_contig; ## replace subject names with taxids
#awk '{if($3 =="NA") print $0}' UBE05_R002_contig | paste - UBE05_R002.5_contig | awk -F "\t" 'BEGIN{OFS="\t"}; {print $1,$2,$8,$4,$5,$6,$7}' | awk -F "\t" 'BEGIN{OFS="\t"}; {if($3 ~ /[a-z]+|[A-Z]+/) print $1,$2,"NA",$4,$5,$6,$7; else print $0}' > UBE05_R002.6_contig
#awk '{if($3 !="NA") print $0}' UBE05_R002_contig | cat - UBE05_R002.6_contig > UBE05_R002.7_contig; ## all missing taxids (RVDB) restored


#awk -F "\t" 'BEGIN{OFS="\t"}; {if($5 != 0.0 && $3 !~ "NA") print $1,$5,$6,"~"$4"~"}' UBE05_R002.7_contig | sort -V -k 1,1 | sed -z 's/^/contig\tviral_blastp_relative_bitscore\tcellular_blastp_relative_bitscore\tNCVOG\n/' | sed 's/\tNA/\tnone/g' > UBE05_R004.1_contig;
#echo "0022;Group_1;;0023;Group_1;;0031;Group_1;;0038;Group_1;;0052;Group_1;;0076;Group_1;;0249;Group_1;;0304;Group_1;;1164;Group_1;;0211;Group_2;;0237;Group_2;;0241;Group_2;;0272;Group_2;;0276;Group_2;;0320;Group_2;;1068;Group_2;;1353;Group_2;;0009;Group_3;;0010;Group_3;;0032;Group_3;;0034;Group_3;;0037;Group_3;;0040;Group_3;;0262;Group_3;;0271;Group_3;;0274;Group_3;;0278;Group_3;;0317;Group_3;;1117;Group_3;;1122;Group_3;;1368;Group_3;;0012;Group_4;;0035;Group_4;;0062;Group_4;;0225;Group_4;;0230;Group_4;;0245;Group_4;;0246;Group_4;;0267;Group_4;;0273;Group_4;;0319;Group_4;;0965;Group_4;;1049;Group_4;;1059;Group_4;;1087;Group_4;;1120;Group_4;;1136;Group_4;;1191;Group_4;;1192;Group_4;;1326;Group_4;;1354;Group_4;;1360;Group_4;;0077;Group_5;;1076;Group_5;;1175;Group_5:" | sed -e 's/;;/\n/g' -e 's/;/\t/g' -e 's/://g' | sed -E 's/^([0-9]+)\t/~\1~\t/g' > UBE05_R004.2_contig;
awk 'FNR==NR { array[$1]=$2; next } { for (i in array) gsub(i, array[i]) }1' UBE05_R004.2_contig UBE05_R004.1_contig | sed -Ee 's/~.+~/none/g' | sed '1d' | sed -E 's/Group_2|Group_3|Group_4|Group_5/Group_2+/g' | sort -V -k 1,4 > UBE05_R004_contig;
awk -F "\t" '{if($2-$3 >= 0.1) print $1"\t"$4}' UBE05_R004_contig | sort -V | uniq -c | sed -E 's/ +([0-9]+) (.+)-----(.+)/\2\t\1\t\3/g' | awk -F"\t" 'BEGIN{OFS="\t"}; {print $3,$2 > "UBE05_R004.3_contig_"$4}';
awk -F "\t" '{if($2-$3 >= 0.1) print $1"\t"$4}' UBE05_R004_contig | sort -V | uniq -c | sed -E 's/ +([0-9]+) (.+)-----(.+)/\2\t\1\t\3/g' | awk -F"\t" '{print $3"\t"0}' | sort -Vu > UBE05_R004.4_contig
cat UBE05_R004.4_contig UBE05_R004.3_contig_none | sort -k 1,1V -k 2,2rn | awk '!a[$1]++' | sort -V -k 1,1 > UBE05_R005.1_contig;
cat UBE05_R004.4_contig UBE05_R004.3_contig_Group_2+ | sort -k 1,1V -k 2,2rn | awk '!a[$1]++' | sort -V -k 1,1 | cut -f 2 > UBE05_R005.2_contig;
cat UBE05_R004.4_contig UBE05_R004.3_contig_Group_1 | sort -k 1,1V -k 2,2rn | awk '!a[$1]++' | sort -V -k 1,1 | cut -f 2 > UBE05_R005.3_contig;
paste UBE05_R005.1_contig UBE05_R005.2_contig UBE05_R005.3_contig | sed -E 's/(.+)(_.ontig[0-9])/\1\t\1\2/g' | sort -k 1,1V -k 3,4nr | sed -z 's/^/## counts of genes with rbitscores virus > cell\ngenome\tcontig\ttotal_virus_genes\tgroup_2+_NCVOGs\tgroup_1_NCVOGs\n/g' > UBE05_R006_contig


