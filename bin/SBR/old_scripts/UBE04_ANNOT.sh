#!/bin/bash
## Make list of all UBE genes with best BLAST hit with virus bias each for NT, VNT, AA, and VAA in 1 row per gene and append the NCVOG only BLAST results
##### CL variable inputs:
## $1 BLAST_species_NCV
## $2 species
#### make NCVOGs def list:
## cut -f7 NCVOG > NCVOG_list; awk -F'\t' 'NR==FNR{a[$1]=$0;next} ($0 in a){print a[$0]}' NCVOGdef NCVOG_list > NCVOGdef2; paste NCVOG NCVOGdef2 > NCVOG_NCVOGdef;


cut -f13 $1 > UBE04_BLAST_NCV_list; sed -i 's/ /|/g' UBE04_BLAST_NCV_list; awk -F'|' '{print $2}' UBE04_BLAST_NCV_list > UBE04_BLAST_NCV_list2; awk -F'\t' 'NR==FNR{a[$1]=$0;next} ($0 in a){print a[$0]}' ../NCVOG_NCVOGdef UBE04_BLAST_NCV_list2 > UBE04_BLAST_NCV2; paste $1 UBE04_BLAST_NCV2 > UBE04_BLAST_NCV3; cut -f1-13 UBE04_BLAST_NCV3 > UBE04_BLAST_NCV4; awk -F'\t' '{print $22 $23}' OFS='|' UBE04_BLAST_NCV3 > UBE04_BLAST_NCV5; paste UBE04_BLAST_NCV4 UBE04_BLAST_NCV5 > /projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/UBE_fastas/BLAST_UBE04_$2_UBE_NCVOG;

sort -V -u -k 1,1 /projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/UBE_fastas/BLAST_UBE04_$2_UBE_NCVOG | cut -f1 > UBE04_NCV005; grep -Ff UBE04_NCV005 UBE03_AA004 | sort -V -k 1,1 > UBE04_NCV006; grep -vF -f UBE04_NCV005 UBE03_AA004 | sort -V -k 1,1 > UBE04_NCV007;   ### extracts UBE genes with NO BLAST hits 

awk -F'\t' '$31 = $31 FS "D_D_D"' OFS='\t' UBE04_NCV007 > UBE04_NCV007A;sort -V -u -k 1,1 /projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/UBE_fastas/BLAST_UBE04_$2_UBE_NCVOG > UBE04_NCV008; awk -F'\t' '$15 = $15 FS "A_A_A"' OFS='\t' UBE04_NCV008 > UBE04_NCV008A; sed -i 's/;START//g' UBE04_NCV008A; sed -i 's/;END//g' UBE04_NCV008A; paste UBE04_NCV006 UBE04_NCV008A > UBE04_NCV009; cat UBE04_NCV009 UBE04_NCV007A | sort -V -k 10,10 > UBE04_NCV010; ## remove all but 1 virus hits per gene; and label them; append UBE genes and BLAST results together; add UBE genes with NO BLAST hits

for file in xx*; do awk '{print FILENAME}' "$file" >> UBE04_NCV014; awk '{print $0}' "$file" >> UBE04_NCV015; done; paste UBE04_NCV014 UBE04_NCV015 > UBE04_NCV016; sed -i 's/;START//g' UBE04_NCV016; sed -i 's/;END//g' UBE04_NCV016; sed -i 's/;/\t/g' UBE04_NCV016; sort -u -V -k 10,10 UBE04_NCV016 > UBE04_NCV017; paste UBE04_NCV010 UBE04_NCV017 > UBE04_NCV018; cut -f 1-33 UBE04_NCV018 > UBE03_FINAL_GENE_NCVOG;


paste UBE03_FINAL_GENE_NT UBE03_FINAL_GENE_VNT UBE03_FINAL_GENE_AA UBE03_FINAL_GENE_VAA UBE03_FINAL_GENE_NCVOG > UBE04_FINAL_GENE_LABEL001; cat /projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/UBE03_HEADINGS_A_A_Ax5 UBE04_FINAL_GENE_LABEL001 > UBE04_FINAL_GENE_LABEL002; awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$67,$5,$6,$10,$14,$16,$132,$18,$51,$84,$117,$150,$19,$52,$85,$118,$151,$20,$53,$86,$119,$152,$27,$60,$93,$126,$159,$28,$61,$94,$127,$160,$29,$62,$95,$128,$161,$31,$64,$97,$130,$162,$32,$65,$98,$131,$164}' UBE04_FINAL_GENE_LABEL002 > /projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/UBE_fastas/UBE04_FINAL_ANNOT_$2;

