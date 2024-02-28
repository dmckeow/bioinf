#!/bin/bash
##### manual steps
### get NCVOG fastas:
## find members of orthogroups from NCVOG_NCVOGdef >> UBE06_001
#awk -F'\t' '{print > "UBE06_NN_"$10"_002"}' UBE06_001; # split by orthogroup names
#for file in UBE06_NN_*_002; do mv "$file" "${file// /_}"; done; for file in UBE06_NN_*; do mv "$file" "${file//-/_}"; done; for file in UBE06_NN_*; do mv "$file" "${file//(/}"; done; for file in UBE06_NN_*; do mv "$file" "${file//)/}"; done; # fix file names 
#for file in UBE06_NN_*_002; do awk -F'\t' '{print $1}' $file > "$file"_003; done; # get gi no. list
#for file in UBE06_NN_*_003; do grep -wf $file "/projet/fr2424/sib/dmckeown/db/virus/fastas/NCVOG.fa" | sed 's/>gi/gi/g' > "$file"_004; done;
#for file in UBE06_NN_*_004; do seqkit grep -n -f $file "/projet/fr2424/sib/dmckeown/db/virus/fastas/NCVOG.fa" > "$file"a; done; 

##### replace bad deflines in 004a:
###gi|9216 284504040 strand[272892..276245]1117_aa | DNA polymerase elongation subunit (family B) [Marseilevirus]
#with 
##gb|0000000|XXX|AVR53042.1| DNA polymerase delta catalytic subunit [Marseillevirus]
###gi|9212 441432435_441432434 gi|441432435|ref|YP_007354477.1| DNA polymerase type B [Acanthamoeba polyphaga moumouvirus] gi|441432434|ref|YP_007354476.1| DNA polymerase type B with an intein Hint (Hedgehog/Intein) domain [Acanthamoeba polyphaga moumouvirus]
#with
##gi|9212 441432435_441432434|ref|YP_007354477.1| DNA polymerase type B [Acanthamoeba polyphaga moumouvirus]
###gi|9213 448935108_448935109 gi|448935108|gb|AGE58659.1| DNA polymerase [Paramecium bursaria Chlorella virus NYs1] gi|448935109|gb|AGE58660.1| DNA polymerase [Paramecium bursaria Chlorella virus NYs1]
#with
##gi|9213 448935108_448935109 |gb|AGE58659.1| DNA polymerase [Paramecium bursaria Chlorella virus NYs1]

## fix 004_005 deflines to >TAXON|SEQ_ID:
#for file in UBE06_NN_*_004a; do sed '/>.*$/s/$/;/g' $file | sed 's/\[/|/g' | sed 's/\]/|/g' | tr -d '\n' | sed 's/>/\n>/g' | awk -F'|' '{print ">"$6,"|"$4,$7}' | sed 's/ /_/g' | sed '/^>_|_$/d' | sed 's/_;/;/g' > "$file"b; done; for file in UBE06_NN_*_004ab; do mv "$file" "${file//004ab/004_005}"; done; ## for multigene phylogeny
#for file in *004_005; do sed 's/;/\n/g' $file > "$file".fa; done; ## prep for single gene phylogeny (return to normal fasta format)

### get brown algal orthologs to NCVOGs fastas:
#for file in UBE06_NN_*_002; do awk -F'\t' '{print $7}' "$file" | uniq > "$file"_006; done;
#for file in UBE06_NN_*_006; do grep -f $file "/projet/fr2424/sib/dmckeown/phaeoex_screen/finalresult/UBE05/UBE05_002" > "$file"_007; done;
#for file in UBE06_NN_*_007; do awk -F'\t' '{print $2}' $file | sed 's/:1-/\t1\t/g'> "$file"_008.gff; done; # prep bed file for brown algal virus marker genes

#find /projet/fr2424/sib/dmckeown/phaeoex_screen/input/04__FINAL_PROTEOMES -name "*proteins.fa" > UBE06_NN_010; # brown algae_proteome_list
#find . -name "UBE06_NN_*_008.gff" | sed 's/.\///g'> UBE06_NN_011; # brown algal_viral marker genes_gff_for aa_list;

GFF1="UBE06_NN_Transcription_factor_S_II_TFIIS_002_006_007_008.gff"
GFF2="UBE06_NN_D5_like_helicase_primase_002_006_007_008.gff"
GFF3="UBE06_NN_DNA_polymerase_elongation_subunit_family_B_002_006_007_008.gff"
GFF4="UBE06_NN_NCLDV_major_capsid_protein_002_006_007_008.gff"
GFF5="UBE06_NN_packaging_ATPase_002_006_007_008.gff"
GFF6="UBE06_NN_Poxvirus_Late_Transcription_Factor_VLTF3_like_002_006_007_008.gff"

NFA1="UBE06_NN_Transcription_factor_S_II_TFIIS_002_003_004_005"
NFA2="UBE06_NN_D5_like_helicase_primase_002_003_004_005"
NFA3="UBE06_NN_DNA_polymerase_elongation_subunit_family_B_002_003_004_005"
NFA4="UBE06_NN_NCLDV_major_capsid_protein_002_003_004_005"
NFA5="UBE06_NN_packaging_ATPase_002_003_004_005"
NFA6="UBE06_NN_Poxvirus_Late_Transcription_Factor_VLTF3_like_002_003_004_005"

BFA1="UBE06_NN_Transcription_factor_S_II_TFIIS_002_006_007_008_009"
BFA2="UBE06_NN_D5_like_helicase_primase_002_006_007_008_009"
BFA3="UBE06_NN_DNA_polymerase_elongation_subunit_family_B_002_006_007_008_009"
BFA4="UBE06_NN_NCLDV_major_capsid_protein_002_006_007_008_009"
BFA5="UBE06_NN_packaging_ATPase_002_006_007_008_009"
BFA6="UBE06_NN_Poxvirus_Late_Transcription_Factor_VLTF3_like_002_006_007_008_009"

#while read file; do bedtools getfasta -fi $file -bed $GFF1 -fo >> "$GFF1".fa; done < UBE06_NN_010; while read file; do bedtools getfasta -fi $file -bed $GFF2 -fo >> "$GFF2".fa; done < UBE06_NN_010; while read file; do bedtools getfasta -fi $file -bed $GFF3 -fo >> "$GFF3".fa; done < UBE06_NN_010; while read file; do bedtools getfasta -fi $file -bed $GFF4 -fo >> "$GFF4".fa; done < UBE06_NN_010; while read file; do bedtools getfasta -fi $file -bed $GFF5 -fo >> "$GFF5".fa; done < UBE06_NN_010; while read file; do bedtools getfasta -fi $file -bed $GFF6 -fo >> "$GFF6".fa; done < UBE06_NN_010; for file in UBE06_NN_*_008.gff.fa; do mv "$file" "${file//.gff.fa/a}"; done;

## fix 008_009 deflines to >TAXON|SEQ_ID:
#for file in UBE06_NN_*_008a; do sed '/>.*$/s/$/;/g' $file | sed 's/contig/|contig/gi' | sed 's/mRNA.//g' | tr -d '\n' | sed 's/>/\n>/g' | awk -F';' '{print $1,";"$2}' | sed 's/ /_/g' | sed '/^_;$/d' | sed 's/_|/|/g' | sed 's/_;/;/g' > "$file"b; done; for file in UBE06_NN_*_008ab; do mv "$file" "${file//008ab/008_009}"; done; ## for multigene phylogeny
#####MANUAL: now remove badly aligned/short hits from brown algal genomes from 008_009
#for file in *008_009; do sed 's/;/\n/g' $file | sed -E 's/:|'|,/_/g' > "$file".fa; done;

##### concatenate (cat) NCVOG and phaeoexplorer NCLDV core genes for single gene phylogeny
#cat $NFA1.fa $BFA1.fa > "$BFA1"_010.fa; cat $NFA2.fa $BFA2.fa > "$BFA2"_010.fa; cat $NFA3.fa $BFA3.fa > "$BFA3"_010.fa; cat $NFA4.fa $BFA4.fa > "$BFA4"_010.fa; cat $NFA5.fa $BFA5.fa > "$BFA5"_010.fa; cat $NFA6.fa $BFA6.fa > "$BFA6"_010.fa; for file in *010.fa; do sed -i 's/_|/|/g' $file | sed -E 's/:|'|,/_/g' ; done;
##### multigene phylogeny will remove multicopy genes per genome #####



