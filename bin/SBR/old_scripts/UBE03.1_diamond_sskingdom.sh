#!/bin/bash

##### RUN IN finalresult/UBE_fastas


##### CL INPUTS:
## $1 BLAST_...UBE_AA
## $2 /projet/fr2424/sib/dmckeown/db/NCBItaxIDsrankedlineage[n].dmp
## $3 UBE_genes_AA.fa input for BLAST
## $4 [species_name] 
##### Ensure diamond BLAST tabular columns are $1-12 default, $13 All Subject Title(s) [salltitles], $14 unique Subject Taxonomy ID(s), separated by a ';' (in numerical order) [staxids]
## if lines counts different then update .dmp file manually and then run again

##### @1 prepare sskingdom list to be pasted to AA Diamond hits vs nr#####

cut -f14 $1 > UBE03.1_TAXID_LIST1; cut -d";" -f1 UBE03.1_TAXID_LIST1 > UBE03.1_TAXID_LIST2;    # make taxid-sskingdom list + remove redundant species ID
awk -F'\t' 'NR==FNR{a[$1]=$0;next} ($0 in a){print a[$0]}' $2 UBE03.1_TAXID_LIST2 > UBE03.1_TAXID1; ##file2 file1 >> output
## the .dmp has a line which is ^\t\t$ so that blank lines are searched for and added to TAXID1

wc -l UBE03.1_TAXID_LIST2; wc -l UBE03.1_TAXID1; ## print to screen, line counts should be the same

##### @2 if line counts not the same: #####
## manually find missing TAXIDS in TAXID_MISSING_LIST and add them with kingdoms to the modified .dmp file (highest n.dmp)
## finally, run script again with the updated .dmp file

awk -F'\t' 'FNR==NR{a[$1];next}!($1 in a)' $2 UBE03.1_TAXID_LIST2 > UBE03.1_TAXID_MISSING; sort -V -u UBE03.1_TAXID_MISSING > UBE03.1_TAXID_MISSING_LIST; ## finds taxids in BLAST results but missing from the .dmp file (error by NCBI) and therefore the TAXID1 file

##### @3 add sskingdom to BLAST results #####
paste $1 UBE03.1_TAXID1 > UBE03.1_FIX4; cut -d$'\t' --complement -f 14 UBE03.1_FIX4 > BLAST_UBE03.1_$4_UBE_AA;


