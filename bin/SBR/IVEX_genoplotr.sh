#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition fast
#SBATCH --cpus-per-task 4
#SBATCH --mem 8GB

### command line
## $1 = reduced version of IVEX004_[genome].gff (see STEP_000)
## $2 = dataframe (df) number (record which genome is which df)
## $3 = coords file with k option (generated by promer - must have ran IVEX_promer.sh)

###### MANUAL STEPS AND PREREQUISITES ######
## The purpose of this script is to convert promer alignments tabular outputs (coords) into inputs for genoplotR (in R) to generate a whole genome alignment visualisation
## this requires prior identification of specific viral regions to be shown aligned with viral reference IVEX000_genomes
## it is recommended that these genoplotR figures only be generated for a small set of viral regions per genome
## there is no point in running this script on genomes lacking long/complex viral sequences
## note that this script will also need to be run on reference genome(s)
## the dataframe (df) number is to be decided by you, and it determines the order of the alignment in the figure, e.g.:
## in this example, I am visualising each algal EVE between 2 reference viral genomes in a repeating pattern
## therefore we will run STEP 001 to generate each df as required by our desired output
## df1 Feldmannia spcies virus 158
## df2 Porterinema fluviatile EVE 7
## df3 Ectocarpus siliculosus virus 1
## df4 Porterinema fluviatile EVE 12
## df5 Feldmannia spcies virus 158
### in addition, STEP 002 will generate a comparison file for every pair of dataframes, e.g.:
### you will need to alter STEP_002 to match the desired figure arrangement - (change variables dfA dfB and the number of comparison files)
## note that STEP_002 generates all comparison file for all df pairs in a single run (delete comparison1 before trying to rerun STEP_002)
## df1 Feldmannia spcies virus 158
## comparison1.tab
## df2 Porterinema fluviatile EVE 7
## comparison2.tab
## df3 Ectocarpus siliculosus virus 1
## comparison3.tab
## df4 Porterinema fluviatile EVE 12
## comparison4.tab
## df5 Feldmannia spcies virus 158


### run from project folder (which contains archives input script tmp finalresult)

##### NCVOG defintion files #####
NCVKEY="/projet/fr2424/sib/dmckeown/phaeoex_screen/input/IVEX_NCVOG_006"

####### STEP_000 ###########
##### before running, in termnial, prepare a copy of the IVEX004.gff to cover only the viral region (so that loci start from 1 in both)
##### repeat this for every viral region - you will be running this script separately on every viral region that you want to be visualised - only use large viral sequences
#s="3498000"; e="3815000"; awk '/^P-fluviatile_contig7\t/' IVEX004_Porterinema-fluviatile.gff | awk -F"\t" -v s=$s -v e=$e 'BEGIN{OFS="\t"};{if($3 ~ "mRNA" && $4 >= s && $5 <= e) print $1,$2,$3,$4-s,$5-s,$6,$7,$8,$9,$10}'

############# STEP_001 ###############
######## PREPARE DATAFRAME THAT WILL BECOME DNA SEGMENT (repeat per genome/molecule involved)
awk -F"\t|;" 'BEGIN{OFS="\t"};{if($3 ~ "mRNA") print $22,$4,$5,$7,"",$18";"$21}' $1 | sed -e 's/ncvog_abbrev=//g' -e 's/\t+\t/\t1\t/g' -e 's/\t-\t/\t-1\t/g' > tmp0_df"$2"

awk -F"\t" '{print $5";"$7}' $NCVKEY | sed '/ncvog_group=Not_assigned/d' | sort -Vu | sort -V -t $'\t' -k 1,1 | sed 's/\t/;/g' > tmp1_df"$2"

##### key genes will be coloured with no label
sed -i -E 's/(ncvog_id=0023;ncvog_group=Group_1)/s\/\1\/#7a0c02\/g/g' tmp1_df"$2" ## D5-like_helicase-primase
sed -i -E 's/(ncvog_id=0031;ncvog_group=Group_1)/s\/\1\/#ac1600\/g/g' tmp1_df"$2" ## DNA_or_RNA_superfamily_II_helicase
sed -i -E 's/(ncvog_id=0076;ncvog_group=Group_1)/s\/\1\/#ac1600\/g/g' tmp1_df"$2" ## DNA_or_RNA_superfamily_II_helicase
sed -i -E 's/(ncvog_id=0038;ncvog_group=Group_1)/s\/\1\/#d23104\/g/g' tmp1_df"$2" ## DNA_polymerase_beta_subunit
sed -i -E 's/(ncvog_id=0052;ncvog_group=Group_1)/s\/\1\/#ec530e\/g/g' tmp1_df"$2" ## Disulfide_thiol_oxidoreductase
sed -i -E 's/(ncvog_id=1164;ncvog_group=Group_1)/s\/\1\/#fb8222\/g/g' tmp1_df"$2" ## Late_transcription_factor_2
sed -i -E 's/(ncvog_id=0022;ncvog_group=Group_1)/s\/\1\/#fdaf35\/g/g' tmp1_df"$2" ## Major_capsid_protein
sed -i -E 's/(ncvog_id=0249;ncvog_group=Group_1)/s\/\1\/#ecd23a\/g/g' tmp1_df"$2" ## Packaging_ATPase
sed -i -E 's/(ncvog_id=0304;ncvog_group=Group_1)/s\/\1\/#caee35\/g/g' tmp1_df"$2" ## Serine_threonine_protein_kinase
sed -i -E 's/(ncvog_id=0241;ncvog_group=Group_2)/s\/\1\/#a0fa3e\/g/g' tmp1_df"$2" ## Proliferating_cell_nuclear_antigen
sed -i -E 's/(ncvog_id=1353;ncvog_group=Group_2)/s\/\1\/#68fa67\/g/g' tmp1_df"$2" ## Ribonucleoside_diphosphate_reductase_alpha_subunit
sed -i -E 's/(ncvog_id=0276;ncvog_group=Group_2)/s\/\1\/#68fa67\/g/g' tmp1_df"$2" ## Ribonucleoside_diphosphate_reductase_beta_subunit
sed -i -E 's/(ncvog_id=0010;ncvog_group=Group_3)/s\/\1\/#32f297\/g/g' tmp1_df"$2" ## BRO_family_protein
sed -i -E 's/(ncvog_id=0262;ncvog_group=Group_3)/s\/\1\/#19ddc1\/g/g' tmp1_df"$2" ## Late_transcription_factor_3
sed -i -E 's/(ncvog_id=0278;ncvog_group=Group_3)/s\/\1\/#28bdea\/g/g' tmp1_df"$2" ## RuvC_Holliday_junction_resolvase
sed -i -E 's/(ncvog_id=0317;ncvog_group=Group_3)/s\/\1\/#4294ff\/g/g' tmp1_df"$2" ## Thioredoxin
sed -i -E 's/(ncvog_id=0225;ncvog_group=Group_4)/s\/\1\/#456be3\/g/g' tmp1_df"$2" ## Class_3_lipase
sed -i -E 's/(ncvog_id=1087;ncvog_group=Group_4)/s\/\1\/#4040a2\/g/g' tmp1_df"$2" ## Papain-like_cysteine_peptidase_Cathepsin_B_group
sed -i -E 's/(ncvog_id=0246;ncvog_group=Group_4)/s\/\1\/#31123b\/g/g' tmp1_df"$2" ## Ulp1_protease
sed -i -E 's/(ncvog_id=1192;ncvog_group=Group_4)/s\/\1\/#999999\/g/g' tmp1_df"$2" ## YqaJ_viral_recombinase
sed -i -E 's/(ncvog_id=1076;ncvog_group=Group_5)/s\/\1\/#000000\/g/g' tmp1_df"$2" ## Integrase_recombinase
sed -i -E 's/(ncvog_id=1175;ncvog_group=Group_5)/s\/\1\/#000000\/g/g' tmp1_df"$2" ## Integrase_resolvase
sed -i -E 's/(ncvog_id=0077;ncvog_group=Group_5)/s\/\1\/#cccccc\/g/g' tmp1_df"$2" ## Sensor_histidine_kinase
sed -i '/^ncvog_id/d' tmp1_df"$2"

sed -f tmp1_df"$2" tmp0_df"$2" > tmp2_df"$2"
sed -i -E -e 's/\([0-9]\)\t/\t/g' -e 's/\tncvog_id=.+/\t#696969/g' tmp2_df"$2"
awk -F"\t" 'BEGIN{OFS="\t"};{if($6 !~ "#696969") gsub(/.*/,"-",$1)}1' tmp2_df"$2" | sed 's/^-\t/\t/g' | sed -z 's/^/name\tstart\tend\tstrand\tcol\tfill\n/1' > df"$2"


##### STEP_002 #############
##### PREPARE COMPARISON FILE FROM PROMER
### see IVEX_promer.sh for prerequisites
### these steps remove duplicated reciprocal matches and realgin to dfA (e.g. df1) vs dfB (e.g. df2) comparison (or df3 df4 and so on)
## dfA dfB comparison and dfB dfA comparison, remove duplicates
## dfA dfB are the query and reference in the promer coords file, ensure you use variables that match exactly the names of your queries/references in the coords file and matching the correct pair ; e.g. if df1 is FsV-158 and df2 is 7, ensure comparison1 is only the promer alignment coords between these 2 molecules
## duplicate and as neeeded and change variables f to .coords file, and dfA dfB to pair of molecular names
if [[ ! -e ./comparison1.tab ]]; then
dfA="FsV-158"; dfB="7"; awk -F"\t" 'BEGIN{OFS="\t"};{if($12 != $13) print $1,$2,$3,$4,"grey",$12,$13,$7}' $3 | grep -P "\t$dfA\t$dfB\t" > tmp1; awk -F"\t" 'BEGIN{OFS="\t"};{if($12 != $13) print $3,$4,$1,$2,"grey",$13,$12,$7}' $f | grep -P "\t$dfA\t$dfB\t" > tmp2; cat tmp1 tmp2 | sort -Vu | cut --complement -f 6,7 | sed -z 's/^/start1\tend1\tstart2\tend2\tcol\tvalues\n/1' > comparison1.tab
dfA="7"; dfB="EsV-1"; awk -F"\t" 'BEGIN{OFS="\t"};{if($12 != $13) print $1,$2,$3,$4,"grey",$12,$13,$7}' $3 | grep -P "\t$dfA\t$dfB\t" > tmp1; awk -F"\t" 'BEGIN{OFS="\t"};{if($12 != $13) print $3,$4,$1,$2,"grey",$13,$12,$7}' $f | grep -P "\t$dfA\t$dfB\t" > tmp2; cat tmp1 tmp2 | sort -Vu | cut --complement -f 6,7 | sed -z 's/^/start1\tend1\tstart2\tend2\tcol\tvalues\n/1' > comparison2.tab
dfA="EsV-1"; dfB="12"; awk -F"\t" 'BEGIN{OFS="\t"};{if($12 != $13) print $1,$2,$3,$4,"grey",$12,$13,$7}' $3 | grep -P "\t$dfA\t$dfB\t" > tmp1; awk -F"\t" 'BEGIN{OFS="\t"};{if($12 != $13) print $3,$4,$1,$2,"grey",$13,$12,$7}' $f | grep -P "\t$dfA\t$dfB\t" > tmp2; cat tmp1 tmp2 | sort -Vu | cut --complement -f 6,7 | sed -z 's/^/start1\tend1\tstart2\tend2\tcol\tvalues\n/1' > comparison3.tab
dfA="12"; dfB="FsV-158"; awk -F"\t" 'BEGIN{OFS="\t"};{if($12 != $13) print $1,$2,$3,$4,"grey",$12,$13,$7}' $3 | grep -P "\t$dfA\t$dfB\t" > tmp1; awk -F"\t" 'BEGIN{OFS="\t"};{if($12 != $13) print $3,$4,$1,$2,"grey",$13,$12,$7}' $f | grep -P "\t$dfA\t$dfB\t" > tmp2; cat tmp1 tmp2 | sort -Vu | cut --complement -f 6,7 | sed -z 's/^/start1\tend1\tstart2\tend2\tcol\tvalues\n/1' > comparison4.tab
dfA="FsV-158"; dfB="22"; awk -F"\t" 'BEGIN{OFS="\t"};{if($12 != $13) print $1,$2,$3,$4,"grey",$12,$13,$7}' $3 | grep -P "\t$dfA\t$dfB\t" > tmp1; awk -F"\t" 'BEGIN{OFS="\t"};{if($12 != $13) print $3,$4,$1,$2,"grey",$13,$12,$7}' $f | grep -P "\t$dfA\t$dfB\t" > tmp2; cat tmp1 tmp2 | sort -Vu | cut --complement -f 6,7 | sed -z 's/^/start1\tend1\tstart2\tend2\tcol\tvalues\n/1' > comparison5.tab
dfA="22"; dfB="EsV-1"; awk -F"\t" 'BEGIN{OFS="\t"};{if($12 != $13) print $1,$2,$3,$4,"grey",$12,$13,$7}' $3 | grep -P "\t$dfA\t$dfB\t" > tmp1; awk -F"\t" 'BEGIN{OFS="\t"};{if($12 != $13) print $3,$4,$1,$2,"grey",$13,$12,$7}' $f | grep -P "\t$dfA\t$dfB\t" > tmp2; cat tmp1 tmp2 | sort -Vu | cut --complement -f 6,7 | sed -z 's/^/start1\tend1\tstart2\tend2\tcol\tvalues\n/1' > comparison6.tab
fi
