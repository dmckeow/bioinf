##### files needed and formatting

## genes.hmm.gz 	A gzip of concatenated HMM profiles.
HMMER3/f [3.1b2 | February 2015]
NAME  Enolase_C
ACC   PF00113.19
DESC  Enolase, C-terminal TIM barrel domain
LENG  296
ALPH  amino
RF    no
MM    no
CONS  yes
CS    yes
MAP   yes
DATE  Tue Oct 13 14:39:38 2015
NSEQ  9
EFFN  0.654785
CKSUM 2360876519
GA    21.50 21.50;
TC    21.50 21.50;
NC    21.40 21.40;
BM    hmmbuild HMM.ann SEED.ann
SM    hmmsearch -Z 11927849 -E 1000 --cpu 4 HMM pfamseq
STATS LOCAL MSV      -11.0260  0.70140
STATS LOCAL VITERBI  -11.7198  0.70140
STATS LOCAL FORWARD   -5.3917  0.70140

## genes.txt 	A TAB-delimited file containing three columns; gene name, accession number, and source of HMM profiles listed in genes.hmm.gz.
gene    accession       hmmsource
Enolase_C       PF00113 Rinke_2013
TIM     PF00121 Rinke_2013
ATP‐synt_C      PF00137 Rinke_2013
PGK     PF00162 Rinke_2013
Ribosomal_S12   PF00164 Rinke_2013
Ribosomal_S7    PF00177 Rinke_2013
Ribosomal_L2    PF00181 Rinke_2013

## kind.txt 	A flat text file which contains a single word identifying what type of profile the directory contains. If this word is ‘singlecopy’, the profile is used to calculate percent completeness and contamination. Otherwise it will only be used to visualize contigs with HMM hits without being utilized to estimate completeness.
singlecopy

## reference.txt 	A file containing source information for this profile to cite it properly.
Rinke et al, http://www.nature.com/nature/journal/v499/n7459/full/nature12352.html

## target.txt 	A file containing the target alphabet and context. See this for more details. For this particular collection the target will be AA:GENE (because it was prepared using amino acid alignments, and we want them to be searched within gene calls stored in contigs databases), however, the the target term could be any combination of AA, DNA, or RNA for the alphabet part, and GENE or CONTIG for the context part (well, except AA:CONTIG, because we can’t translate contigs).
AA:GENE

## noise_cutoff_terms.txt 	A file to specify how to deal with noise. See this comment for more information on the contents of this file.
-E 1e-12


########### What is in the RVDB hmm.gz?
HMMER3/f [3.3.1 | Jul 2020]
NAME  FAM010643
LENG  238
ALPH  amino
RF    no
MM    no
CONS  yes
CS    no
MAP   yes
DATE  Thu Jan  6 16:14:08 2022
NSEQ  4
EFFN  0.519531
CKSUM 1500063752
STATS LOCAL MSV      -10.6429  0.70341
STATS LOCAL VITERBI  -11.6154  0.70341
STATS LOCAL FORWARD   -5.2230  0.70341


#### what is in the RVDB annot files
cat 23296.txt
LENGTH  96
LCA     Viruses::unclassified viruses::unclassified DNA viruses::unclassified dsDNA viruses::Pandoravirus
NBSEQ   4
KEYWORDS:
--
KEYWORDS FROM SEQUENCES:
hypothetical    4
protein 4
pandoravirus    4
quercus 1
dulcis  1
salinus 1
celtis  1
SEQUENCES:
acc|REFSEQ|YP_009483803.1|REFSEQ|NC_037667|hypothetical protein [Pandoravirus quercus]
acc|REFSEQ|YP_008319930.2|REFSEQ|NC_021858|hypothetical protein [Pandoravirus dulcis]
acc|REFSEQ|YP_009430139.1|REFSEQ|NC_022098|hypothetical protein [Pandoravirus salinus]
acc|GENBANK|QBZ81709.1|GENBANK|MK174290|hypothetical protein [Pandoravirus celtis]




########### prepare genes.txt
### one entry per protein
for f in annot/*.txt; do awk -v f=$f '{print f";"$0}' $f | awk '/\|/' | sed 's/;/|/g' | cut -d "|" -f 1,4 | sed -e 's/annot\//000000000/g' -Ee 's/\.txt//g' -e 's/.*([0-9][0-9][0-9][0-9][0-9][0-9]\|.*)/FAM\1/g' | sed 's/|/\t/g' ; done > tmp_all_annots_RVDB
