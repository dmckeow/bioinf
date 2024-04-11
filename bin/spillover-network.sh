#!/bin/bash

umask 002


####################### SOFTWARE, DATABSES, ETC #####################################
### LOAD software available via shared environment on server:
eval "$(conda shell.bash hook)"
conda activate bioinftools

########################################################################
########################################################################

step="a4"
## first run phylogeny
## run in: /home/dcschroe/dmckeow/data/Spillover_FINAL/phylogeny

####################### THE SCRIPT #####################################

#################################################################


########################################## run mmseqs to clsuter
if [[ "${step}" =~ "A" ]] || [[ "$step" =~ "a1" ]]; then
###########################################
rm -fr mmseqs; mkdir mmseqs
for f in *.fa; do
    N=$(basename $f .fa)
    mmseqs easy-cluster $f mmseqs/${N}.mmseqs mmseqs/tmp.${N}.mmseqs --min-seq-id 0.3 -c 0.5

### give mmseqs clusters group names
awk -F "\t" '{if($1 == $2) print $1}' mmseqs/${N}.mmseqs_cluster.tsv > mmseqs/${N}.mmseqs.clusters
counter=1
while IFS=$'\t' read -r line; do
    group=$(printf '%03d' $counter)
    echo "/^$line\\t/ s/$/\t$group/g"
    ((counter++))
done < mmseqs/${N}.mmseqs.clusters | sed -f - mmseqs/${N}.mmseqs_cluster.tsv > mmseqs/tmp.${N}.mmseqs.clusters && mv mmseqs/tmp.${N}.mmseqs.clusters mmseqs/${N}.mmseqs.clusters

cut -f 3 mmseqs/${N}.mmseqs.clusters | sort -V | uniq -c | sed -E 's/^ +(.+) (.+)/\/\\t\2\$\/ s\/$\/\\t\1\/g/g' | sed -f - mmseqs/${N}.mmseqs.clusters | cut -f 2-4 > mmseqs/tmp.${N}.mmseqs.clusters && mv mmseqs/tmp.${N}.mmseqs.clusters mmseqs/${N}.mmseqs.clusters

#sed -i -z 's/^/Contig\tMmseqCluster\tNumInMmseqCluster\n/g' mmseqs/${N}.mmseqs.clusters

done

#################################################
fi
#################################################
########################################## run fastANI
if [[ "${step}" =~ "A" ]] || [[ "$step" =~ "a2" ]]; then
###########################################

##### now split fasta into mmseqs cluster by files
rm -fr mmseqs/*.mmseqclustered

for f in mmseqs/*.mmseqs.clusters; do
    while IFS=$'\t' read -r Contig MmseqCluster NumInMmseqCluster; do
        N=$(basename $f .mmseqs.clusters)
        if [ "$NumInMmseqCluster" -gt 1 ]; then
            mkdir -p mmseqs/${N}.${MmseqCluster}.mmseqclustered
            seqkit grep -p $Contig ${N}.fa -o mmseqs/${N}.${MmseqCluster}.mmseqclustered/${Contig}.fa
        fi
    done < $f
done

#################################################
fi
#################################################
########################################## run fastANI for mmseqs clustered sequences
if [[ "${step}" =~ "A" ]] || [[ "$step" =~ "a3" ]]; then
###########################################
rm -fr fastANI; mkdir fastANI
rm -f mmseqs/*.mmseqclustered.list
for d in mmseqs/*.mmseqclustered; do
    N=$(basename $d .mmseqclustered)
    realpath $d/*.fa > mmseqs/${N}.mmseqclustered.list
done

## fastANI
for f in mmseqs/*.mmseqclustered.list; do
    fastANI --ql $f --rl $f -o fastANI/$(basename $f .mmseqclustered.list) --matrix --fragLen 200
done

#################################################
fi
#################################################
########################################## run fastANI without mmseqs clustered sequences
if [[ "${step}" =~ "A" ]] || [[ "$step" =~ "a4" ]]; then
###########################################
rm -fr fastANI/*.NoMmseqs*

for f in *.fa; do
    seqkit split -i --by-id-prefix "" -O fastANI/$(basename $f .fa).NoMmseqs $f
done

for f in fastANI/*.NoMmseqs; do
    realpath $f/*.fa > ${f}.list
done

## fastANI
for f in fastANI/*.NoMmseqs.list; do
        fastANI --ql $f --rl $f -o fastANI/$(basename $f .NoMmseqs.list).NoMmseqs --matrix --fragLen 200
done

#################################################
fi
#################################################
