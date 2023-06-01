#!/bin/bash

##### example for running script:
### you must provide your MSI login ($1)
### e.g.: bash /data/datasync1.sh dmckeow@agate.msi.umn.edu

a="$1"
date=`date +%d_%b_%Y`
time=`date +%k:%M:%S`
SYNCLOG1="/data/gridion_datasynclog1"
touch $SYNCLOG1

red='\033[0;31m'
NC='\033[0m'

########## 000 - gzip all .fastq files in data folder and synchronise gridion data folder (applies to EVERYTHING in this folder) with MSI server scratch folder (TEMPORARY storage)  ##############


find /data/gridion_* -name "*.fastq" > /home/grid/tmp1_datasync
for f in $(cat /home/grid/tmp1_datasync)
do
	pigz -p 12 $f
done

echo -e "###########################\n GridION data sync started (GridION to MSI scratch storage):\n $date at $time ###########################\n" >> $SYNCLOG1

echo -e "############ folders synced (GridION to MSI scratch storage): ###########\n" >> $SYNCLOG1
realpath /data/gridion_* >> $SYNCLOG1

realpath /data/gridion_* | sed 's/^\/data\///g' > /data/datasynclog_tmp1
rm -f /data/datasynclog_tmp2

rsync --log-file /data/datasynclog_tmp2 -Praz --files-from /data/datasynclog_tmp1 /data/ $a:/scratch.global/dcschroe/data/

echo -e "############ rsync tranfer log (GridION to MSI scratch storage): ###########\n" >> $SYNCLOG1
cat /data/datasynclog_tmp2 >> $SYNCLOG1
cat $SYNCLOG1
rm -f /data/datasynclog_tmp*
echo "Once logged in, do: sbatch /home/dcschroe/dmckeow/projects/DWV/script/datasync2.sh"
ssh -Y ${a}



