#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=300GB
#SBATCH --partition ag2tb
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

echo -e "COMMAND LINE JOB SUBMISSSION:\n\tsbatch /home/dcschroe/dmckeow/projects/DWV/script/datasync2.sh $@"


########## 002 - make copy of all fastq_pass folders to backup shared MSI folder #################
date=`date +%d_%b_%Y`
time=`date +%k:%M:%S`
SYNCLOG2="/home/dcschroe/shared/data/gridion_datasynclog2"
touch $SYNCLOG2

echo -e "###########################\n MSI data sync started (MSI scratch storage to shared storage): $date at $time ###########################\n" >> $SYNCLOG2

find /scratch.global/dcschroe/data -name "fastq_pass" | sed 's/.*data\///g' > /scratch.global/dcschroe/datasynclog_tmp1
find /scratch.global/dcschroe/data -name "sequencing_summary*" | sed 's/.*data\///g' >> /scratch.global/dcschroe/datasynclog_tmp1
find /scratch.global/dcschroe/data -name "report_*.html" | sed 's/.*data\///g' >> /scratch.global/dcschroe/datasynclog_tmp1

rsync --exclude Install_logs --exclude intermediate --log-file /scratch.global/dcschroe/datasynclog_tmp2 -PrazR --files-from /scratch.global/dcschroe/datasynclog_tmp1 /scratch.global/dcschroe/data/ /home/dcschroe/shared/data/

cat /scratch.global/dcschroe/datasynclog_tmp2 >> $SYNCLOG2
rm -f /scratch.global/dcschroe/datasynclog_tmp*

chmod 777 -R /home/dcschroe/shared/data/


########## 002 - test fastq files for not corrupted #################

find /scratch.global/dcschroe/data/ -path "*fastq_pass*.fastq.gz" -type f -exec md5sum {} + | sed -E 's/ +/\t/g' | sed 's/\t.*\/data\//\tdata\//g' | sort -k 2,2V > /scratch.global/dcschroe/datasync_check1
find /scratch.global/dcschroe/data/ -path "*fastq_pass*.fastq.gz" -type f -exec md5sum {} + | sed -E 's/ +/\t/g' | sed 's/\t.*\/data\//\tdata\//g' | sort -k 2,2V > /scratch.global/dcschroe/datasync_check2

echo -e "############ \n MSI data sync md5 checksum comparison (MSI scratch storage to shared storage): ###########\nNumber mismatching md5checksums between new data in scratch vs new data transferred to shared (fastq.gz only):\n" >> $SYNCLOG2
diff /scratch.global/dcschroe/datasync_check1 /scratch.global/dcschroe/datasync_check2 | grep -c '^<' >> $SYNCLOG2
echo -e "############ \n files that may not have transferred properly:\n" >> $SYNCLOG2
diff /scratch.global/dcschroe/datasync_check1 /scratch.global/dcschroe/datasync_check2 | grep '^<' >> $SYNCLOG2

########## 003 - synchronise ALL DATA with 2ndtier storage #################


s3cmd sync --progress --no-delete-removed /scratch.global/dcschroe/data s3://dcschroe



E=$(diff /scratch.global/dcschroe/datasync_check1 /scratch.global/dcschroe/datasync_check2 | grep -c '^<')
if [[ $E -le 0 ]]; then
rm -fr /scratch.global/dcschroe/data/
fi
