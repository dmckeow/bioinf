#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition long
#SBATCH --cpus-per-task 8
#SBATCH --mem 32GB

################## run from PC

###############################
###############################

##### SBR_cluster > PC internal SSD (ALL)
#rsync -Praz --delete dmckeown@slurm0.sb-roscoff.fr:/projet/fr2424/sib/dmckeown/phaeoex_screen/ /mnt/c/Users/Dean\ Mckeown/Documents/PC_iSSD/SBR_backup/

##### SBR_cluster < PC internal SSD (ALL)
#rsync -Praz --delete /mnt/c/Users/Dean\ Mckeown/Documents/PC_iSSD/SBR_backup/ dmckeown@slurm0.sb-roscoff.fr:/projet/fr2424/sib/dmckeown/phaeoex_screen/

##### SBR_cluster < PC internal SSD (only for specific folders)
#rsync -Praz --delete /mnt/c/Users/Dean\ Mckeown/Documents/PC_iSSD/SBR_backup/script/ dmckeown@slurm0.sb-roscoff.fr:/projet/fr2424/sib/dmckeown/phaeoex_screen/script/

###############################
###############################

##### PC internal SSD to PC external HDD
#rsync -Praz --delete /mnt/c/Users/Dean\ Mckeown/Documents/PC_iSSD/ /mnt/d/eHDD_backup/

##### personal to shared folder
rsync -Praz --delete /projet/fr2424/sib/dmckeown/phaeoex_screen/ /shared/projects/phaeoexplorer_virus/phaeoex_screen/
