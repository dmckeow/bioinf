#!/bin/bash
#
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#
#
#SBATCH --partition amdsmall
#SBATCH --cpus-per-task 2
#SBATCH --mem 4GB

####################### INPUTS #####################################
###### load software ######

##### command line arguments required #####
## $1 = your plain file of PCR results after STEP 000 (MANUAL)
## $2 = sample group column header name that will be the x axis labels - this is category by which your data will be grouped and percentages calculated

####################################################

#################### STEP 000 (MANUAL) #################################
####### BEFORE RUNNING SCRIPT DO THESE STEPS:
## 1. In Q-Rex, manually adjust bins and ensure appropriate genotypes are categorised as bins
## 2. Ctrl-a copy paste all fields from Absolute Quantification followed by Melt Curve (append them together left to right)
## 3a. Manually check all samples with the Melt Curve Flag of "peaks without bin". If sample actually does have a peak within or near a bin, replace Flag with "--" and add to the field Comment "manually corrected flag". The Flag field indicates positive or negative PCR result: "--" = positive, "Multiple peaks" = positive, with multiple variants, "No peak" = negative, "Peaks without bin" = negative.
## 3b. Manually check the Bin peaks and Add. peaks fields for multiple peaks per sample. If multiple peaks look real, change Flag to "Multiple peaks" and change the genotype - e.g. "A_and_B"
## 3c. Manually check for false negatives/positives. Any sample with low take-off and flat amplification curve should be negative, even if it has a Cq and/or genotype. Also check samples that have a Cq value but no genotype and vice versa. Change all false Cq, genotype, or multiple peak temperatures to "--".
## 4. In spreadsheet, find and replace all "\n" with ", " (regular expression needs to be on)
## 5. Copy this data onto the corresponding info for the samples.
## 6. If needed, make any new columns needed by editnig or by merging column(s) together (excel formula =D2&" "&E2) according to the sample groups required for figure - e.g. merge location and colony if you want figure to summarise the data in these groups, or leave it as it is if you plan to summarise the data by location or colony only.
## 7. Copy this data into a plain file (including all sample and PCR info).

##### fix any non-unix formatting and replace all characters with underscores that aren't letters, digits, periods, tabs, dash

a="file"; dos2unix $a; sed -z 's/$/\n/1' $a > "$a"_rtpcr1; sed -i '/^$/d' "$a"_rtpcr1; sed -i -E 's/[^A-Za-z0-9\t\.]/_/g' "$a"_rtpcr1; sed -i -E 's/_+/_/g' "$a"_rtpcr1; sed -i -E 's/_$|^_//g' "$a"_rtpcr1; sed -i -E 's/\t_|_\t/\t/g' "$a"_rtpcr1;
