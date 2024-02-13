# bioinf
* A basic pipeline composed of custom scripts that run various publicly available tools.
* Main purpose is for common bioinformatics protocols including assembly and binning of sequence data, phylogeny, and other genomic analyses.
* Currently the pipeline is focused on using Oxford Nanopore Sequencing data, and on viral metagenomics

## INSTALLATION STEP 1 - install ANVI'O
* Anvi'o is the tool we will use for binning - the developers are the Merenlab, and they have very good guides for [anvi'o](https://anvio.org/).
* First thing is to create a conda environment with Anvi'o 8 installed within it
* Please follow the exact instructions for [installing anvio](https://anvio.org/install/linux/stable/), BUT with one **IMPORTANT EXCEPTION** - name your environment `bioinftools` instead of `anvio-8`. The pipeline will look for a conda environment named bioinftools:
* IMPORTANT - during your installation, be sure to follow their step that details how to setup key resources for anvi'o - currently this is under: (4.1) Setup [key resources](https://anvio.org/install/)

## INSTALLATION STEP 2 - get other tools and packages
* Copy the following into your terminal:
```shell
conda activate bioinftools ## your conda environment must be active to use it

mamba install -y -c bioconda kaiju
mamba install -y -c bioconda seqkit
mamba install -y -c bioconda hyphy
mamba install -y -c bioconda minimap2
mamba install -y -c bioconda porechop
mamba install -y -c conda-forge pigz
mamba install -y -c bioconda taxonkit
mamba install -y -c bioconda mafft
mamba install -y -c bioconda cd-hit
mamba install -y -c conda-forge -c bioconda mmseqs2
mamba install -y bioconda::canu
mamba install bioconda::bedtools

## R packages
mamba install -y -c bioconda r-svglite
mamba install -y -c conda-forge r-reshape
```
* Install packages within R environment
```shell
conda activate bioinftools
R
install.packages("RColorBrewer")
install.packages("dplyr")
install.packages("forcats")
install.packages("ggalluvial")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("plotrix")
install.packages("reshape")
install.packages("svglite")
install.packages("tidyr")
install.packages("tidyverse")
install.packages("treemapify")
install.packages("viridisLite")
install.packages("viridis")
```
* Possible issue with CANU and conda - if it does not work, then install as binary, as detailed [here](https://github.com/marbl/canu )
* Once it is installed you must add canu to your path by adding `export PATH="path/to/your/canu-2.2/bin:$PATH"` to your ~/.bashrc using an editor like vim or nano, save the file changes and then do `source ~/.bashrc`

## INSTALLATION STEP 3 - download and setup bioinf pipelines, scripts, and databases
`FOR a NEW SETUP and bioinfdb: run 3a. AND 3b.`
`If your system already has the databases setup (i.e. someone has already ran step 3a.), then run 3b. ONLY`

### 3a. For a fresh setup (creates new databases from scratch)
* If you are creating a fresh setup for bioinf, copy the following, but replace /path/to/temporary_directory with your system's path to the temporary_directory storage:
```shell
git clone https://github.com/dmckeow/bioinf.git
chmod +x bioinf/bin/*
sbatch --cpus-per-task=24 --time=96:00:00 --mem=240GB --partition your_partition -o slurm.%N.%j.out -e slurm.%N.%j.err bioinf/bin/bioinf-setup.sh -t /path/to/temporary_storage -d /where/to/create/your/database/directory
```
* The setup will now run, and it will probably take overnight to complete, because it is building some large databases such as nr for Kaiju. Once it is finished, do:
```shell
source ~/.bashrc ## OR logout and log back in
```
* Note that this example is for submitting to a SLURM cluster - the parameters before bioinf-setup.sh are SLURM and system specific - change them according to your own requirements.

### 3b. Create a directory that links to a preexisting database
* If you have access to a bioinfdb that has already been setup, then you should just run 3b - it will simply configure your version of bioinf to run using the pre-existing setup. This will complete within minutes and save 100s in GB of storage space.
* Creates symbolic links to pre-existing databases. This allows us to share a core bioinfdb that multiple users can use. Each user can then add/remove databases from their own database without affecting anyone elses database setup (EXCEPT that deleting files in the original database will remove them from all linked db copies)
* Copy the following code, but with your own paths:
```shell
git clone https://github.com/dmckeow/bioinf.git
chmod +x bioinf/bin/*
bioinf/bin/bioinf-setup.sh -t /path/to/temporary_storage -d /where/to/create/your/database/directory -E /path/to/existing/bioinfdb ## replace with your actual paths
source ~/.bashrc ## OR logout and log back in
```
* Temporary storage is specified because some scripts such as bioinf-assembly-canu.sh will generate quite large temporary files - so you should choose somewhere with enough capacity to handle 10s to 100s of GB of data. If you are on a SLURM cluster, use the scratch directory!

## And that is the installation and setup done!

___

# USAGE
* Within the bioinf/bin folder is a collection of scripts, written in bash, python, and R
* They should now be available through your PATH: in your terminal type bioinf and press TAB twice and you should see all of the scripts available
* All scripts send their output to the current directory from which the script was run
* Every script explains how it can be run if you give it the --help flag, e.g.:
```shell
bioinf-assembly-canu.sh -h
```

## RUNNING on a Computing Cluster, specifically SLURM
* Hopefully you are using a computing cluster, as many of these scripts use software that need a lot of computing resources
* If using a SLURM cluster using sbatch - you must provide the SLURM parameters
* The SLURM parameters below are only examples - you must provide ones which are suited to your available resources and system configuration
* You SHOULD specify a number of threads to use through the SLURM flag --cpus-per-task, otherwise the scripts will use just 1 thread, and some will take a million years to finish
* Use the script flag --threads to specify thread numbers on non-SLURM systems

#### WITH a bash (.sh) script in SLURM
* In this example we have a bash script called script_for_submitting_slurm_jobs.sh, which contains the following code:
```shell
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=240GB
#SBATCH --partition your_partition
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

## for a bash script:

bioinf-binning.sh -i /path/to/samplelistfile -p project_name

## for a python script:

eval "$(conda shell.bash hook)" ## this allows access to conda in the batch job
conda activate bioinftools ## the conda environment must be activated BEFORE the python script is submitted
bioinf-phylogeny.py -i /path/to/input_file -p project_name
```
* Then we submit the batch job script:
```shell
sbatch script_for_submitting_slurm_jobs.sh
```

___

#### WITHOUT a bash (.sh) script in SLURM
* **ALTERNATIVELY**, you can sbatch without a script, and do it through the command line instead, which will produce exactly the same result as using the script e.g.:
```shell
sbatch --time=24:00:00 --cpus-per-task=12 --mem=64GB --partition your_partition -o slurm.%N.%j.out -e slurm.%N.%j.err bioinf-binning.sh -i /path/to/samplelistfile -p project_name
```
* **NOTE** that this must be done differently for a python script (they end with .py), where we use the "--wrap" option of sbatch to activate conda before running the script
```shell
sbatch --time=24:00:00 --cpus-per-task=12 --mem=64GB --partition your_partition -o slurm.%N.%j.out -e slurm.%N.%j.err --wrap='eval "$(conda shell.bash hook)"; conda activate bioinftools; bioinf-phylogeny.py -i /path/to/input_file -p project_name'
```

___

### RUNNING without a Computing Cluster

* You can also run any script locally, i.e. through your local system
* Use the script flag --threads to specify thread numbers on non-SLURM systems
* Note that some of these pipelines require 10s of GB of RAM and multiple threads to complete in a reasonable time
* Running a script on a local system:
```shell
conda activate bioinftools ## the conda environment must be activate in your interactive session, for any script

## for a bash script:
bioinf-binning.sh -i /path/to/samplelistfile -p project_name --threads 12

## for a python script:
python bioinf-phylogeny.py -i /path/to/input_file -p project_name

```
* Please note that the bioinf-assembly-canu.sh script will NOT work on a non-SLURM system

___

# WHAT TO RUN
* **To see what scripts are available, in your terminal** `type bioinf and press TAB twice`
* For the most part, these scripts are modular and do not depend on you having run any other particular script. Exceptions to this will share common name structure - e.g. scripts beginning with `bioinf-binning` are part of a set workflow, i.e. you must run `bioinf-binning.sh` before running `bioinf-binning-summarise.sh`
### bioinf-assembly-canu.sh
### bioinf-binning.sh
### bioinf-binning-summarise.sh (WORK IN PROGRESS)
### bioinf-phylogeny.py

___

## OTHER NOTES
### Viruses and Anvi'o
* Anvi'o is primarily made for bacterial metagenomics, so there are several features and associated outputs that will not be useful to virology. Do not be concerned if any default files or summaries from anvi'o contain little useful information
* Are your viruses not in nice, neat bins? Try this - in the anvi-interactive interface, under the "Main tab", click "Items Order", and select "Sequence Composition (D: Euclidean; L: Ward)"

### Customising the databases
* All of the database directories are customisable - meaning that whatever you place there will be used by the pipelines (so long as it is the correct format and in the correct directory)
* No specific databases are required by the pipelines, but having NO databases in any of these database directories will cause some scripts to complain. Some scripts may not work if certain database categories were not used at all.
```shell
bioinfdb
    |------BLAST ## nr database files built by makeblastdb (see below)
    |
    |------DMND ## .dmnd database files built by diamond makedb (see below)
    |
    |------HMM ## subdirectories each containing a single hmmer profile built from a set of reference proteins (see below)
    |
    |------KAIJU ## subdirectories each containing a single kaiju database built by kaiju-makedb (see below)
    |
    |------VOGDB ## this is simply a copy of VOGDB. There is no reason to customise this

```
* Use soft links to databases files that already exist on your system to avoid space wastage - e.g. `ln -s /path/to/shared/nr.dmnd /path/to/my/bioinfdb/DMND/nr.dmnd`

#### BLAST
* Simple. By default, bioinftools does not build any BLAST databases (because Diamond is faster), but whatever blast database files or links you place in the BLAST folder will be used by the pipeline.
* See BLAST's instructions for making a database (makeblastdb)
```shell
ln -s /common/bioref/blast/latest/nr.* /home/dcschroe/dmckeow/bioinfdb/BLAST ## example of making a symbolic link to an existing blastdb
```

#### DMND
* Easy! A diamond database file is just a single .dmnd file. See [diamond](https://github.com/bbuchfink/diamond). Also:
```shell
conda activate bioinftools
diamond makedb -h
## OR
bioinf-setup.sh -h ## see STEP B1 for diamond makedb e.g.:

sbatch --cpus-per-task=12 --time=96:00:00 --mem=64GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err bioinf-setup.sh -t /scratch.global -d /home/dcschroe/dmckeow -s B1 --dmnd /path/to/protein.fasta
```
#### HMM
* Annoying! Unfortunately, anvi'o requires a specific and stupid set of 6 files to run a custom HMM as detailed [here](https://merenlab.org/2016/05/21/archaeal-single-copy-genes/):
```shell
genes.hmm.gz ## hmm profile generated from a sequence set using HMMER3
genes.txt ## 3 columns: gene, accession, hmmsource - the gene names must match the NAME info in genes.hmm.gz
kind.txt ## the database name. This WILL appear in your anvi'o outputs
noise_cutoff_terms.txt ## this sets the evalue for the HMM searches
reference.txt ## this is just a citation or whatever
target.txt ## sequence type e.g. AA:GENE
```
* Have a look at bioinfdb/HMM/VOGDB created by the bioinf-setup script and copy that format
* Each HMM db must have its own subdirectory within bioinfdb/HMM e.g. bioinfdb/HMM/VOGDB) that contains these 6 files and NOTHING ELSE. Seriously, anvi'o will not run the HMMs if it finds an unexpected file there.

#### KAIJU
* Easy! A kaiju database file is just 3 files. See [here](https://github.com/bioinformatics-centre/kaiju). Also:
```shell
conda activate bioinftools
kaiju-makedb -h
## OR
bioinf-setup.sh -h ## see STEP B2 for kaiju makedb e.g.:

sbatch --cpus-per-task=12 --time=96:00:00 --mem=64GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err bioinf-setup.sh -t /scratch.global -d /home/dcschroe/dmckeow -s B1 --kaiju fungi
```
* Kaiju can build predetermined databses OR build from sequences provided by the user: [here](https://github.com/bioinformatics-centre/kaiju)
* Each kaiju db must have its own subdirectory within bioinfdb/KAIJU e.g. bioinfdb/KAIJU/rvdb) that contains these 3 files:
```shell
names.dmp
nodes.dmp
<database_name>.fmi
```

#### VOGDB
* Unless you want to work with the viral gene info in this database, you should just ignore this database