# bioinf
* A basic pipeline composed of custom scripts that run various publicly available tools.
* Main purpose is for common bioinformatics protocols including assembly and binning of sequence data, phylogeny, and other genomic analyses.
* Currently the pipeline is focused on using Oxford Nanopore Sequencing data, and on viral metagenomics

## INSTALLATION STEP 1 - get anvi'o
* Anvi'o is the tool we will use for binning - the developers are the Merenlab, and they have very good guides for [anvi'o](https://anvio.org/), which is capable of many useful analyses

* First, you need to create a conda environment with Anvi'o 7.1 installed within it
* If you are on a server, it might be better to install conda (miniconda) on your own user profile (and mamba if you want faster installation - use conda install for this)
* Please follow the exact instructions for [installing anvio](https://anvio.org/install/), BUT with one IMPORTANT EXCEPTION - name your environment bioinftools instead of anvio-7.1. The pipeline will look for a conda environment named bioinftools:
* e.g. INSTEAD of
```
conda create -y --name anvio-7.1 python=3.6
```
* do
```
conda create -y --name bioinftools python=3.6
```
* If you have a problem with installation, the answer is probably on the anvi'o installation page.
* If certain packages install very slowly or timeout, try "mamba install" instead of "conda install"

* IMPORTANT - during your installation, be sure to follow their step that details how to setup key resources for anvi'o - currently this is under: (4.1) Setup [key resources](https://anvio.org/install/)

## INSTALLATION STEP 2 - get other tools
* Copy paste the following into your terminal
* Note that conda environments must be activated in order to install, change or use them and the packages within them:
```
conda activate bioinftools
```
* Copy paste the following into your terminal:
```
conda install -y -c bioconda kaiju
conda install -y -c bioconda seqkit
conda install -y -c bioconda hyphy
conda install -y -c bioconda minimap2
conda install -y -c bioconda porechop
conda install -y -c conda-forge pigz
conda install -y -c bioconda taxonkit
conda install -y -c bioconda canu
```


## INSTALLATION STEP 3 - download and setup bioinf pipelines, scripts, and databases

### For a new setup
* If you are creating a fresh setup for bioinf, copy the following, but replace /path/to/temporary_directory with your system's path to the temporary_directory storage:
```
cd $HOME
git clone https://github.com/dmckeow/bioinf.git
cd bioinf
sbatch --cpus-per-task=24 --time=96:00:00 --mem=240GB --partition ag2tb -o slurm.%N.%j.out -e slurm.%N.%j.err bioinf-setup.sh -t /path/to/temporary_directory
```
* The setup will now run, and it will probably take overnight to complete, because it is building some large databases such as nr for Kaiju
* Note that this example is for submitting to a SLURM cluster - the parameters before bioinf-setup.sh are SLURM and system specific - change them according to your own requirements.

### For a pre-existing setup
* If you or someone within your group has already setup bioinf, then maybe do this instead - it will simply configure your version of bioinf to run using the pre-existing setup
* This will be done in minutes, because it will only need to set the paths to the databases
* It will also save 100s of GB in storage space by avoiding having separate databases built for each user
* Creates soft links to files - so each user will work using their own bioinf directory copy. This means they can still add their own custom databases and they will not appear in anyone elses
* Copy the following code, but replace /path/to/temporary_directory with your system's path to the temporary_directory storage, AND replace /path/to/bioinf with the path to a pre-existing bioinf repository:
```
cd $HOME
git clone https://github.com/dmckeow/bioinf.git
cd bioinf
bash bioinf-setup.sh -t /path/to/temporary_directory -P /path/to/preexisting/bioinf/db
```
* Temporary storage is specified because some scripts such bioinf-assembly-canu.sh will generate very large temporary files - so you should choose somewhere with enough capacity to handle 10s to 100s of GB of data. If you are on a SLURM cluster, use the scratch directory!
* This version of the setup process can be run locally, so we just use "bash" to run it - no SLURM sbatch required

## And that is the installation and setup done!


# HOW TO RUN
* Within the bioinf/bin folder is a collection of scripts, written in bash, python, and R
* They should now be available through your PATH: in your terminal type bioinf and press TAB twice and you should see all of the scripts available
* All scripts send their output to the current directory from which the script was run
* Every script explains how it can be run if you give it the --help flag, e.g.:
```
bioinf-assembly-canu.sh -h
```

## RUNNING on a Computing Cluster, specifically SLURM
* Hopefully you are using a computing cluster, as many of these scripts use software that need a lot of computing resources
* There is no guidance for other systems other than local or SLURM - sorry, it is all that we know
* All scripts can be submitted on a SLURM cluster using sbatch - but you must provide the SLURM parameters
* These SLURM parameters are only examples - you must provide ones which are suited to your available resources and system configuration
* You SHOULD specify a number of threads to use through the SLURM flag --cpus-per-task, otherwise the scripts will use just 1 thread, and some will take a million years to finish
* Use the script flag --threads to specify thread numbers on non-SLURM systems

* In this example we have a bash script called script_for_submitting_slurm_jobs.sh, which contains the following code:
```
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=240GB
#SBATCH --partition ag2tb
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err

bioinf-binning.sh -i /path/to/samplelistfile -p project_name
```
* Then we submit the batch job script:
```
sbatch script_for_submitting_slurm_jobs.sh
```
* ALTERNATIVELY, you can sbatch without a script, and do it through the command line instead, which will produce exactly the same result as using the script e.g.:
```
sbatch --time=24:00:00 --cpus-per-task=12 --mem=240GB --partition long -o slurm.%N.%j.out -e slurm.%N.%j.err bioinf-binning.sh -i /path/to/samplelistfile -p project_name
```

### RUNNING PIPELINES WITHOUT A COMPUTING CLUSTER 

* You can also run any script locally, i.e. through your local machine or interactive terminal. Trying to run most of these scripts locally on a SLURM cluster will probably fail - use sbatch!
* Use the script flag --threads to specify thread numbers on non-SLURM systems
* Note that these pipelines are likely too resource intensive for most home computers
* running a script on a local machine:
```
bioinf-binning.sh -i /path/to/samplelistfile -p project_name --threads 24
```
* Please note that the bioinf-assembly-canu.sh script will NOT work on a non-SLURM system

# WHAT TO RUN
* In your terminal type bioinf and press TAB twice and you should see all of the scripts available
* For the most part, these scripts are modular and do not depend on you having run any other particular script. Exceptions to this will share common name structure - e.g. scripts beginning with "bioinf-binning" are part of a set workflow, i.e. you must run bioinf-binning.sh before running bioinf-binning-summarise.sh
## bioinf-assembly-canu.sh
## bioinf-binning.sh
## bioinf-binning-summarise.sh
## bioinf-phylogeny.py

## OTHER NOTES
### Viruses and Anvi'o
* Anvi'o is primarily made for bacterial metagenomics, so there are several features and associated outputs that will not be useful to virology. Do not be concerned if any default files or summaries from anvi'o contain little useful information
* Are your viruses not in nice, neat bins? Try this - in the anvi-interactive interface, under the "Main tab", click "Items Order", and select "Sequence Composition (D: Euclidean; L: Ward)"

### Customising the databases
* The database directories bioinf/db/CUSTOM_DMND and bioinf/db/CUSTOM_HMMS are customisable. If you place correctly formatted database files in these directories, then the pipeline will use them. You can even delete any databases built by the setup script
* Use soft links to databases files/folders that already exist on your system to avoid space wastage
* Please see the bioinf-setup script - it has a step that can build any custom DIAMOND database you want
* For custom HMMs you are unfortunately on your own, because anvi'o requires a very specific and annoying set of files to run a custom HMM as detailed [here](https://merenlab.org/2016/05/21/archaeal-single-copy-genes/). Alternatively, you can look at the files in bioinf/db/CUSTOM_HMMS created by the bioinf-setup script and copy them

#### Customising the kaiju databases
* Please see the bioinf-setup script - it has an extra step that allows you to build other kraken databases
* Note that is simply telling Kaiju to build one of the databases from it predetermined databses as detailed [here](https://github.com/bioinformatics-centre/kaiju)
* Any correctly formatted kaiju database in the bioinf/db directory will be used by the pipeline
* Use soft links to databases files/folders that already exist on your system to avoid space wastage
* Note that the setup script creates a kaiju database for viruses and bacteria/archaea/eukaryota

