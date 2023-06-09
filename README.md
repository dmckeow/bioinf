# bioinf
* A basic pipeline composed of custom scripts that run various publicly available tools.
* Main purpose is for common bioinformatics protocols including assembly and binning of sequence data, phylogeny, and other genomic analyses.
* Currently the pipeline is focused on using Oxford Nanopore Sequencing data, and on viral metagenomics

## INSTALLATION STEP 0 - get CONDA and MAMBA
* Before you begin, you need [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html). Even if conda is already installed through your server system, **it will probably be better if you install your own copy**
* In addition to conda, you should also install [mamba](https://mamba.readthedocs.io/en/latest/installation.html). This allows you you to use `mamba install` instead of `conda install` which **much faster**

## INSTALLATION STEP 1 - get anvi'o
* Anvi'o is the tool we will use for binning - the developers are the Merenlab, and they have very good guides for [anvi'o](https://anvio.org/).
* First thing is to create a conda environment with Anvi'o 7.1 installed within it
* Please follow the exact instructions for [installing anvio](https://anvio.org/install/), BUT with one **IMPORTANT EXCEPTION** - name your environment `bioinftools` instead of `anvio-7.1`. The pipeline will look for a conda environment named bioinftools:
```shell
## INSTEAD of:
conda create -y --name anvio-7.1 python=3.6
 ## do:
conda create -y --name bioinftools python=3.6
```
* `The anvio page uses conda install, but you should use mamba install instead`
* If you have a problem with installation, the answer is probably on the anvi'o [installation guide](https://anvio.org/install/)
* IMPORTANT - during your installation, be sure to follow their step that details how to setup key resources for anvi'o - currently this is under: (4.1) Setup [key resources](https://anvio.org/install/)

* Here is a copy of the anvi'o packages listed on the anvi'o installation page, but with `mamba install` and with specific verions for all packages. This for convenience sake, and also in case any package conflicts occur due to the wrong package version being installed
```shell
mamba install -y -c bioconda "sqlite>=3.31.1"
mamba install -y -c bioconda prodigal=2.6.3
mamba install -y -c bioconda mcl=14.137
mamba install -y -c bioconda muscle=3.8.1551
mamba install -y -c bioconda hmmer=3.3.2
mamba install -y -c bioconda diamond=2.1.6
mamba install -y -c bioconda blast=2.5.0
mamba install -y -c bioconda megahit=1.2.9
mamba install -y -c bioconda spades=3.15.5
mamba install -y -c bioconda bowtie2 tbb=2019.8
mamba install -y -c bioconda bwa=0.7.17
mamba install -y -c bioconda samtools=1.9
mamba install -y -c bioconda centrifuge=1.0.4_beta
mamba install -y -c bioconda trimal=1.4.1
mamba install -y -c bioconda iqtree=2.2.2.3
mamba install -y -c bioconda trnascan-se=2.0.9
mamba install -y -c bioconda r-base=4.2.2
mamba install -y -c bioconda r-stringi=1.7.12
mamba install -y -c bioconda r-tidyverse=2.0.0
mamba install -y -c bioconda r-magrittr=2.0.3
mamba install -y -c bioconda r-optparse=1.7.3
mamba install -y -c bioconda bioconductor-qvalue=2.30.0
mamba install -y -c bioconda fasttree=2.1.11
mamba install -y -c bioconda vmatch=2.3.0

# this last one may cause some issues. if it doesn't install,
# don't worry, you will still be fine:
mamba install -y -c bioconda fastani=1.33
```

## INSTALLATION STEP 2 - get other tools
* Copy the following into your terminal:
```shell
conda activate bioinftools ## your conda environment must be active to use it

mamba install -y -c bioconda kaiju=1.9.2
mamba install -y -c bioconda seqkit=2.4
mamba install -y -c bioconda hyphy=2.5
mamba install -y -c bioconda minimap2=2.2
mamba install -y -c bioconda porechop=0.2
mamba install -y -c conda-forge pigz
mamba install -y -c bioconda taxonkit=0.14
mamba install -y -c bioconda mafft=7.5

### packages for isONcorrect:
pip install isONcorrect
mamba install -c bioconda spoa
pip install isONclust
mamba install -c bioconda "pychopper>=2.0"
```

* Sadly, installing canu through conda is a bit dodgy - so install version 2.2 of [CANU](https://github.com/marbl/canu ) by `downloading a binary release as detailed in the CANU git`
* Once it is installed you must add canu to your path as follows:
```shell
sed -i -z "s|$|\nexport PATH=\"path\/to\/your\/canu-2.2\/bin:\$PATH\"\n|g" $HOME/.bashrc ## replace path\/to\/your\/ with your actual path - you MUST use \/ instead of / in your path
source $HOME/.bashrc
canu --help ## if you see canu's help message then you are good to go
```
* Alternatively, you can add canu to your path by just copy pasting `export PATH="path/to/your/canu-2.2/bin:$PATH"` into your ~/.bashrc using an editor like vim or nano, save the file changes and then do `source ~/.bashrc`
* If this absolutely isn't working you can also try to install canu using conda

## INSTALLATION STEP 3 - download and setup bioinf pipelines, scripts, and databases
`CHOOSE BETWEEN 3a. 3b. or 3c. :`
### 3a. For a fresh setup
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

### 3b. IF your computing group has a pre-existing bioinfdb that you can share:
* If someone has already setup bioinf by running the full bioinf-setup.sh script (and you have access permissions to it), then consider doing this instead - it will simply configure your version of bioinf to run using the pre-existing setup. This will complete within minutes and save 100s in GB of storage space.
* Creates symbolic links to pre-existing databases. Each user still has their own database directory. Each user can then add/remove databases without affecting anyone elses database setup (EXCEPT that deleting files in the original database will remove them from all linked db copies)
* Copy the following code, but with your own paths:
```shell
git clone https://github.com/dmckeow/bioinf.git
chmod +x bioinf/bin/*
bioinf/bin/bioinf-setup.sh -t /path/to/temporary_storage -d /where/to/create/your/database/directory -E /path/to/existing/bioinfdb ## replace with your actual paths - if bioinf is already setup (step A1 ahs put it on your PATH), then you can do bioinf-setup.sh instead of bioinf/bin/bioinf-setup.sh
source ~/.bashrc ## OR logout and log back in
```
* Temporary storage is specified because some scripts such as bioinf-assembly-canu.sh will generate quite large temporary files - so you should choose somewhere with enough capacity to handle 10s to 100s of GB of data. If you are on a SLURM cluster, use the scratch directory!
* This version of the setup process can be run locally - no SLURM sbatch required

### 3c. IF you already have your own personal bioinfdb setup and just need to update the locations of the databases and script files
* This is useful if you need to change the location of your bioinf, bioinfdb, or bioinftmp
* Intended for a user to re-connect to their own database
* BEWARE: if you leave out -s A1, then a fresh setup will begin and your existing db will get deleted
```shell
git clone https://github.com/dmckeow/bioinf.git
chmod +x bioinf/bin/*
bioinf/bin/bioinf-setup.sh -t /path/to/temporary_storage -d /path/to/directory/containing/your/bioinfdb -s A1 ## replace with your actual paths
source ~/.bashrc ## OR logout and log back in
```


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
    |------DMND ## .dmnd database files built by diamond makedb (see below)
    |
    |------HMM ## subdirectories each containing a single hmmer profile built from a set of reference proteins (see below)
    |
    |------KAIJU ## subdirectories each containing a single kaiju database built by kaiju-makedb (see below)
    |
    |------VOGDB ## this is simply a copy of VOGDB. There is no reason to customise this

```
* Use soft links to databases files that already exist on your system to avoid space wastage - e.g. `ln -s /path/to/shared/nr.dmnd /path/to/my/bioinfdb/DMND/nr.dmnd`

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