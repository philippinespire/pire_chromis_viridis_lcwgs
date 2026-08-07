## Overview

This pipeline calculates diversity estimates within the ANGSD frameowork. It is written in Nextflow DSL2 and designed to be run on Old Dominion WAHAB cluster. The pipeline starts from BAM files.

**Important for WAHAB HPC:** For nextflow to correctly install conda environments, a personal conda installation is necessary. Installation instruction for miniconda (my personal favorite) and how to use it can be found [here](https://www.anaconda.com/docs/getting-started/miniconda/main)

## Folder structure

## Installation

The pipeline can be installed directly from git

```bash

git clone https://github.com/mariannedehasque/nf-pipelines.git

```
This GitHub repository also contains the `nf-mapping-ancient-merged` and `nf-trim-merged-unmerged` nextflow pipelines.

## Required files

The pipeline requires the following inputfiles:

* Reference (fasta or fna format)
* BED file with masked repeat regions
* BAM files

In addition, we need to create the following files too (see instructions below):
* Samplesheet
* ANGSD sites file
* Bam inputfile
* Contig file

### Samplesheet

This is a samplesheet unique to this pipeline. For an example, see ./inputfiles/samplesheet.csv

The samplesheet is a comma-separated file (csv) and consists of the following columns:

| sample     | bam                                               | pop     | era    | region    |
|------------|---------------------------------------------------|---------|--------|-----------|
| TzoCMal001 | /path/to/bam/TzoCMal001.merged.realn.bam          | TzoCMal | modern | Malampaya |
| TzoABol002 | /path/to/bam/TzoABol002.merged.realn.bam          | TzoABol | historic | Malampaya |
| TzoCMal003 | /path/to/bam/TzoCMal003.merged.realn.bam          | TzoCMal | modern | Malampaya |

Important:
* `era` takes one of two arguments: modern or historic
* `region` is used to define which historic and modern sites should be directly compared. Use the same value here if you want to directly compare two sites.
* `pop` is used to define the populations for which to calculate diversity matrixes

### ANGSD sites file
The ANGSD sites file contains information on the sites that will be analyzed. We will use this file to remove repeat regions from the analysis. This files corresponds to the `-sites` flag in [ANGSD](https://www.popgen.dk/angsd/index.php/Sites).

To create the sites file from the GenErode output (after loading all modules):

```bash
awk '{print $1"\t"$2+1"\t"$3}' ./path/to/reference/reference.repma.bed > ./path/to/reference/reference.repma.angsd.txt

angsd sites index ./path/to/reference/reference.repma.angsd.txt
```

### Bam inputfile
The filelist is a file containing the full path for each bam file with one filename per row.

This corresponds to the input given with the `-bam` flag in [ANGSD](https://www.popgen.dk/angsd/index.php/Input).

### Contig file
Specify the contigs/regions for which to run the pipeline. One contig/region per line. This corresponds to the `-r` flag in [ANGSD](https://www.popgen.dk/angsd/index.php/Input). The pipeline will submit a job per contig/region.

## Configuration

In the `main.nf` file, adjust the parameters. Make sure all files are in the correct directory.

In addition, the following parameters should be defined:
* maxdepth: I use 10X the expected total coverage here. The expected coverage can be calculated for the `dpstats.txt` files
* minind: I use 70-80% number of total individuals here, though this also depends on the dataset.

```bash
// Default Parameters
params.samplesheet = "${projectDir}/inputfiles/samplesheet.csv"
params.contigs     = "${projectDir}/data/reference/reference.contigs.txt"
params.outdir      = "${projectDir}/results"
params.reference   = "${projectDir}/data/reference/reference.fasta"
params.bed_file    = "${projectDir}/data/reference/reference.repma.angsd.txt"
params.species     = "Sor"
params.maxdepth    = 1000
params.minind      = 20

```

## To Run

I like to run nextflow pipelines in tmux so that I can keep it running in the background even when logging of.

To create a new tmux screen:

```bash
tmux new -s nextflow
```

!Good to know! Leaving a tmux sessions is notoriously difficult. To leave and reopen, type the following commands:

```bash
# To leave, literraly press these keys
Ctrl + B, followed by D

# To re-enter a session
tmux a -t nextflow

```

Then start the nextflow pipeline as follows (while in tmux):

```bash
# Load the nextflow module
module load container_env
module load nextflow

# Activate bash
bash

# Start the run
nextflow run main.nf -profile standard -resume

```

If everything went well, the pipeline will run and submit jobs to the queue. 
On the first run, a new conda environment will be created. This can take some time.
