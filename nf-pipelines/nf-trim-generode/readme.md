# nf-trim-generode: Historic & Modern DNA Mapping Pipeline

## Overview
This is a Nextflow DSL2 workflow designed for processing high-throughput sequencing data across two eras (historical and modern cohorts). It is an extension of nf-trim-merged-unmerged that incorporates GenErode features and takes both historical and modern reads. It does masking of repeats in the reference genome, adapter trimming, overlapping read merging, modern read trimming, mapping, and optional fastqc quality evaluation, AMBER evaluation, and mapdamage rescaling of historical bam files.

The pipeline maps historical samples with a user-chosen algorith to extract their average mapped read length distribution, uses this value to trim modern reads (limiting temporal length biases), and selects an alignment algorithm (`bwa aln` vs `bwa mem`) based on the historical average read length. It re-maps the historical reads if the mapper is different than initially done.

It is built to run on Old Dominion University's WAHAB cluster. 

---

## Installation & Environment Setup

The pipeline can be installed directly via git:

```bash
git clone https://github.com/philippinespire/nf-pipelines.git
cd ./nf-trim-merged-unmerged
```
This GitHub repository also contains other nextflow pipelines.


### Personal Conda Prerequisite
**Important for WAHAB HPC** For nextflow to correctly install conda environments, a personal conda installation is necessary. Installation instruction for miniconda (my personal favorite) and how to use it can be found [here](https://www.anaconda.com/docs/getting-started/miniconda/main)

## Folder structure
Before running, this is the setup:
```
nf-trim-merged-unmerged/  
├── data/  
│   ├── reference/         # Reference genome FASTA, .fai, .dict, and all BWA index layers  
│   └── symlinks/          # Symbolic links pointing to raw paired-end sequence files  
├── inputfiles/  
│   ├── samplesheet.csv    # Target metadata sheet (sample,era)  
│   └── fastq_filenames.txt # Legacy tracker single-column fallback sheet  
├── scripts/  
│   ├── split_reads.sh     # Custom length splitting script  
│   └── samremovedup.py    # Custom python script for duplicate removal  
├── main.nf  
├── mapping_modules.nf  
└── nextflow.config  
```

## Parameters

### Inputfiles

The pipeline requires the following inputfiles:

* Reference (fasta or fna format)
* Reference index files (.bwt, .ann, .sa, .pac, .ann, .amb )
* BED file with masked repeat regions
* File with the name and era of all fastq files to be processed (`inputfiles/samplesheet.csv`)
* Directory containing all fastq files

All inputfiles can be copied or generated from the GenErode directory. Below is code that can help to generate the inputfiles.

```bash
# Activate bash shell

bash

# Change directory to the nextflow pipeline
cd ./nf-trim-merged-unmerged

# Create new directories
mkdir data
mkdir ./data/reference
mkdir ./data/symlinks
mkdir inputfiles

# Create softlinks to the raw fastq files
ln -s /Generode/data/raw_reads_symlinks/modern/*fastq.gz ./data/symlinks

# Create fastq filenames file. Manually adjust the file if necessary (e.g. if not all samples from GenErode are to be used)
ls ./data/symlinks/*fastq.gz | xargs -n1 basename | cut -d "_" -f1,2,3 | uniq > ./inputfiles/fastq_filenames.txt

# Create softlinks to reference and repma bed file
# Adjust the path to the reference if necessary
ln -s /Generode/reference/<reference>.fasta ./data/reference/
ln -s /Generode/reference/<reference>.fasta.* ./data/reference/
ln -s /Generode/reference/<reference>.repma.bed ./data/reference/

```

**Important** The pipeline assumes that fastq files are named in the following way: 
`<sampleID>_<index>_<flowcellID>_R1.fastq.gz` eg. `TzoCMta031_1_22CVWFLT3_R1.fastq.gz`
Make sure that all fastq file names are unique and follow this structure. Otherwise the pipeline will fail.

### Sample CSV Structure Example 
This could be, for example, `inputfiles/samplesheet.csv`
```
sample,era
TzoCMta031_1_22CVWFLT3,historical
CviAPal001_Ex1_L4,modern
```

### Configuration

In the `main.nf` file, adjust the parameters to fit your project:

```bash
// --- Default Parameters ---
params.samplesheet = "${projectDir}/inputfiles/samplesheet.csv" // This is the preferred way to provide sample metadata, including sample IDs and eras (modern or historical)
params.indir        = "${projectDir}/data/symlinks" // This is the directory where the raw FASTQ files are expected to be located. The pipeline will look for files named <sample_id>_R1.fastq.gz and <sample_id>_R2.fastq.gz in this directory.
params.outdir       = "${projectDir}/results" // This is the directory where all output files will be written. The pipeline will create subdirectories for different types of output (e.g., fastq, bam, stats).
params.reference    = "${projectDir}/data/reference/<reference>.fasta" // This is the path to the reference genome FASTA file that will be used for read mapping. The <reference> placeholder should be replaced with the actual reference name (e.g., hg19, mm10).
params.bed_file     = "${projectDir}/data/reference/<reference>.repma.bed" // input bed file of repeats and CpG sites for downstream ANGSD analyses. Only used if params.run_repeatmasking is set to false.
mm10). This file is used for filtering reads during mapping and for calculating depth statistics.
params.historical_mapper        = "aln" // Default starting point for historical read mapping
params.run_repeatmasking        = true // Whether to run RepeatModeler and RepeatMasker on the reference genome (true/false). If false, it needs an input bed file of repeats and CpG sites for downstream ANGSD analyses.
params.run_historical_fastqc    = false
params.run_historical_mapdamage = false
params.use_historical_rescaled  = false
params.run_historical_amber     = false
params.run_modern_amber         = false
```
Note that running repeat masking is a fairly slow process (a few hours).

Other parameters that you are less likely to adjust can be found in the header of `main.nf`.

## Pipeline Features
1. Two metadata options: The configuration ingests sample information from one of two channels:
   1. Samplesheet (Preferred): Parses a comma-separated file (`inputfiles/samplesheet.csv`) containing sample identifiers and era tags (modern vs historical).
   2. Legacy Mode (Fallback): Scans a single-column flat file (`inputfiles/fastq_filenames.txt`) and sends all targets to the modern processing track.
2. Two-pass historical mapping
   1. First-pass mapping: Historical segments are cleaned via fastp and aligned using a baseline algorithm specified by the user (Default: `bwa aln`).
   2. Mapped length calculation: High-quality alignments are used to calculate average read length.
   3. Mapper choice: If the calculated read length is $\le$ 80 bp, the pipeline sets `bwa aln` (optimal for short, degraded molecules). If the calculated true mapped length is > 80 bp, the pipeline upgrades to `bwa mem`.
   4. Second-pass mapping: A second mapping pass is executed only if the chosen mapper differs from the first pass. 
3. Modern trimming: Long modern reads are split down the middle or truncated to match the average length of the historical set, reducing length-bias signatures in downstream analyses.


## To Run

It is helpful to run nextflow pipelines in tmux so that you can keep it running in the background even when logging off.

To create a new tmux screen:

```bash
tmux new -s nextflow
```

!Good to know! Leaving a tmux sessions is notoriously difficult. To leave and reopen, type the following commands:

```bash
# To leave, literally press these keys
Ctrl + B, followed by D

# To re-enter a session
tmux a -t nextflow

```

Then start the nextflow pipeline as follows (while in tmux):

```bash
# Load the nextflow module
module load container_env
module load nextflow

# Start the run
nextflow run main.nf -profile standard -resume
```

If everything went well, the pipeline will run and submit jobs to the queue. 
On the first run, a new conda environment will be created. This can take some time.

## Output
Results are organized as
```
results
├── data
│   ├── bam          # Merged and realigned bam and bam index files
│   ├── bam_rescaled # Rescaled historical bam and index files
│   ├── fastq        # Trimmed historical reads (R1 and R2), and trimmed (R1 and R2) and trimmed+merged modern reads
│   └── reference    # Repeat-masking files: .combined_mask.bed, .cleaned.regions, .chr, .regions, repma.angsd.txt/bin/idx (angsd sites files)
│       └── modeler  # Output from RepeatModeler
├── amber            # amber output
├── depth            # Average sequencing depth for the regions in the BED file, per sample
├── fastp            # fastp output
├── fastqc           # fastqc output
├── mapdamage        # Mapdamage output plots and txt files
└── stats            # historical_trimlength.txt (calculated ave historical read length)
```

## Software Stack
The pipeline uses:

* fastp (v0.23.4) — Quality-trimming and overlapping mate pairing/merging.
* seqtk (v1.4) — Down-sampling and structural end-read truncation mapping.
* bwa (v0.7.17) — Dual aln/mem short-read genome matching engine.
* samtools (v1.19) — Indexing, sorting, track filtering, and deep map coverage summaries.
* fastqc (v0.11.9) — Per-pass data quality assurance metrics.
* mapDamage (v2.2) — Postmortem historical damage assessment and base-quality score rescaling.
* AMBER (v2.0) — Target alignment quality extraction evaluations.
* angsd (v0.94) — Multi-individual genotype processing database assembly.
