# nf-paralog: Test for paralogs

## Overview
This is a Nextflow DSL2 workflow designed for processing high-throughput sequencing data to test for paralogs. It is an extension of nf-trim-generode that only maps modern reads, maps them with very few quality tests, runs an ANGSD test for HWE and depth outliers, and runs ngsParalog. It optionally runs the dupHMM step of ngsParalog.

It is built to run on Old Dominion University's WAHAB cluster. 

---

## Installation & Environment Setup

The pipeline can be installed directly via git:

```bash
git clone https://github.com/philippinespire/nf-pipelines.git
cd ./nf-paralog
```
This GitHub repository also contains other nextflow pipelines.


### Personal Conda Prerequisite
**Important for WAHAB HPC** For nextflow to correctly install conda environments, a personal conda installation is necessary. Installation instruction for miniconda (my personal favorite) and how to use it can be found [here](https://www.anaconda.com/docs/getting-started/miniconda/main)

If you run dupHMM (an optional step in this pipeline), you will need to install expm, truncnorm, and docopt packages in your home directory on Wahab, since the base R installation on R doesn't have them. The pipeline expects them specficially in `~/.R/library`:
```
bash
module load container_env R
mkdir -p ~/.R/library
crun Rscript -e '
  user_lib <- file.path(Sys.getenv("HOME"), ".R", "library")
  pkgs <- c("expm", "truncnorm", "docopt")
  install.packages(pkgs, lib = user_lib, repos = "https://cloud.r-project.org")
'
```

## Folder structure
Before running, this is the setup:
```
nf-paralog/  
├── data/  
│   ├── reference/         # Reference genome FASTA, .fai, .dict, and all BWA index layers  
│   └── symlinks/          # Symbolic links pointing to raw paired-end sequence files  
├── inputfiles/  
│   ├── samplesheet.csv    # Target metadata sheet (sample,era)  
├── scripts/  
│   └── samremovedup.py    # Custom python script for duplicate removal  
├── main.nf  
└── nextflow.config  
```

## Parameters

### Inputfiles

The pipeline requires the following inputfiles:

* Reference (fasta or fna format)
* Reference index files (.bwt, .ann, .sa, .pac, .ann, .amb )
* File with the name and era of all fastq files to be processed (`inputfiles/samplesheet.csv`)
* Directory containing all fastq files

All inputfiles can be copied or generated from the GenErode directory. Below is code that can help to generate the inputfiles.

```bash
# Activate bash shell

bash

# Change directory to the nextflow pipeline
cd ./nf-paralog

# Create new directories
mkdir data
mkdir ./data/reference
mkdir ./data/symlinks
mkdir inputfiles

# Create softlinks to the raw fastq files
ln -s /Generode/data/raw_reads_symlinks/modern/*fastq.gz ./data/symlinks

# Create softlinks to reference files
# Adjust the path to the reference if necessary
ln -s /Generode/reference/<reference>.fasta ./data/reference/
ln -s /Generode/reference/<reference>.fasta.* ./data/reference/

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
params.modern_era   = "modern"  // label in the "era" column of the samplesheet that identifies modern individuals
params.min_ind_ratio = 0.5      // Require coverage in at least 50% of samples
params.high_depth_quantile = 0.995 // target high depth percentile cutoff for flagging problematic high depth regions
params.lr_quantile  = 0.999     // Target percentile cutoff for ngsParalog log-likelihood values (e.g., 0.999 = top 0.1% highest LR sites)
params.rmdup_script = "${projectDir}/scripts/samremovedup.py"
params.ngsparalog_bin = "/archive/carpenterlab/pire/softwares/ngsParalog/ngsParalog"  // path to ngsParalog binary
params.run_duphmm   = true      // Toggle dupHMM execution
params.duphmm_script = "/archive/carpenterlab/pire/softwares/ngsParalog/dupHMM.R"     // Path to dupHMM.R
params.duphmm_emit  = 0         // 0 = LR only, 1 = both LR and Coverage
```

Other parameters that you are less likely to adjust can be found in the header of `main.nf`.

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
nextflow run main.nf -profile wahab -resume
```

If everything went well, the pipeline will run and submit jobs to the queue. 
On the first run, a new conda environment will be created. This can take some time.

## Output

### Directory Structure

```text
results/
├── angsd/
│   ├── high_depth_regions.bed
│   └── hwe_excess_het.bed
├── fastp/
│   ├── <sample>_fastp_report.html
│   └── <sample>_fastp_report.json
├── paralogs/
│   ├── ngsparalog_duphmm_regions.bed
│   ├── ngsparalog_sites.bed
│   └── ngsparalog_threshold.txt
└── large_data/
    ├── angsd/
    │   ├── cvi_angsd.arg
    │   ├── cvi_angsd.counts.gz
    │   ├── cvi_angsd.depthGlobal
    │   ├── cvi_angsd.depthSample
    │   ├── cvi_angsd.hwe.gz
    │   └── cvi_angsd.mafs.gz
    ├── bam/
    │   ├── <sample>.merged.realn.bam
    │   └── <sample>.merged.realn.bam.bai
    └── paralogs/
        ├── cvi_avg_depth.tsv
        └── cvi_ngsparalog.lr.txt

```

### Output File Descriptions

* **`angsd/hwe_excess_het.bed`**: A BED file listing genomic coordinates that significantly deviate from Hardy-Weinberg Equilibrium due to excess heterozygosity ($F < 0$ and $p < 10^{-3}$). Used downstream to mask suspected paralogous or duplicated regions.
* **`angsd/high_depth_regions.bed`**: BED file of contiguous genomic regions exceeding the empirical high-coverage quantile threshold (params.high_depth_quantile, default 99.5th percentile).
* **`fastp/`**: Quality control and read-preprocessing reports generated by `fastp` for each processed sample in HTML and JSON formats.
* **`paralogs/ngsparalog_sites.bed`**: A BED file containing candidate paralogous positions identified by `ngsParalog`, filtered using the dynamically calculated high-LR quantile threshold (`params.lr_quantile`).
* **`paralogs/ngsparalog_threshold.txt`**: Text file recording the calculated numeric LR threshold value derived from `params.lr_quantile` (default 99.9th percentile).
* **`paralogs/ngsparalog_duphmm_regions.bed`**: BED file containing genomic intervals classified as duplicated by dupHMM.R HMM state inference (only created when `params.run_duphmm = true`).

**`large_data/` (Raw Output Files)**

* **`large_data/bam/`**: Final processed alignments for modern samples. These BAMs include merged single-end and unmerged paired-end reads mapped without MAPQ filtering ($q=0$), marked for duplicates, merged by sample ID, and realigned around indels.
* **`large_data/angsd/`**: Unfiltered, raw calculation files output by ANGSD, including compressed Hardy-Weinberg test scores (`.hwe.gz`), major/minor allele frequencies (`.mafs.gz`), sample/global depth distributions (`.depthSample`, `.depthGlobal`), and base counts (`.counts.gz`).
* **`large_data/paralogs/cvi_ngsparalog.lr.txt`**: The unfiltered likelihood-ratio table aggregated across all contigs/scaffolds from `ngsParalog calcLR`.
* **`large_data/paralogs/cvi_avg_depth.tsv`**: Per-site average depth matching the `cvi_ngsparalog.lr.txt` site list (only generated when `params.duphmm_emit = 1`).

## Software Stack
The pipeline uses:

* Nextflow (DSL2): Pipeline execution, concurrency, and task management.
* `fastp`: Read cleaning, adapter/poly-G removal, quality filtering, and paired-end read merging.

* `bwa` (`bwa-mem`): Reference genome read mapping.
* `samtools`: Alignment sorting, indexing, merging, coordinate collation, markdup, faidx indexing, and mpileup generation.
* `GATK 3.x` (`GenomeAnalysisTK.jar`): Indel target identification (`RealignerTargetCreator`) and local realignment (`IndelRealigner`).
* `ANGSD`: Hardy-Weinberg Equilibrium (`-doHWE`), major/minor allele frequency estimation (`-doMaf`), depth calculations, and allele counting.
* `ngsParalog`: Likelihood-ratio calculations (`calcLR`) across mpileups to identify duplicated/paralogous regions.
* `R` (`v4.x+`) & `dupHMM.R`: Hidden Markov Model segmentation for paralog/duplication state inference using R packages `expm`, `truncnorm`, and `docopt`.
* Python 3: In-stream duplicate removal via `samremovedup.py`.
* Java Runtime Environment (JRE): Dependency for running GATK.
* POSIX Utilities (`awk`, `zcat`, `bash`): Genomic window splitting, BED conversion, filtering, and text handling.