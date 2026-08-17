## Overview of nf-angsd-selection

This pipeline runs (optional) chi-squared selection scans from data across two time-points, does (optional) ld-pruning, calculates diversity estimates within the ANGSD framework on all putatively neutral sites (including monomorphic sites), and then makes PCA and admixture plots and calculates Fst with ANGSD on putatively neutral ld-pruned sites. It only considers sites that meet minimum individual and maximum depth limits in each of the population-eras (set as a fraction and a multiplier, respectively).

It is written in Nextflow DSL2 and designed to be run on Old Dominion WAHAB cluster (use `-profile wahab`). The pipeline starts from BAM files, such as those output by nf-trim-generode. It was developed from nf-angsd-diversity by Malin Pinsky (August 2026) with plenty of help from Gemini.

**Important for WAHAB HPC:** For nextflow to correctly install conda environments, a personal conda installation is necessary. Installation instruction for miniconda (my personal favorite) and how to use it can be found [here](https://www.anaconda.com/docs/getting-started/miniconda/main)

## Installation

The pipeline can be installed along with the others in this repo directly from git

```bash

git clone https://github.com/mariannedehasque/nf-pipelines.git

```

## Required files

The pipeline requires the following inputfiles:

* Reference (fasta or fna format)
* BED file with masked repeat regions (eg, output by Generode or by nf-trim-generode)
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
* `region` is used to define which historic and modern sites should be directly compared. Use the same value here if you want to directly compare sites across eras. Diversity is calculated for each region-era combination.
* `pop`: not used.

If you want to use mapdamage-rescaled BAM files for historical samples (eg, from Generode or nf-trim-generode), point to those files in this input file.


### ANGSD sites file
The ANGSD sites file contains information on the sites that will be analyzed. We will use this file to remove repeat regions from the analysis. This files corresponds to the `-sites` flag in [ANGSD](https://www.popgen.dk/angsd/index.php/Sites).

To create the sites file from the GenErode output (after loading all modules):

```bash
awk '{print $1"\t"$2+1"\t"$3}' ./path/to/reference/reference.repma.bed > ./path/to/reference/reference.repma.angsd.txt

angsd sites index ./path/to/reference/reference.repma.angsd.txt
```

Alternatively, this file is created automatically as part of the nf-trim-generode pipeline and output in `results/data/reference/*.repma.angsd.txt`.

### Bam inputfile
The filelist is a file containing the full path for each bam file with one filename per row.

This corresponds to the input given with the `-bam` flag in [ANGSD](https://www.popgen.dk/angsd/index.php/Input).

If you want to use mapdamage-rescaled BAM files for historical samples (eg, from Generode or nf-trim-generode), point to those files in this input file.

### Contig file
Specify the contigs/regions for which to run the pipeline. One contig/region per line. This corresponds to the `-r` flag in [ANGSD](https://www.popgen.dk/angsd/index.php/Input). The pipeline will submit a job per contig/region.

## Configuration

In the `main.nf` file, adjust the file paths for your project's reference and bed file.

In addition, the following parameters should be defined:
* `min_ind_frac`: The fraction of individuals that have to be represented at a locus/site to use it. This is applied per region-era and the intersection of the passing loci in each region-era is used (default 0.7)
* `max_depth_mult`: The multiple of the average depth to use as a ceiling in order to include a locus/site. This filter is applied per region-era and the intersection of the passing loci in each region-era is used (default 10)
* `ld_prune`: Whether to prune out loci in linkage disequilibrium (often a good idea). The PCA, admixture, and FST calculations will then use the LD-pruned set of loci (default true).
* `max_kb_dist`: LD pruning uses a sliding window to test for disequilibrium. This sets the window width. Narrower will run faster, wider will run more slowly (default 50)
* `min_weight`: Loci that are correlated (r2) more than this threshold will be trimmed out (default 0.2)
* `run_selection`: Whether to run chi-squared selection scan (default true)
* `generations`: Number of generations to use for Ne calculations (default 114 for Albatross to contemporary populations with a 1 year generation time)
* `fdr_cutoff`: False detection rate (FDR) cutoff for iterative selection scan (default 0.05)
* `max_rounds`: Maximum number of rounds to run for iterative selection scan (default 20)
* `n_boot`: Number of bootstrap replicates to run for Ne calculations (default 1000)
* `min_ind`: Minimum number of individuals to include in Ne calculations (default 4)
* `fst_window`: Window size for sliding window FST calculations (default 50000)
* `fst_step`: Step size for sliding window FST calculations (default 10000)


```bash
// Default Parameters
params.samplesheet = "${projectDir}/inputfiles/samplesheet.csv" // sample sheet csv with columns for sample name, bam path, pop, era, and region. era is historic or modern. region matches historic and modern populations for diversity calculations. pop is not used
params.contigs     = "${projectDir}/inputfiles/contig_list.txt" // text file with one contig name per line. Contig names must match those in the reference genome FASTA file.
params.outdir      = "${projectDir}/results" // output directory for all results
params.reference   = "${projectDir}/data/reference/reference.fasta"
params.bed_file    = "${projectDir}/data/reference/reference.repma.angsd.txt"
params.species     = "Cvi" // 3-letter species code for output file naming. used in ANGSD output files and ACER results.
params.min_ind_frac   = 0.70  // Minimum percentage of individuals required per site (default 0.70)
params.max_depth_mult = 10    // Max depth threshold as a multiplier of the average total expected coverage (default 10)
params.ld_prune    = true   // Set to true to enable LD pruning (default true)
params.max_kb_dist = 50     // Maximum pairwise distance in kb to test for LD if pruning
params.min_weight  = 0.2    // Minimum r2 threshold for pruning filter (default 0.2)
params.run_selection = true  // Set to false to skip ACER by default (default true)
params.generations = 114    // Number of generations to use for Ne calculations. Default is 114 for Albatross to contemporary populations with a 1 year generation time. Adjust for other species as needed.
params.fdr_cutoff  = 0.05   // FDR cutoff for iterative selection scan (default 0.05)
params.max_rounds  = 20     // Maximum number of rounds to run for iterative selection scan (default 20)
params.n_boot      = 1000   // Number of bootstrap replicates to run for Ne calculations (default 1000)
params.min_ind     = 4      // Minimum number of individuals to include in Ne calculations (default 4)
params.fst_window = 50000  // Window size in bp for sliding window FST calculations (default 50000)
params.fst_step   = 10000  // Step size in bp for sliding window FST calculations (default 10000)


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
nextflow run main.nf -profile wahab -resume

```

If everything went well, the pipeline will run and submit jobs to the queue. 
On the first run, a new conda environment will be created. This can take some time.

## Output
## Pipeline Output Directory Structure

All pipeline outputs are written to the directory specified by `--outdir` (default: `./results/`). Large files are isolated in the `large_data/` directory so they can be easily managed or excluded via `.gitignore`.

results/
├── inputfiles/                # Input BAM file lists
├── sites/                     # Summary metrics across generated site set libraries
├── selection/                 # ACER iterative selection scan outputs
├── PCAngsd/                   # Population structure (PCA & Admixture)
├── diversity/                 # Summary nucleotide diversity statistics and comparison plots
├── fst/                       # Global FST values and Manhattan plots
└── large_data/                # Large binary, intermediate, and site library files (Ignored by git)
    ├── sites/                 # Base coordinate indices, region mappings, and processed site libraries
    ├── beagle/                # Individual genotype likelihood matrices
    ├── mafs/                  # Point allele frequency tables
    ├── saf/                   # Population sample allele frequency distributions
    ├── thetas/                # Per-site nucleotide diversity estimates
    └── fst/                   # Binary FST indices, 2D SFS matrices, and windowed FST tables


### Output Directories Tracked by Git

#### 1. `inputfiles/`
Contains formatted list files of BAM paths passed to ANGSD processes
* **`bamlist.txt`**: Plain text list of all BAM files used across the full population dataset.
* **`${pop}.${era}.bamlist.txt`**: BAM file paths partitioned by specific population/region and sampling era.

#### 2. `sites/`
Contains site library summary metrics.
* **`site_counts.tsv`**: Summary table detailing total site counts across all processed site subsets (`all_callable`, `all_snps`, `selected_loci`, `callable_neutral`, `snps_neutral`, `snps_pruned`, `snps_pruned_neutral`) generated by `PREPARE_SITE_SETS`.

#### 3. `selection/`
Contains output tables and plots from the ACER iterative selection scan.
* **`iteration_summary.tsv`**: Round-by-round summary of effective population size ($N_e$) estimation and selection iterations.
* **`final_test_results.tsv`**: Per-SNP selection scan statistics, including raw $\chi^2$ test $p$-values and FDR-adjusted values.
* **`ne_bootstrap.tsv`**: Bootstrap replicates used to estimate baseline $N_e$.
* **`*.png`**: Manhattan plots and diagnostic visualizations generated during the selection scan.

#### 4. `PCAngsd/`
Contains population structure and individual admixture results derived from LD-pruned neutral loci.
* **`${params.species}.cov`**: Estimated individual covariance matrix for PCA calculation.
* **`${params.species}.*.Q`**: Individual ancestry proportion matrix (Q-matrix) generated by PCAngsd admixture analysis.
* **`${params.species}.pcangsd.plot.pdf`**: Publication-ready PCA scatter plot colored by population and era.
* **`${params.species}.admixture.pdf`**: Stacked bar plot displaying individual admixture proportions.

#### 5. `diversity/`
Summary nucleotide diversity files calculated across all callable neutral loci (monomorphic + polymorphic).
* **`${region}_${era}.sfs`**: Folded site frequency spectrum (SFS) calculated across callable neutral sites.
* **`${region}_${era}.pestPG`**: Windowed summary statistics generated by `thetaStat`, including pairwise nucleotide diversity ($\pi$), Watterson's $\theta$, and Tajima's $D$.
* **`pi_summary_${region}.tsv`**: Summary table of mean nucleotide diversity ($\pi / \theta_P$) with 95% bootstrap confidence intervals across eras.
* **`pi_historic_vs_modern_${region}.png`**: Comparison plot of historic vs. modern nucleotide diversity per region.

#### 6. `fst/`
$F_{ST}$ metrics between historical and modern populations derived from LD-pruned neutral loci.
* **`${region}_hist_vs_mod.global_fst.txt`**: Global unweighted and weighted $F_{ST}$ values, in that order.
* **`${region}_fst_manhattan.png`**: Manhattan plot displaying genome-wide $F_{ST}$ sliding window calculations across chromosomes.

---

### `large_data/` Directory (Untracked / Ignored by Git)

Contains heavy intermediate data, binary index files, and site set libraries (excluded by adding `results/large_data/` to `.gitignore`).

#### 1. `large_data/sites/`
Base genomic coordinate indices, region mappings, and processed site set libraries.
* **`all_sites.pos`** (`.bin`, `.idx`, `.regions`): Reference BED genomic targets (all 1-based genomic positions extracted directly from the input BED file prior to sequence read coverage or quality filtering) generated by `INDEX_BED_SITES`.
* **`${params.species}.sites.txt`** (`.bin`, `.idx`, `.regions`, `.chrs`): Polymorphic site indices, region mappings, and chromosome lists output by `ANGSD_EXTRACT_SITES`. Not filtered by minind or maxdepth.
* **`all_callable.pos`**: QC-filtered shared callable sites (the intersected subset of reference BED positions passing `minind_frac` and `maxdepth_mult` thresholds across every region-era subset) extracted by `PREPARE_SITE_SETS`. Includes monomorphic sites.
* **`all_snps.pos`**: Full set of polymorphic SNPs passing quality and $p$-value filters.
* **`selected_loci.pos`**: Position coordinates of candidate loci identified under selection by ACER.
* **`callable_neutral.pos`** (`.bin`, `.idx`): All callable positions minus candidate selected loci.
* **`snps_neutral.pos`**: Polymorphic SNPs minus candidate selected loci.
* **`snps_pruned.pos`**: LD-pruned polymorphic SNP coordinates.
* **`snps_pruned_neutral.pos`** (`.bin`, `.idx`): Polymorphic SNPs that are both LD-pruned and neutral.

#### 2. `large_data/beagle/`
Full-genome individual genotype likelihood matrices in Beagle format.
* **`${params.species}.beagle.gz`**: Concatenated individual genotype likelihood matrix generated by `ANGSD_COLLECT_OUTPUT`. The positions are those in `all_snps.pos`.

#### 3. `large_data/mafs/`
Point allele frequency tables across populations and subsets.
* **`${params.species}.mafs.gz`**: Combined SNPS frequency table across all individuals. The positions are those in `all_snps.pos`.
* **`${region}_${era}.mafs.gz`**: Major/minor allele frequency outputs per population/era subset, including monomorphic. The positions are those in `all_callable.pos`.

#### 4. `large_data/saf/`
Population sample allele frequency distributions used for downstream site frequency spectra.
* **`${region}_${era}.saf.gz`** (`.idx`, `.pos.gz`): Binary Sample Allele Frequency distributions generated by `ANGSD_GL_POP`. The positions are those in `all_callable.pos`.

#### 5. `large_data/thetas/`
Per-site binary likelihoods for nucleotide diversity estimates.
* **`${region}_${era}.thetas.gz`** (`.idx`): Binary per-site theta likelihood output generated by `realSFS saf2theta`. The positions are those in `all_callable.pos`.

#### 6. `large_data/fst/`
Binary $F_{ST}$ indices, 2D-SFS matrices, and sliding window tables.
* **`${region}_hist_vs_mod.fst.gz`** (`.idx`): Binary per-site $F_{ST}$ index files generated by `realSFS fst index`. The positions are those in `snps_pruned_neutral.pos`.
* **`${region}_hist_vs_mod.2dsfs`**: 2D Site Frequency Spectrum matrix estimated between historical and modern populations. The positions are those in `snps_pruned_neutral.pos`.
* **`${region}_windowed_fst.txt`**: Genome-wide sliding window $F_{ST}$ table output by `realSFS fst stats2`. The positions are those in `snps_pruned_neutral.pos`.
