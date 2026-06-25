# Directory Structure

| directory-name | date-created | dataset |
| -------------- | ------------ | ------- |
| 1st_sequencing_run | 2023-04-28 | test lane |
| 2nd_sequencing_run | 2023-06-16 | Alb & Contemp sequences from Palapag, Northern Samar|
| angsd_analysis | 2023? | Jem's angsd analysis |
| scripts | 2026-05 | scripts for analysis |
| output | 2026-05 | output from analyses |
| nf-pipelines | 2026-05 | nextflow pipelines for trimming & angsd |
| temp | 2026-05 | files not tracked by git |

# Chromis viridis lcWGS

---
Jem Baldisimo and Malin Pinsky
---
This repository outlines the roadmap we followed to move *Chromis viridis* through the [Low Coverage Whole Genome Sequencing Pipeline](https://github.com/philippinespire/pire_lcwgs_data_processing). This provides steps taken & analysis we did to gain insight on the historical population demography of this very popular aquarium fish, also known as the green chromis.

---

## 1. Completed fq.gz pre-processing
Working in `/archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/`

The [pire_fq_gz_processing](https://github.com/philippinespire/pire_fq_gz_processing) instructions and scripts (Garcia et al 2021) were used. The purpose of the preprocessing pipeline was to Trim, deduplicate, decontaminate, and repair the raw `fq.gz` files.

<details><summary><i>Checked quality of data using Multi_FASTQC</i></summary>
<p>


Multi_FASTQC [Report](https://github.com/philippinespire/pire_lcwgs_data_processing/blob/main/chromis_viridis/2nd_sequencing_run/fq_raw/fqc_raw_report.html)
```
Potential issues:
 * duplication - moderate to high in Albatross; low to modeterate in Contemporary --- Higher in Albatross than in  Contemporary
   * Alb: 22.8-79.4%, Contemp: 4.3-32%%
 * gc content - reasonable
   *41%-57%
 *number of reads
   *Albatross - 3.7-131.2 M seq; Contemporary - 1.8-56 M seq
```
</p>
</details>

<details><summary><i>1st trim</i></i></summary>
<p>

 1st FASTP [Report](https://github.com/philippinespire/pire_lcwgs_data_processing/blob/main/chromis_viridis/2nd_sequencing_run/fq_fp1/1st_fastp_report.html)
```
Potential issues:
 * duplication - Low to moderate in Albatross, low in Contemp
   * Alb: 14.4 to 61.8%, Contemp: 2.4-16%
 * gc content - Reasonable
   * Alb: 39.6-56.8%, Contemp: 41-44.7%
 * Passing filter
   * Alb: 89.8-98.4%, Contemp: 75.4-98.6%
 * % adapter 
   * Alb: 7-86.6% (Only 2 individuals had >77%), Contemp: 9.7-59.7% (Only 5 had >55%)%
 * number of reads
   * Albatross - 3.8-252 M, at least 10 M for majority; Contemp - 3-110 M, at least 7M reads for Contemporary
```
</p>
</details>

<details><summary><i>Clumpify</i></summary>
<p>

 Clumpify [Report](https://github.com/philippinespire/pire_lcwgs_data_processing/blob/main/chromis_viridis/2nd_sequencing_run/fq_fp1_clmp/fqc_clmp_report.html)

Clumpify worked well so I moved on to the next step
</p>
</details>

<details><summary><i>2nd Trim</i></summary>
<p>
Fastp2 [Report](https://github.com/philippinespire/pire_lcwgs_data_processing/blob/main/chromis_viridis/2nd_sequencing_run/fq_fp1_clmp_fp2/2nd_fastp_report.html)

```
Potential issues:
  * % duplication -low to moderate in Alb; low in Comtemp
    * Alb: 2.3-36.0%, Contemp: 0.4-4.6%
  * GC content - reasonable
    * Alb: 39.6-56.8%, Contemp: 41-44.7%
  * passing filter - great
    * Alb: 98.7-99.6%, Contemp: 99.4-99.6%
  * % adapter -
    * Alb: 0.2-1.1%, Contemp: 0.1-1.1%
  * number of reads -
    * Alb: 5-158 M, mostly at least 15 M, Contemp: 3-94 M, mostly at least 5M
```
</p>
</details>

<details><summary><i>Checked Fastqscreen files</i></summary>
<p>

After checking whether runFQSCRN_6.bash worked, I found that there was an error:

```
slurm-fqscrn.1735299.2.out:No reads in Cvi-APal_002-Ex1-1B-lcwgs-1-2.clmp.fp2_r1.fq.gz, skipping
slurm-fqscrn.1735299.5.out:No reads in Cvi-APal_003-Ex1-1C-lcwgs-1-2.clmp.fp2_r2.fq.gz, skipping
```

So I ran FQSCRN again for the 2 individual sequences:
```
bash /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 1 Cvi-APal_002-Ex1-1B-lcwgs-1-2.clmp.fp2_r1.fq.gz
bash /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 1 Cvi-APal_003-Ex1-1C-lcwgs-1-2.clmp.fp2_r2.fq.gz
```

There were no more errors or failed individuals so I ran Multiqc for all samples again

 MultiQC [Report](https://github.com/philippinespire/pire_lcwgs_data_processing/blob/main/chromis_viridis/2nd_sequencing_run/fq_fp1_clmp_fp2_fqscrn/fastq_screen_report.html)

```
Potential issues:
  * one hit, one genome, no ID -
     Alb: 66-92%, Contemp: 89-97%
  * no one hit, one genome to any potential contaminators (bacteria, virus, human, etc) -
    * Alb: 8-44%, Contemp: 3-11%% 
```
</p>
</details>


<details><summary><i>Reviewed Fastqscreen repaired files</i></summary>
<p>
MultiQC [Report](https://github.com/philippinespire/pire_lcwgs_data_processing/blob/main/chromis_viridis/2nd_sequencing_run/fq_fp1_clmp_fp2_fqscrn_rprd/fqc_rprd_report.html)

 ```
Potential issues:
  * % duplication - wide range of duplication for Albatross
     Alb: 3-34.6%, Contemp: 1.2-8.3%%
  * GC content -
    * Alb: 39-53%, Contemp: 40-43%
  * number of reads -
    * Alb: 2-73 M, Contemp: 1.4-45 M
```
</p>
</details>

---

## 2. Copied the Reference Genome into the species folder
Working in `/archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/`

I copied the best assembly from the ssl pipeline (see link below) and renamed it to reference.ssl.Cvi.fasta that was copied into the refGenome folder
```
/home/e1garcia/shotgun_PIRE/pire_ssl_data_processing/chromis_viridis/SPAdes_Cvi-CPal-A_decontam_R1R2_noIsolate/scaffolds.fasta
```

---
## 3. Mapped & Filtered BAM files
_June 2026 note: this step is no longer used. Replaced by [Step 7](#7-nextflow-trimming) that trims and maps reads_.  

Working in `/archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/`

<details><summary><i>mkBAM step</i></summary>
Created mkBAM folder & linked fq.gz files from fq_fp1_clmp_fp2_fqscrn_repaired:

```
cd /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/2nd_sequencing_run
mkdir mkBAM
ln fq_fp1_clmp_fp2_fqscrn_rprd/*fq.gz mkBAM
```

Copied dDocent config file & edited it to suit the lcwgs data
```
cp /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/sphaeramia_nematoptera/1st_sequencing_run/mkBAM/config.6.lcwgs .
```

After changing the settings, I executed this command:
```
sbatch /home/e1garcia/shotgun_PIRE/dDocentHPC/dDocentHPC2.sbatch mkBAM config.6.lcwgs
```

</p>
</details>

<details><summary><i>fltrBAM step</i></summary>

To filter resulting BAM files, I executed the following command:
```
sbatch /home/e1garcia/shotgun_PIRE/dDocentHPC/dDocentHPC2.sbatch fltrBAM config.6.lcwgs
```

</p>
</details>

## 4. Generated Mapping stats using mappedReadStats
_June 2026 note: this step is no longer used. Replaced by [Step 7](#7-nextflow-trimming) that trims and maps reads_.   

Working in `/archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/`

```
#navigate to mkBAM folder
# mappedReadStats.sbatch <Path to BAM file dir> <coverageMappedReads>
sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/mappedReadStats.sbatch . coverageMappedReads
```

---
## 5. Visualised results using the Process Sequencing Metadata Repo
Working in `/archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/`

I followed Kevin's repo for [Process Sequencing Metadata](https://github.com/philippinespire/process_sequencing_metadata).

Make sure  you do the following:  
-If cloning a repo and working in your own folder, make sure to always git pull before doing anything!
-Move output files to the species repo after you work in your own folder
-Remember to change line 306 of the visualizeSslCsslLcwgsMETADA.R script. For lcwgs, I changed this line based on what my mkBAM_T setting was and the 2nd to last setting for fltrmBAM:
```
  mutate(mkBAM_T_setting = ceiling(25 * fqc_avg_sequence_length_rprd / 100))
```

Output files are under the process_sequencing_metadata folder

Graphs showed information on depth of coverage for both Albatross & Contemporary

---
## 6. GenErode
_June 2026 note: this step is no longer used. Replaced by [Step 7](#7-nextflow-trimming) that trims and maps reads_.   
Working in `/archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/`

Jem ran this, and then ANGSD. No details available right now.
```
find /archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/GenErode/modern -maxdepth 1 -type f -name 'Cvi-CPal_*' -printf '%f\n' | cut -c 10-12 | sort | uniq | wc -l

ls /archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/GenErode/results/modern/mapping/reference.ssl.Cvi20k_rename/*.merged.rmdup.merged.realn.bam | wc -l

ls /archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/GenErode/results/modern/mapping/reference.ssl.Cvi20k_rename/*.merged.rmdup.merged.realn.bai | wc -l

find /archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/GenErode/historical -maxdepth 1 -type f -name 'Cvi-APal_*' -printf '%f\n' | cut -c 10-12 | sort | uniq | wc -l

ls /archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/GenErode/results/historical/mapping/reference.ssl.Cvi20k_rename/*.merged.rmdup.merged.realn.rescaled.bam | wc -l

ls /archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/GenErode/results/historical/mapping/reference.ssl.Cvi20k_rename/*.merged.rmdup.merged.realn.rescaled.bam.bai | wc -l

ls /archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/GenErode/results/gerp/reference.ssl.Cvi20k_rename.ancestral.rates.gz | wc -l
```

## 7. NextFlow Trimming
Malin Pinsky 2026 May. Working in `/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/`

This step applied Marianne's Nextflow trimming script to _all_ individuals and mapped them with bwa mem. This is not a standard application, since typically only the modern individuals are trimmed and they are mapped with bwa aln. By doing this, I didn't rescale the historical reads based on damage patterns, which is done in the generode pipeline. This is probably ok, since the reads have very little damage.

Cloned the nf-piplines repo and removed its status as a git repo (removed .git/ and .gitignore).
```
git clone https://github.com/mariannedehasque/nf-pipelines.git
rm -rf nf-pipelines/.git
rm nf-pipelines/.gitignore
```
From within `nf-trim-merged-unmerged/`, Made directories per instructions in [nf-trim-merged-unmerged](https://github.com/mariannedehasque/nf-pipelines/tree/main/nf-trim-merged-unmerged).

Made symlinks for NextFlow to the raw fastq files and renamed the symlinks with ```bash fix_symlinks.sh --apply``` (drop the ```--apply``` for a dry run). Script is now in `scripts/`. This scripts checks for corrupt files in ```/archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/fq_raw/``` (see ```2026-05-12_corrupt_fastq_report.txt```, now in `output/`), symlinks to them if not corrupt, checks for the equivalent in ```fq_fp1_clmp_fp2_fqscrn_rprd``` if they are corrupted, and links there instead if possible. It also checks that the R1 and R2 symlinks both point to the same directory (either raw or repaired). 

Created the list of filenames from within `nf-trim-merged-unmerged/`:
```
ls ./data/symlinks/*fastq.gz | xargs -n1 basename | cut -d "_" -f1,2,3 | uniq > ./inputfiles/fastq_filenames.txt
```

Created softlinks to reference and repma file
```
ln -s /archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/GenErode/reference/reference.ssl.Cvi20k_rename.fasta ./data/reference/
ln -s /archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/GenErode/reference/reference.ssl.Cvi20k_rename.fasta.* ./data/reference/
ln -s /archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/GenErode/reference/reference.ssl.Cvi20k_rename.dict ./data/reference/
ln -s /archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/GenErode/reference/reference.ssl.Cvi20k_rename.repma.bed ./data/reference/
```
Downloaded some files missing from clone of repo (not sure why were they missing):
```
wget https://raw.githubusercontent.com/mariannedehasque/nf-pipelines/refs/heads/main/nf-trim-merged-unmerged/main.nf
wget https://raw.githubusercontent.com/mariannedehasque/nf-pipelines/refs/heads/main/nf-trim-merged-unmerged/environment.yml
wget https://raw.githubusercontent.com/mariannedehasque/nf-pipelines/refs/heads/main/nf-trim-merged-unmerged/nextflow.config
```
Added entries from `nf-pipelines/.gitignore` to this repo's `.gitignore`.

Edited main.nf with reference name and length of historical reads (121 bp). Length derived from [Jem's MultiQC report](https://github.com/philippinespire/pire_chromis_viridis_lcwgs/blob/main/2nd_sequencing_run/fq_fp1_clmp_fp2_fqscrn_rprd/fqc_rprd_report.html)

Edited `main.nf` to use `bwa mem` for reads >80bp, and existing `bwa aln` for shorter reads. Edited `nextflow.config` and `environment.yml` to run mapdamage and amber.

Added mapdamage step to `main.nf` and turned on amber.

Ran nf-trim-merged-unmerged pipeline from within `nf-trim-merged-unmerged/`:
```
tmux new -s nextflow
bash
module load container_env
module load nextflow
nextflow run main.nf -profile standard -resume
```
Type `Ctrl-B` and then `D` to leave tmux.   
To rejoin tmux:
```
tmux a -t nextflow
```

### 7.1 Diagnosing excessive soft-clipping
Mapdamage plots show up to 30% softclipping on read ends. Does not appear to adapter, based on FastQC reports. Wrote a script to summarize the amount of soft-clipping in each scaffold of an individual:
```
scripts/softclip_by_scaffold_primary.sh /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-merged-unmerged/results/data/bam/CviAPal014.merged.L121.realn.bam output

scripts/softclip_by_scaffold_primary.sh /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-merged-unmerged/results/data/bam/CviCPal001.merged.L121.realn.bam output
```
Inspecting reveals three scaffolds in CviAPal014 (25 in CviCPal001) with softclipping in >60% of reads and >30% of bases. CviCPal001 has up to 56% bases softclipped on a scaffold. Put output in `output/`.

Manual inspection with IGV suggests lots of softclipping in areas of high depth. Calculate depth and softclipping in 500 bp windows:
```
scripts/softclip_by_window_primary.sh /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-merged-unmerged/results/data/bam/CviAPal014.merged.L121.realn.bam 500 output

scripts/softclip_by_window_primary.sh /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-merged-unmerged/results/data/bam/CviCPal001.merged.L121.realn.bam 500 output
```
Manual plots of the output don't reveal any obvious problems with high-depth regions. Some have a high fraction of soft-clipping, but low-depth regions do, too. Maybe misplaced bases during library prep caused by "over-chewing" by the KAPA end repair enzyme or by single-stranded overhangs from problems during enzymatic fragmentation. Unclear. Soft-clipping seems like a reasonable solution in any case.

### 7.2 Clean up
Removed the `nf-pipelines/nf-trim-merged-unmerged/work` directory.

## 8. ANGSD structure and diversity
Malin, 2026 June. Working in `/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-angsd-diversity`.

Wrote `scripts/make_samplesheet_from_bam.sh` to create `inputfiles/samplesheet.csv`:
```
scripts/make_samplesheet_from_bam.sh /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-merged-unmerged/results/data/bam inputfiles/samplesheet.csv
```

Created the ANGSD sites file next to the reference in `nf-trim-merged-unmerged`. It needs an input file with the starting position shifted by 1 (first command), then load angsd and run the sites command:
```
awk '{print $1"\t"$2+1"\t"$3}' /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-merged-unmerged/data/reference/reference.ssl.Cvi20k_rename.repma.bed > /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-merged-unmerged/data/reference/reference.sslCvi20k_rename.repma.angsd.txt

module load container_env
module load angsd/0.940

crun angsd sites index /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-merged-unmerged/data/reference/reference.sslCvi20k_rename.repma.angsd.txt
```

Create the BAM inputfile with a custom script:
```
scripts/list_bam_paths.sh /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-merged-unmerged/results/data/bam inputfiles/bam_list.txt
```

Create the contigs file by reading the contig names from the reference genome and removing the > character:
```
grep '^>' /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-merged-unmerged/data/reference/reference.ssl.Cvi20k_rename.fasta | sed 's/^>//' > inputfiles/contig_list.txt
```

Update the two main conda paths in test and standard from `nextflow.config` to specify `conda = "/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/environments/nf-angsd-diversity.yml"` instead of the existing path to marianne's yml file.

Calculate the expected coverage from the dpstats files output by amber in the nf-trim-mergd-unmerged pipeline (outputs 202.6):
```
awk '{s+=$1} END{print s}' /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-merged-unmerged/results/results/stats/Cvi*.bam.dpstats.txt
```

Update the parameters in `main.nf` for this run:
- locations of samplesheet.csv, contig_list.txt, reference fasta, bed file
- species code Cvi
- set maxdepth to 10x the expected depth = 2026
- set minind to 61, which is 70% of the 87 individuals we have

Odd that there is no parameter for the bam list.

Start nextflow in my existing tmux window, which already has bash activated and the container_env and nextflow modules loaded (see [Step 7](#7-nextflow-trimming)):
```
tmux a -t nextflow
cd ../nf-angsd-diversity # switch to the new pipeline
nextflow run main.nf -profile standard
```

The [PCA](nf-pipelines/nf-angsd-diversity/results/PCAngsd/Cvi.pcangsd.plot.pdf) shows substantial divergence along PC1 and PC2, but the log file reveals the SVD algorithm didn't converge.

### 8.1 PCA and admixture
Run pcangsd with more iterations and ask it to calculate admixture proportions, then plot:
```
module load container_env ngsTools
crun pcangsd --iter 500 -b /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-angsd-diversity/results/GL/Cvi.beagle.gz -t 8 -o output/Cvi.pcangsd --admix

crun Rscript nf-pipelines/nf-angsd-diversity/scripts/pcangsd.R --cov output/Cvi.pcangsd.cov --samplesheet nf-pipelines/nf-angsd-diversity/inputfiles/samplesheet.csv --species Cvi --out output/Cvi.pcangsd.plot.pdf
```
The [new PCA](output/Cvi.pcangsd.plot.pdf) has converged. It is effectively the same as before. Suggests K=3.

Plot the admixture proportions with a new script:
```
crun Rscript scripts/plot_admixture.R output/Cvi.pcangsd.admix.3.Q nf-pipelines/nf-angsd-diversity/inputfiles/samplesheet.csv output/Cvi.pcangsd.admix.pdf
```

The [admixture plot](output/Cvi.pcangsd.admix.pdf) mostly separates by era, as expected, but 8 contemporary individuals group with the historical ones.

### 8.2 FST historical-modern
Wrote a couple scripts to calc genome-wide and windowed (50kb windows, 10kb steps) fst historical vs. modern: scripts/calc_fst_modern_historic.sbatch, which calls scripts/calc_fst_modern_historic.sh
```
sbatch scripts/calc_fst_modern_historic.sbatch
```

Created `output/fst_historic_vs_modern/` with the output files, including 2D sfs and windowed fsts. The weighted global FST is 0.110033. High!

### 8.3 Plot ANGSD diversity
Plot the mean pi values by historical vs. modern with whiskers for the 95% CIs. Uses a new custom script that calculates per-site pi and bootstraps to get 95% CIs:
```
bash
module load container_env R
crun Rscript scripts/plot_tp_historic_modern.R nf-pipelines/nf-angsd-diversity/results/angsd_pop_theta/CviAPal_historic.pestPG nf-pipelines/nf-angsd-diversity/results/angsd_pop_theta/CviCPal_modern.pestPG output/tp_historic_vs_modern_mean_ci.png
```

The [plot of pi](output/tp_historic_vs_modern_mean_ci.png) suggests higher diversity in the modern samples. Before we think too hard on this, let's check for species identity. If modern is mixing two species (see the admixture plot), that would explain higher diversity.

## 9. MitoZ
Malin Pinsky, June 2026. Working in `/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/`

Put the mitoz scripts in the `scripts/` directory:
```
wget https://github.com/philippinespire/pire_lcwgs_data_processing/raw/refs/heads/main/scripts/MitoZ_wahab/runMitoZ_array_lcwgs.bash
wget https://github.com/philippinespire/pire_lcwgs_data_processing/raw/refs/heads/main/scripts/MitoZ_wahab/runMitoZ_array_lcwgs.sbatch
```

Bash script sets up an array of sbatch calls, one for each individual. I modified the filename pattern in .bash and .sbatch scripts to match Cvi: *.clmp.fp2_r1.fq.gz. Also modified scripts to take OUTDIR as the 2nd argument and the number of individuals to run as a 4th argument (0 means all individuals). Modifed the sbatch script to check for an existing output directory; skip the individual if one exists. 

Run on one individual to set up the mitoz database:
```
bash
scripts/runMitoZ_array_lcwgs.bash /archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/fq_fp1_clmp_fp2 /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/mitoz 32 1
```

Then run on the full 2nd sequencing run:
```
scripts/runMitoZ_array_lcwgs.bash /archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/fq_fp1_clmp_fp2 /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/mitoz 32 0
```

Cancelled some jobs that were running for >7 hrs, removed their output directories in mitoz/, and ran this again to try again. Ran into corrupted .fq.gz format issues for 26 of them and didn't try to solve.

Next, move the sbatch .out files into MitoZ as well, then move mitoz into output/
```
mkdir mitoz/logs
mv MitoZ-*.out mitoz/logs
mv mitoz temp/
```

Summarize the COX1 sequences across all individuals into `output/MitoZ_output.fasta`:
```
grep 'COX1' temp/mitoz/Cvi*/Cvi*.result/*.cds | sed 's/_MitoZ.*//g' | sed 's/^/>/g' > MitoZ_labels
grep -A1 'COX1' temp/mitoz/Cvi*/Cvi*.result/*.cds | grep 'cds-' | sed 's/.*cds-//g' > MitoZ_seqs
paste -d '\n' MitoZ_labels MitoZ_seqs > output/MitoZ_output.fasta
rm MitoZ_labels
rm MitoZ_seqs
```

Have 32 COX1 sequences (18 A, 14 C). Blasted them with a script:
```
sbatch --job-name=blastn_remote_top2 --cpus-per-task=4 --mem=16G --time=24:00:00 --wrap="bash scripts/blastn_remote_top2.sh output/MitoZ_output.fasta output/MitoZ_output_vs_nt"
```
Note this returned warnings from NCBI about many requests from this IP address, which likely slowed down response times. Could use megablast instead of blastn in the future? Or download the nt database locally.

Detailed results in `output/MitoZ_output_vs_nt_all.tsv` with columns  qseqid (query sequence ID from the FASTA header), sacc (subject accession), pident (percent identity), length (alignment length), qcovs (query coverage per subject), evalue (E-value), bitscore (BLAST bit score), staxids (subject taxonomic IDs), sscinames (subject scientific names), stitle (subject title/description). 

Note that APal_002 had two queries, and they both have the same qseqid. A better labeling command would have been `grep 'COX1' temp/mitoz/Cvi*/Cvi*.result/.cds | sed 's/_MitoZ.//g' | awk '{k[$0]++; print ">"$0"_hit"k[$0]}' > MitoZ_labels`, which would have produced unique qseqids.

Trimmed to top two in `MitoZ_output_vs_nt_top2.tsv` and with headers and labeled more clearly in `MitoZ_output_vs_nt_top2_compact.tsv`.

1 of 17 APal individuals matches Cvi. Others match bacteria.  
7 of 14 CPal individuals match Cvi (CPal_002, CPal_005, CPal_030, CPal_031, CPal_032, CPal_052, and CPal_093). Others match _Chromis atripectoralis_.

Bingo! The CPal individuals that matched Cvi are the "purple" individuals in the [admixture plot](output/Cvi.pcangsd.admix.pdf). The others are _C. atripectoralis_. We mostly collected _C. atripectoralis_.

## 10. ANGSD structure and diversity with only viridis


## Future
- Sliding window FST
- ACER selection scan
