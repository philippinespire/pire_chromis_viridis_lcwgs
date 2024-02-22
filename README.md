# Directory Structure

| directory-name | date-created | dataset |
| -------------- | ------------ | ------- |
| 1st_sequencing_run | 2023-04-28 | test lane |
| 2nd_sequencing_run | 2023-06-16 | Alb & Contemp sequences from Palapag, Northern Samar|

# Chromis viridis lcWGS

---
Jem Baldisimo
---
This repository outlines the roadmap we followed to move *Chromis viridis* through the [Low Coverage Whole Genome Sequencing Pipeline](https://github.com/philippinespire/pire_lcwgs_data_processing). This provides steps taken & analysis we did to gain insight on the historical population demography of this very popular aquarium fish, also known as the green chromis.

---

## 1. Completed fq.gz pre-processing

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

I copied the best assembly from the ssl pipeline (see link below) and renamed it to reference.ssl.Cvi.fasta that was copied into the refGenome folder
```
/home/e1garcia/shotgun_PIRE/pire_ssl_data_processing/chromis_viridis/SPAdes_Cvi-CPal-A_decontam_R1R2_noIsolate/scaffolds.fasta
```

---
## 3. Mapped & Filtered BAM files

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

```
#navigate to mkBAM folder
# mappedReadStats.sbatch <Path to BAM file dir> <coverageMappedReads>
sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/mappedReadStats.sbatch . coverageMappedReads
```

---
## 5. Visualised results using the Process Sequencing Metadata Repo
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

