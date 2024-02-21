# *Chromis viridis* lcWGS, test lane
### Kevin Labrador
#

## 1. Pre-processing


This follows the instructions from [pire_fq_gz_processing](https://github.com/philippinespire/pire_fq_gz_processing/blob/main/README.md).


<details>
        <summary>0. Set up directory</summary>

Run by klabrador on 2023-04-28

- Go to PIRE lcwgs directory on Wahab and create the appropriate subdirectory for the species.

```bash
cd /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing
mkdir -p chromis_viridis/1st_sequencing_run/
```

- Go to the assigned working directory and create the following subdirs:

```bash
mkdir fq_raw fq_fp1 fq_fp1_clmp fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn fq_fp1_clmp_fp2_fqscrn_rprd
```

</details>


<details>
        <summary>1. Download data from TAMUCC grid</summary>

Run by klabrador on 2023-04-28

- dataset location @ https://gridftp.tamucc.edu/genomics/20230425_PIRE-Cvi-lcwgs-testlane/

```bash
cd ./fq_raw
sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/gridDownloader.sh . https://gridftp.tamucc.edu/genomics/20230425_PIRE-Cvi-lcwgs-testlane/
```

</details>


</details>


<details>
        <summary>2. Proofread the decode file(s)</summary>

Run by klabrador on 2023-04-28

- No. of samples = 104
- All files look good.


</details>


<details>
        <summary>3. Edit the decode file</summary>

Run by klabrador on 2023-04-28

- No need to edit the decode file.


</details>



<details>
        <summary>4. Make a copy of the `fq_raw` files prior to renaming </summary>

Run by klabrador on 2023-04-28

```
mkdir /RC/group/rc_carpenterlab_ngs/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run
mkdir /RC/group/rc_carpenterlab_ngs/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run/fq_raw

# on wahab replace <yourPireDirPath> with /home/e1garcia/shotgun_PIRE
cd /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run/fq_raw
screen cp ./* /RC/group/rc_carpenterlab_ngs/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/fq_raw 

# `ctrl-a`  and then `d` to detach from the `screen` terminal

# look at your screen jobs running
screen -ls
```

</details>


<details>
        <summary>5. Perform a renaming dry run</summary>

Run by klabrador on 2023-04-28

```
cd /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run/fq_raw

salloc
bash /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/renameFQGZ.bash Cvi_LCWGS-TestLane_SequenceNameDecode.tsv
```
- Looks good.

</details>


<details>
        <summary>6. Rename the files for real</summary>

Run by klabrador on 2023-04-28

```
bash /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/renameFQGZ.bash Cvi_LCWGS-TestLane_SequenceNameDecode.tsv rename
```

</details>


<details>
        <summary>7. Check data quality</summary>

Run by klabrador on 2023-04-28

```
cd /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run

#sbatch Multi_FASTQC.sh "<indir>" "<mqc report name>" "<file extension to qc>"
sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/Multi_FASTQC.sh "fq_raw" "fqc_raw_report"  "fq.gz"

```

- job submitted; job ID not recorded.
- job finished successfully

</details>


<details>
        <summary>8. First trim</summary>

Run by klabrador on 2023-04-28

```
cd /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run
sbatch ../../../pire_fq_gz_processing/runFASTP_1st_trim.sbatch fq_raw fq_fp1
```

- job submitted: 1508250
- job finished successfully

</details>


<details>
        <summary>9a. Remove duplicates</summary>

Run by klabrador on 2023-04-28

```
cd /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run
bash ../../../pire_fq_gz_processing/runCLUMPIFY_r1r2_array.bash fq_fp1 fq_fp1_clmp /scratch/hpc-0289 20
```

- job submitted: 1508297
- job finished successfully

</details>


<details>
        <summary>9b. Check duplicate removal success</summary>

Run by klabrador on 2023-04-29

```
cd /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run

salloc
enable_lmod
module load container_env mapdamage2

crun R < /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/checkClumpify_EG.R --no-save
exit
```

- Clumpify succesfully worked on all samples.

</details>


<details>
        <summary>9c. Generate metadata on deduplicated FASTQ files</summary>

Run by klabrador on 2023-04-29

```
cd /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run

sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/Multi_FASTQC.sh "fq_fp1_clmp" "fqc_clmp_report"  "fq.gz"

```

- job submitted: 1510092
- job finished successfully

</details>


<details>
        <summary>10. Second Trim</summary>

Run by klabrador on 2023-05-01

```
cd /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run

# Run cssl script for lcwgs.
sbatch ../../../pire_fq_gz_processing/runFASTP_2_cssl.sbatch fq_fp1_clmp fq_fp1_clmp_fp2
```

- job submitted: 1511955
- job finished successfully

</details>

<details>
        <summary>11. Decontaminate files</summary>

Run by klabrador on 2023-05-01

```
cd /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run

bash ../../../pire_fq_gz_processing/runFQSCRN_6.bash fq_fp1_clmp_fp2 fq_fp1_clmp_fp2_fqscrn 20
```

- job submitted: 1511991

Confirm if all files were successfully completed. 
```
# Check file count.
ls fq_fp1_clmp_fp2_fqscrn/*tagged.fastq.gz | wc -l
ls fq_fp1_clmp_fp2_fqscrn/*tagged_filter.fastq.gz | wc -l 
ls fq_fp1_clmp_fp2_fqscrn/*screen.txt | wc -l
ls fq_fp1_clmp_fp2_fqscrn/*screen.png | wc -l
ls fq_fp1_clmp_fp2_fqscrn/*screen.html | wc -l
## n = 208 for all outfiles

# Check for errors
grep 'error' slurm-fqscrn.*out
grep 'No reads in' slurm-fqscrn.*out
## No errors/no-reads found.
```

- Run MultiQC for fqscrn output.
```
cd /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run
sbatch ../../../pire_fq_gz_processing/runMULTIQC.sbatch fq_fp1_clmp_fp2_fqscrn fastq_screen_report
```

- job submitted: 1572125
- job finished successfully

</details>


<details>
        <summary>12. Execute `runRepair`</summary>

Run by klabrador on 2023-05-02

```
cd /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run

sbatch ../../../pire_fq_gz_processing/runREPAIR.sbatch fq_fp1_clmp_fp2_fqscrn fq_fp1_clmp_fp2_fqscrn_rprd 40
```

- job submitted: 1573171
- job finished successfully


Run `Multi_FASTQC.sh`
```
cd /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run
sbatch /home/e1garcia/shotgun_PIRE/pire_fq_gz_processing/Multi_FASTQC.sh "./fq_fp1_clmp_fp2_fqscrn_rprd" "fqc_rprd_report" "fq.gz"
```

- job submitted: 1590138
- job finished successfully

</details>


<details>
        <summary>13. Clean Up</summary>

```
mkdir logs
cd /home/e1garcia/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run
mv *out logs/
```

Run by klabrador on 2023-08-01

I have to clean up the repo to make more space. I followd the **CLEANING UP** section in the [pire_ssl_data_processing](https://github.com/philippinespire/pire_ssl_data_processing).

1. Check initial file size

```
du -h | sort -rh > Cvi_lcwgs_beforeDeleting_InterimFiles

```
- File size: 262 Gb

2. Make a copy of important files in the RC.

```
# Get inside a node
salloc

# Confirm that the RC directory for the species exist
ls /RC/group/rc_carpenterlab_ngs/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run

# Contaminated file
cp -R fq_fp1_clmp_fp2 /RC/group/rc_carpenterlab_ngs/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run

# Decontaminated file
cp -R fq_fp1_clmp_fp2_fqscrn_rprd /RC/group/rc_carpenterlab_ngs/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run

# Check if file size is equal between the origin and destination directories.
du -sh *
du -sh /RC/group/rc_carpenterlab_ngs/shotgun_PIRE/pire_lcwgs_data_processing/chromis_viridis/1st_sequencing_run/*

```
- File sizes are equal; files were copied successfully.

3. Delete unneeded files

```
# create log file before removing
ls -ltrh *raw*/*fq.gz > deleted_files_log
ls -ltrh *fp1/*fq.gz >> deleted_files_log
ls -ltrh *clmp/*fq.gz >> deleted_files_log
ls -ltrh *fqscrn/*fastq.gz >> deleted_files_log

# remove files
rm *raw*/*fq.gz
rm *fp1/*fq.gz
rm *clmp/*fq.gz
rm *fqscrn/*fastq.gz
```

4. Document the new size of directories
```
du -h | sort -rh > Cvi_lcwgs_afterDeleting_InterimFiles
```
- File size = 71 Gb

Freed space: 262 Gb - 71 Gb =  191 Gb 


5. Move the cleaning files into the logs dir.

```
mv deleted_files_log logs
``` 
</details>





