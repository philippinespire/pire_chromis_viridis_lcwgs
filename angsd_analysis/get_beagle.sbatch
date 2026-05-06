#!/bin/bash -l

#SBATCH --job-name=angsd_beagle
#SBATCH -o angsd_beagle-%j.out 
#SBATCH -p main 
#SBATCH -n 1  
#SBATCH -N 1  
#SBATCH --cpus-per-task=8 

enable_lmod
module load container_env angsd
export SINGULARITY_BIND=/home/e1garcia

OUTDIR=$1


cd ${OUTDIR}

crun angsd \
        -b bam_list_all_fullpath.txt \
        -anc /archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/GenErode/reference/reference.ssl.Cvi20k_rename.fasta \
        -out angsd_depth1_15_notrans_snplist \
        -doSaf 1 -noTrans 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -dumpCounts 1 \
        -P 8 \
	-sites global_snp_list_depth1_15_notrans.txt -rf global_snp_list_depth1_15_notrans.chrs 
