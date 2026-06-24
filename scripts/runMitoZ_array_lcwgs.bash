#!/bin/bash 

#Script for running MitoZ over multiple lcwgs libraries
#Use contam libraries created after second trim (fp2)
#just needs four arguments:
#1)location of fq_fp1_clmp_fp2 folder
#2)output directory for MitoZ results
#3)number of nodes to use (set to 32 for a test run and it cranked through ~100 libraries in ~2 hrs)
#4)number of individuals to run (0 = all individuals)

#will execute MitoZ in parallel through the companion script runMitoZ_array.sbatch

#Pass in the maximum number of nodes to use at once
nodes=$3    # eg 32 - running 4 at a time?
num_inds=$4 # 0 = all individuals

if [[ -z ${num_inds} ]]; then
	echo "Usage: $0 <INDIR> <OUTDIR> <NODES> <NUM_INDIVIDUALS_OR_0_FOR_ALL>"
	exit 1
fi

INDIR=$1                 #example= /archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/fq_fp1_clmp_fp2
OUTDIR=$2                 #example= /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/mitoz
FQPATTERN=*.clmp.fp2_r1.fq.gz        #forward reads. it will automatically find the reverse reads as well, as long as they are in the same folder and have the same prefix naming style

SCRIPTPATH=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

all_samples=$(ls $INDIR/$FQPATTERN | sed -e 's/.clmp.fp2_r1\.fq\.gz//' -e 's/.*\///g')
all_samples=($all_samples)

total_samples=${#all_samples[@]}

if [[ ${num_inds} -eq 0 || ${num_inds} -ge ${total_samples} ]]; then
	array_max_idx=$((${total_samples}-1))
else
	array_max_idx=$((${num_inds}-1))
fi

sbatch --array=0-${array_max_idx}%${nodes} $SCRIPTPATH/runMitoZ_array_lcwgs.sbatch ${INDIR} ${OUTDIR}
