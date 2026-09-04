#!/bin/bash
#SBATCH --job-name=prlgs_angsd
#SBATCH --output=logs/prlgs_angsd_%j.out
#SBATCH --error=logs/prlgs_angsd_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# ==============================================================================
# Script Name   : paralogs_01_angsd.sbatch
# Description   : Population-level paralog screening for Chromis viridis (Cvi).
#                 Uses ANGSD to evaluate site statistics, detect severe excess
#                 heterozygosity (F < -0.3, p < 1e-4), and flag high-coverage 
#                 depth outliers (>1.75x mean depth).
# Environment   : Wahab HPC Cluster
# Pipeline Step : Step 1 of 3 (Paralog Filtering Workflow)
#
# Malin Pinsky, August 2026, with support from Gemini
#
# Inputs        :
#   - Reference FASTA (GCA_051013605.1_ASM5101360v1_genomic_20kb.fna)
#   - BAM Filelist   (output/paralogs/bam.filelist containing paths to modern & rescaled BAMs)
#
# Outputs       :
#   - TEMP_DIR/cvi_angsd_qc.*         Large raw ANGSD stats & depth outputs
#   - OUT_DIR/hwe_excess_het.bed      1-indexed, closed interval, merged BED of excess het loci (ANGSD format)
#   - OUT_DIR/depth_outliers.bed      1-indexed, closed interval, merged BED of high-depth loci (ANGSD format)
#
# Dependencies  : module load container_env angsd
# ==============================================================================

module load container_env angsd

REF="/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/data/GCA_051013605.1_ASM5101360v1_genomic_20kb.fna" # reference genome
OUT_DIR="output/paralogs" # where to write results
TEMP_DIR="temp/paralogs"   # directory for large intermediate ANGSD files
BAM_LIST="${OUT_DIR}/bam.filelist" # list of BAM files to process
OUT_PREFIX="${TEMP_DIR}/cvi_angsd_qc" # prefix for output files

mkdir -p logs ${OUT_DIR} ${TEMP_DIR} # create output directories if they don't exist

# Step 1: Run ANGSD to calculate site statistics, HWE, and read depth
crun.angsd angsd \
    -b ${BAM_LIST} \
    -ref ${REF} \
    -out ${OUT_PREFIX} \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
    -minMapQ 25 -minQ 20 \
    -GL 1 -doMajorMinor 1 -doMaf 1 -doCounts 1 \
    -doDepth 1 -dumpCounts 1 \
    -doHWE 1 -doSnpStat 1 -HWE_pval 1 \
    -nThreads ${SLURM_CPUS_PER_TASK}

# Step 2: Extract loci with excess heterozygosity (F < -0.3 and p < 1e-4)
# Combines adjacent positions into 1-indexed intervals [start, end]
zcat ${OUT_PREFIX}.hwe.gz | awk '
    NR>1 && $5 < -0.3 && $6 < 1e-4 {
        chr = $1; pos = $2
        if (chr == cur_chr && pos == end + 1) {
            end = pos
        } else {
            if (cur_chr != "") print cur_chr "\t" start "\t" end
            cur_chr = chr; start = pos; end = pos
        }
    }
    END {
        if (cur_chr != "") print cur_chr "\t" start "\t" end
    }' > ${OUT_DIR}/hwe_excess_het.bed

# Step 3: Extract coverage depth outliers (>1.75x mean coverage across loci)
grep -i "depth" ${OUT_PREFIX}.arg || true # Log check

# Step 4: Compute genome-wide mean depth (pos.gz: [Chr] [Pos] [Total_Depth])
MEAN_DEPTH=$(zcat ${OUT_PREFIX}.pos.gz | awk '{sum+=$3} END {print sum/NR}')
CUTOFF=$(awk -v mean="${MEAN_DEPTH}" 'BEGIN {print mean * 1.75}')

echo "Genome mean depth: ${MEAN_DEPTH}x. Outlier cutoff set at: ${CUTOFF}x."

# Step 5: Extract high-depth loci, combining adjacent positions into 1-indexed intervals
zcat ${OUT_PREFIX}.pos.gz | awk -v cutoff="${CUTOFF}" '
    $3 > cutoff {
        chr = $1; pos = $2
        if (chr == cur_chr && pos == end + 1) {
            end = pos
        } else {
            if (cur_chr != "") print cur_chr "\t" start "\t" end
            cur_chr = chr; start = pos; end = pos
        }
    }
    END {
        if (cur_chr != "") print cur_chr "\t" start "\t" end
    }' > ${OUT_DIR}/depth_outliers.bed