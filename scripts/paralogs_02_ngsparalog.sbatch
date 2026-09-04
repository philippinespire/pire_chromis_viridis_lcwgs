#!/bin/bash
#SBATCH --job-name=ngsParalog
#SBATCH --output=logs/ngsParalog_%j.out
#SBATCH --error=logs/ngsParalog_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=36:00:00

# ==============================================================================
# Script Name   : paralogs_02_ngsparalog.sbatch
# Description   : Aligns modern FASTQ reads without MAPQ filtering to generate
#                 unfiltered BAMs (retaining MAPQ=0 multi-mappers). Runs 
#                 ngsParalog calcLR on candidate excess heterozygosity loci.
#                 Likelihood-ratio testing for paralogous genomic loci in 
#                 Chromis viridis (Cvi). Streams an unfiltered samtools mpileup
#                 (q=0, Q=0) into ngsParalog calcLR to detect multi-mapping read
#                 signatures and extracts candidate paralogs (LR > 10).
# Environment   : Wahab HPC Cluster
# Pipeline Step : Step 2 of 3 (Paralog Filtering Workflow)
#
# Malin Pinsky, August 2026, with support from Gemini
#
# Inputs        :
#   - Reference FASTA (GCA_051013605.1_ASM5101360v1_genomic_20kb.fna)
#   - FASTQ Directory (Directory containing Nextflow trimmed FASTQs)
#   - Compiled Binary(Path to local ngsParalog executable)
#   - BED file of candidate loci (hwe_excess_het.bed) from running paralogs_01_angsd.sbatch
#
# Outputs       :
#   - $OUT_DIR/cvi_ngsparalog.lr.txt   (Per-site likelihood ratio scores)
#   - $OUT_DIR/ngsparalog_sites.bed    (BED file of candidate loci, LR > 10)
#
# Dependencies  : module load container_env samtools bwa
#                 parallel
#                 Compiled ngsParalog binary executable
# ==============================================================================

module load container_env samtools bwa

ulimit -s unlimited

OUT_DIR="output/paralogs"
TEMP_DIR="temp/paralogs"
UNFILTERED_BAM_DIR="${TEMP_DIR}/unfiltered_bams" # where to put the unfiltered bam files
FASTQ_DIR="/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-generode/results-iridian/data/fastq" # Directory containing Nextflow trimmed FASTQs
REF="/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/data/GCA_051013605.1_ASM5101360v1_genomic_20kb.fna"
NGSPARALOG_BIN="/archive/carpenterlab/pire/softwares/ngsParalog/ngsParalog"

mkdir -p logs ${OUT_DIR} ${TEMP_DIR} ${UNFILTERED_BAM_DIR}

# ------------------------------------------------------------------------------
# 1. Check for candidate loci in hwe_excess_het.bed
# ------------------------------------------------------------------------------
INPUT_BED="${OUT_DIR}/hwe_excess_het.bed"

if [ ! -s "${INPUT_BED}" ]; then
    echo "NOTICE: ${INPUT_BED} is empty or missing. No candidate loci to test."
    touch ${OUT_DIR}/cvi_ngsparalog.lr.txt ${OUT_DIR}/ngsparalog_sites.bed
    exit 0
fi

# Convert 1-based closed ANGSD intervals [start, end] to standard 0-based BED [start-1, end]
# required for samtools mpileup -l
CONVERTED_BED="${TEMP_DIR}/hwe_excess_het_0based.bed"
awk 'OFS="\t" {print $1, $2-1, $3}' ${INPUT_BED} > ${CONVERTED_BED}

# ------------------------------------------------------------------------------
# 2. Parallel Alignment (4 jobs in parallel x 4 threads each = 16 CPUs)
# ------------------------------------------------------------------------------
align_single() {
    fq=$1; ref=$2; out_dir=$3
    sample=$(basename "$fq" | sed 's/_trimmed_merged.*//')
    out_bam="${out_dir}/${sample}_unfiltered.sorted.bam"
    
    if [ ! -f "${out_bam}" ]; then
        # 1. Map and keep all alignments (including secondary/supplemental, filtering only unmapped)
        # 2. Name-collate and fixmate to prepare for markdup
        # 3. Position sort and remove duplicates (-r flag in markdup)
        crun bwa mem -t 4 "${ref}" "${fq}" | \
            crun.samtools samtools view -b -F 4 -@ 2 - | \
            crun.samtools samtools collate -O -@ 2 - | \
            crun.samtools samtools fixmate -m -@ 2 - - | \
            crun.samtools samtools sort -m 3G -@ 2 - | \
            crun.samtools samtools markdup -r -@ 2 - "${out_bam}"
            
        crun.samtools samtools index "${out_bam}"
    fi
    echo "${out_bam}"
}
export -f align_single

# Find FASTQs and run 4 alignments concurrently
FASTQ_LIST=(${FASTQ_DIR}/*_trimmed_merged*.fastq.gz)
parallel -j 4 align_single {} "${REF}" "${UNFILTERED_BAM_DIR}" ::: "${FASTQ_LIST[@]}" > "${UNFILTERED_BAM_LIST}"

NUM_BAMS=$(wc -l < "${UNFILTERED_BAM_LIST}")
MIN_IND=$(( NUM_BAMS / 2 ))

# ------------------------------------------------------------------------------
# 3. Parallelized mpileup | ngsParalog across Chromosomes
# ------------------------------------------------------------------------------
# Split BED file into scaffold/chromosome chunks
SPLIT_DIR="${TEMP_DIR}/bed_splits"
mkdir -p "${SPLIT_DIR}" "${TEMP_DIR}/lr_splits"
awk -v dir="${SPLIT_DIR}" '{print > (dir "/" $1 ".bed")}' "${CONVERTED_BED}"

# Process each chromosome chunk in parallel across 16 CPUs
run_chunk() {
    bed_chunk=$1; bam_list=$2; ref=$3; min_ind=$4; bin=$5; temp_dir=$6
    chr=$(basename "${bed_chunk}" .bed)
    
    crun samtools mpileup -b "${bam_list}" -f "${ref}" -l "${bed_chunk}" -q 0 -Q 0 | \
        ${bin} calcLR -infile - -outfile - -minQ 20 -minind ${min_ind} > "${temp_dir}/lr_splits/${chr}.lr.txt"
}
export -f run_chunk

parallel -j ${SLURM_CPUS_PER_TASK} run_chunk {} "${UNFILTERED_BAM_LIST}" "${REF}" "${MIN_IND}" "${NGSPARALOG_BIN}" "${TEMP_DIR}" ::: ${SPLIT_DIR}/*.bed

# Combine parallel LR outputs into final output
cat ${TEMP_DIR}/lr_splits/*.lr.txt > ${OUT_DIR}/cvi_ngsparalog.lr.txt

# Extract high LR loci (LR > 10) in 1-based closed format
awk '$3 > 10 {print $1 "\t" $2 "\t" $2}' ${OUT_DIR}/cvi_ngsparalog.lr.txt > ${OUT_DIR}/ngsparalog_sites.bed