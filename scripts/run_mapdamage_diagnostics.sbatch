#!/bin/bash
################################################################################
# SCRIPT NAME: run_mapdamage_diagnostics.sh
# Malin Pinsky, 2026, with assistance from Gemini
#
# DESCRIPTION: Loops through BAM files matching 'CviCPal*.bam' to run mapDamage
#              diagnostics (without quality score rescaling). It isolates the
#              outputs into a temporary directory before renaming and moving 
#              them to a flat directory with the sample name prepended.
#
# REQUIREMENTS:
#   - mapDamage2 (and its R dependencies) must be available in your environment.
#   - A valid reference genome FASTA file used for the original mapping.
#   - Slurm Workload Manager (for cluster submission).
#
# USAGE:
#   1. Open this script and edit the following variables:
#      - 'REF_FASTA'            : Set to the path of your reference genome.
#   2. Submit the job to the Slurm queue using:
#      sbatch run_mapdamage_diagnostics.sbatch
#
# OUTPUTS:
#   - Results will be saved in directory specied by OUTPUT_BASE (eg, output/mapdamage/)
#   - Files will be named like: CviCPal001_Fragmisincorporation_plot.pdf
#   - Slurm log will be saved as: logs/mapdamage_<jobID>.out
################################################################################

#SBATCH -n 1
#SBATCH --time 06:00:00
#SBATCH --cpus-per-task=2            # mapDamage doesn't heavily multithread, but 2 cores helps with I/O
#SBATCH --job-name=mapd_diag
#SBATCH --output=/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/logs/mapdamage_%j.out
#SBATCH --error=/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/logs/mapdamage_%j.err
#SBATCH --mem=16G 

# Exit immediately if a command exits with a non-zero status
set -e

# Change this to the reference FASTA file path!
#REF_FASTA="nf-pipelines/nf-trim-generode/data/reference/reference.ssl.Cvi20k_rename.fasta" # first run with our in-house genome
REF_FASTA="/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/data/GCA_051013605.1_ASM5101360v1_genomic_20kb.fna" # second run with the iridian genome

# Define directories
# INPUT_DIR="nf-pipelines/nf-trim-generode/results/data/bam" # first run with our in-house genome
INPUT_DIR="nf-pipelines/nf-trim-generode/results-iridian/data/bam" # second run with the iridian genome
# OUTPUT_BASE="output/mapdamage" # first run with our in-house genome
OUTPUT_BASE="output/mapdamage-iridian" # second run with the iridian genome

# Load mapdamage2 module
module load container_env mapdamage2

# Ensure the main output directory exists
mkdir -p "$OUTPUT_BASE"

# Loop through the matching BAM files
for bam in "${INPUT_DIR}"/CviCPal*.bam; do
    # Gracefully handle the loop if no matching files are found
    [ -e "$bam" ] || { echo "No files matching CviCPal*.bam found."; exit 1; }

    # Extract the sample name (filename without .bam extension)
    sample_name=$(basename "$bam" .bam)
    
    # Create a temporary subdirectory for this specific run
    tmp_dir="${OUTPUT_BASE}/tmp_${sample_name}"
    mkdir -p "$tmp_dir"
    
    echo "----------------------------------------"
    echo "Running mapDamage for: ${sample_name}"
    echo "----------------------------------------"
    
    # Run mapDamage into the temporary folder
    crun mapDamage -i "$bam" -r "$REF_FASTA" -d "$tmp_dir" --no-stats
    
    # Loop through the generated files, prepend the sample name, 
    # and move them into the main output folder
    for file in "$tmp_dir"/*; do
        if [ -f "$file" ]; then
            filename=$(basename "$file")
            mv "$file" "${OUTPUT_BASE}/${sample_name}_${filename}"
        fi
    done
    
    # Clean up the empty temporary directory
    rmdir "$tmp_dir"

done

echo "========================================"
echo "All files processed and renamed successfully!"
echo "Results are inside: $OUTPUT_BASE"
echo "========================================"