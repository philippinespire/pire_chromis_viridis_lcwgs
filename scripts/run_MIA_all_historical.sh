#!/usr/bin/env bash

# submits an sbatch job for each historical individual in
# /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-generode/data/symlinks/CviAPal*_Ex1_L4_R1.fastq.gz
# Malin Pinsky 2026 with help from Google Gemini

# Define paths
FASTQ_DIR="/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-generode/data/symlinks"
SBATCH_SCRIPT="scripts/run_MIA.sbatch"

# Double-check that the sbatch script actually exists
if [[ ! -f "$SBATCH_SCRIPT" ]]; then
    echo "ERROR: SBATCH script not found at $SBATCH_SCRIPT"
    echo "Please run this script from the directory containing your 'scripts' folder."
    exit 1
fi

echo "Scanning for forward reads in: $FASTQ_DIR"
echo "--------------------------------------------------"

# Loop through each forward read file matching the pattern
for R1_FILE in "${FASTQ_DIR}"/CviAPal*_Ex1_L4_R1.fastq.gz; do
    # Ensure the glob matched actual files before trying to submit
    if [[ -f "$R1_FILE" ]]; then
        FILE_NAME=$(basename "$R1_FILE")
        echo "Submitting job for: $FILE_NAME"
        
        # Submit the job to Slurm, passing the current forward read file
        sbatch "$SBATCH_SCRIPT" "$R1_FILE"
    else
        echo "No files matching 'CviAPal*_Ex1_L4_R1.fastq.gz' were found."
    fi
done

echo "--------------------------------------------------"
echo "All jobs submitted!"