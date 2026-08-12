#!/bin/bash
# ==============================================================================
# Script Name: submit_acer.sh
# Author: Malin Pinsky, August 2026 (with Gemini assistance)
# Description: SLURM batch submission wrapper for the ACER selection scan (run_acer.R).
#              This script dynamically passes all command-line arguments directly
#              to the underlying R script using "$@".
#
# Usage: 
#   sbatch submit_acer.sh [arguments for run_acer.R]
#
# Required Arguments:
#   --hist_mafs     Comma-separated list of historical maf files
#   --mod_mafs      Comma-separated list of modern maf files
#   --region_names  Comma-separated list of region names
#
# Optional Arguments (with defaults):
#   --out_dir       Output directory (default: ./results/selection)
#   --helpers       Path to acer_helpers.R (default: acer_helpers.R)
#   --ne_generations Number of generations for Ne calculation (default: 114)
#   --test_gen_start Start generation for testing (default: 0)
#   --test_gen_end   End generation for testing (default: 113)
#   --fdr_cutoff    FDR cutoff (default: 0.05)
#   --max_rounds    Max iterations (default: 20)
#   --n_boot        Number of bootstrap replicates (default: 1000)
#   --min_ind       Min individuals per group (default: 4)
#
# Example Submission:
#   sbatch submit_acer.sh \
#     --hist_mafs=/path/hist1.mafs.gz \
#     --mod_mafs=/path/mod1.mafs.gz \
#     --region_names=Pop1 \
#     --out_dir=./results/Pop1
# ==============================================================================

#SBATCH --job-name=acer                  # Name of your job
#SBATCH --output=logs/acer_%j.out        # Standard output log (%j will be replaced by the job ID)
#SBATCH --error=logs/acer_%j.err         # Standard error log
#SBATCH --time=23:59:00                  # Maximum time requested (HH:MM:SS)
#SBATCH --ntasks=1                       # Number of tasks (usually 1 for a single R script)
#SBATCH --cpus-per-task=1                # Number of CPU cores to allocate
#SBATCH --mem=16G                        # Amount of RAM to allocate



# Load the required container environment and R module for Wahab
module load container_env R 

# Run the R script using crun
# "$@" automatically forwards any arguments passed during the sbatch submission
crun Rscript scripts/run_acer.R "$@"