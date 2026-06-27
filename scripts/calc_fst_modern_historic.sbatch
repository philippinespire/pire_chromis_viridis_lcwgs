#!/usr/bin/env bash
#SBATCH --job-name=cvi_fst_modhist
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G

# Shell script to calculate Fst between modern and historic populations of Chromis viridis using ANGSD and realSFS.
# The script submits a slurm job to run the calculations, and handles command line arguments for specifying input and output directories, population prefixes, number of threads, and whether to use folded or unfolded site frequency spectra.
#
# Malin Pinsky 2026
# This script was written with AI assistance using GitHub Copilot
# with the GPT-5.3-Codex model.
#
# Usage:
#   bash scripts/calc_fst_modern_historic.sbatch --results-dir PATH --outdir PATH
#
# Default argument values:
#   --results-dir  nf-pipelines/nf-angsd-diversity/results/angsd_pop
#   --outdir       output/fst_historic_vs_modern
#   --historic     CviAPal_historic
#   --modern       CviCPal_modern
#   --threads      SLURM_CPUS_PER_TASK (default: 12)
#   --folded       enabled by default; use --unfolded to change it
#   --realSFS-cmd  crun realSFS

set -euo pipefail

# Run from repository root regardless of submit location.
repo_root="/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs"
cd "$repo_root"

# Default values for arguments
results_dir="${repo_root}/nf-pipelines/nf-angsd-diversity/results/angsd_pop"
outdir="${repo_root}/output/fst_historic_vs_modern"
historic_prefix="CviAPal_historic"
modern_prefix="CviCPal_modern"
threads="${SLURM_CPUS_PER_TASK:-12}"
folded=1
realsfs_cmd="crun realSFS"

# Alternatively, read from command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-dir)
      results_dir="$2"
      shift 2
      ;;
    --historic)
      historic_prefix="$2"
      shift 2
      ;;
    --modern)
      modern_prefix="$2"
      shift 2
      ;;
    --outdir)
      outdir="$2"
      shift 2
      ;;
    --threads)
      threads="$2"
      shift 2
      ;;
    --folded)
      folded=1
      shift
      ;;
    --unfolded)
      folded=0
      shift
      ;;
    --realSFS)
      realsfs_cmd="$2"
      shift 2
      ;;
    --realSFS-cmd)
      realsfs_cmd="$2"
      shift 2
      ;;
    *)
      echo "Error: Unknown argument '$1'" >&2
      exit 1
      ;;
  esac
done

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  logdir="${outdir}/logs"
  mkdir -p "$outdir" "$logdir"

  sbatch \
    --output="${logdir}/%x_%j.out" \
    --error="${logdir}/%x_%j.err" \
    "$0" \
    --results-dir "$results_dir" \
    --historic "$historic_prefix" \
    --modern "$modern_prefix" \
    --outdir "$outdir" \
    --threads "$threads" \
    --realSFS-cmd "$realsfs_cmd" \
    $(if [[ "$folded" -eq 1 ]]; then printf '%s ' --folded; else printf '%s ' --unfolded; fi)
  exit 0
fi

module load container_env angsd

bash scripts/calc_fst_modern_historic.sh \
  --results-dir "$results_dir" \
  --historic "$historic_prefix" \
  --modern "$modern_prefix" \
  --outdir "$outdir" \
  --threads "$threads" \
  --realSFS-cmd "$realsfs_cmd" \
  $(if [[ "$folded" -eq 1 ]]; then printf '%s ' --folded; else printf '%s ' --unfolded; fi)
