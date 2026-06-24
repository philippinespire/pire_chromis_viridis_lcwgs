#!/usr/bin/env bash
#SBATCH --job-name=cvi_fst_modhist
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --output=output/fst_historic_vs_modern/logs/%x_%j.out
#SBATCH --error=output/fst_historic_vs_modern/logs/%x_%j.err

set -euo pipefail

# Run from repository root regardless of submit location.
repo_root="/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs"
cd "$repo_root"

mkdir -p output/fst_historic_vs_modern
mkdir -p output/fst_historic_vs_modern/logs

module load container_env angsd

bash scripts/calc_fst_modern_historic.sh \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --realSFS-cmd "crun realSFS"
  --results-dir "${repo_root}/nf-pipelines/nf-angsd-diversity/results/angsd_pop"
  --historic "CviAPal_historic"
  --modern "CviCPal_modern"
  --outdir "${repo_root}/output/fst_historic_vs_modern"
  --folded
