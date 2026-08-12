#!/usr/bin/env bash
#SBATCH --job-name=cvi_fst_modhist
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G

# Shell script to calculate Fst between modern and historic populations of Chromis viridis using ANGSD and realSFS.
# Supports Slurm self-submission and optional position filtering via --sites.
#
# Malin Pinsky 2026
# Written with AI assistance from GitHub Copilot with the GPT-5.3-Codex model
# and from Google Gemini

set -euo pipefail

# Preserve original command line arguments for Slurm self-submission
SAVED_ARGS=("$@")

# Safely determine repository root in both interactive shell and Slurm job contexts
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  repo_root="$SLURM_SUBMIT_DIR"
else
  script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  repo_root="$(cd "${script_dir}/.." && pwd)"
fi
cd "$repo_root"

# Function to safely expand paths to absolute paths relative to repo_root
to_abs_path() {
  local p="$1"
  if [[ -z "$p" ]]; then
    echo ""
  elif [[ "$p" = /* ]]; then
    echo "$p"
  else
    echo "${repo_root}/${p}"
  fi
}

usage() {
  cat <<'EOF'
Usage:
  bash calc_fst_modern_historic.sbatch [options]

Compute weighted FST between historical and modern populations from ANGSD SAF outputs.
If run outside of a Slurm job, this script automatically submits itself via sbatch.

Options:
  --results-dir PATH   Directory containing ANGSD per-population SAF files
                       (default: nf-pipelines/nf-angsd-diversity/results/angsd_pop)
  --historic NAME      Historic SAF prefix (default: CviAPal_historic)
  --modern NAME        Modern SAF prefix (default: CviCPal_modern)
  --sites PATH         Path to indexed sites file to restrict analysis (e.g., ld_pruned_sites.txt)
                       Requires accompanying .bin and .idx files generated via 'angsd sites index'
  --outdir PATH        Output directory (default: output/fst_historic_vs_modern)
  --threads INT        Number of threads for realSFS -P (default: SLURM_CPUS_PER_TASK or 12)
  --folded             Use folded SFS/FST (default)
  --unfolded           Use unfolded SFS/FST
  --realSFS PATH       Path to realSFS executable
  --realSFS-cmd CMD    Full command used to invoke realSFS (default: "crun realSFS")
  -h, --help           Show this help message and exit

Examples:
  # Run standard whole-genome FST calculation
  bash scripts/calc_fst_modern_historic.sbatch --threads 12

  # Restrict FST calculation to LD-pruned sites
  bash scripts/calc_fst_modern_historic.sbatch --sites data/ld_pruned_sites.txt
EOF
}

# Default values for arguments
results_dir="${repo_root}/nf-pipelines/nf-angsd-diversity/results/angsd_pop"
historic_prefix="CviAPal_historic"
modern_prefix="CviCPal_modern"
outdir="${repo_root}/output/fst_historic_vs_modern"
threads="${SLURM_CPUS_PER_TASK:-12}"
folded=1
realsfs_cmd="crun realSFS"
sites_file=""

# Argument parsing
while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-dir)
      results_dir="$(to_abs_path "$2")"
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
    --sites)
      sites_file="$(to_abs_path "$2")"
      shift 2
      ;;
    --outdir)
      outdir="$(to_abs_path "$2")"
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
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Error: Unknown argument '$1'" >&2
      usage >&2
      exit 1
      ;;
  esac
done

# Slurm self-submission wrapper
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  logdir="${outdir}/logs"
  mkdir -p "$outdir" "$logdir"

  sbatch \
    --output="${logdir}/%x_%j.out" \
    --error="${logdir}/%x_%j.err" \
    --cpus-per-task="$threads" \
    "$0" "${SAVED_ARGS[@]}"
  exit 0
fi

# Environment initialization
module load container_env angsd 2>/dev/null || true

read -r -a realsfs_parts <<< "$realsfs_cmd"
if [[ "${#realsfs_parts[@]}" -eq 0 ]]; then
  echo "Error: Empty realSFS command." >&2
  exit 1
fi

run_realsfs() {
  "${realsfs_parts[@]}" "$@"
}

historic_idx="${results_dir}/${historic_prefix}.saf.idx"
modern_idx="${results_dir}/${modern_prefix}.saf.idx"

if [[ ! -s "$historic_idx" ]]; then
  echo "Error: Missing or empty historic SAF index: $historic_idx" >&2
  exit 1
fi
if [[ ! -s "$modern_idx" ]]; then
  echo "Error: Missing or empty modern SAF index: $modern_idx" >&2
  exit 1
fi

# Validate sites file and index sidecars if supplied
sites_args=()
if [[ -n "$sites_file" ]]; then
  if [[ ! -s "$sites_file" ]]; then
    echo "Error: Specified sites file does not exist or is empty: $sites_file" >&2
    exit 1
  fi
  if [[ ! -s "${sites_file}.bin" || ! -s "${sites_file}.idx" ]]; then
    echo "Error: Missing ANGSD index files (${sites_file}.bin / ${sites_file}.idx)." >&2
    echo "Please index your sites file first using: crun angsd sites index $sites_file" >&2
    exit 1
  fi
  sites_args=("-sites" "$sites_file")
fi

mkdir -p "$outdir"

pair_tag="${historic_prefix}_vs_${modern_prefix}"
sfs_file="${outdir}/${pair_tag}.ml.sfs"
fst_prefix="${outdir}/${pair_tag}"
stats_file="${outdir}/${pair_tag}.fst.stats.txt"
windows_file="${outdir}/${pair_tag}.fst.win50kb.step10kb.txt"
summary_file="${outdir}/${pair_tag}.weighted_fst.txt"

fold_args=()
if [[ "$folded" -eq 1 ]]; then
  fold_args+=("-fold" "1")
fi

echo "[1/4] Estimating 2D SFS"
run_realsfs "$historic_idx" "$modern_idx" -P "$threads" "${fold_args[@]}" "${sites_args[@]}" > "$sfs_file"

echo "[2/4] Building FST index"
run_realsfs fst index "$historic_idx" "$modern_idx" -sfs "$sfs_file" -fstout "$fst_prefix" -P "$threads" "${fold_args[@]}" "${sites_args[@]}"

echo "[3/4] Computing genome-wide FST stats"
run_realsfs fst stats "${fst_prefix}.fst.idx" > "$stats_file"

echo "[4/4] Computing sliding-window FST (50kb windows, 10kb step)"
run_realsfs fst stats2 "${fst_prefix}.fst.idx" -win 50000 -step 10000 > "$windows_file"

weighted_fst="$(awk 'NF>=2 {print $2; exit}' "$stats_file")"

{
  echo "pair=${historic_prefix},${modern_prefix}"
  echo "folded=${folded}"
  if [[ -n "$sites_file" ]]; then
    echo "sites_file=${sites_file}"
  fi
  if [[ -n "$weighted_fst" ]]; then
    echo "weighted_fst=${weighted_fst}"
  else
    echo "weighted_fst=NA"
    echo "warning=Could_not_parse_weighted_FST_from_stats_file"
  fi
  echo "sfs_file=${sfs_file}"
  echo "fst_index=${fst_prefix}.fst.idx"
  echo "stats_file=${stats_file}"
  echo "windows_file=${windows_file}"
} > "$summary_file"

echo "Done. Summary written to: $summary_file"
if [[ -n "$weighted_fst" ]]; then
  echo "Weighted FST (${historic_prefix} vs ${modern_prefix}) = $weighted_fst"
else
  echo "Weighted FST could not be parsed automatically. Inspect: $stats_file"
fi