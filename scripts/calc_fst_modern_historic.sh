#!/usr/bin/env bash

# Malin Pinsky 2026
# This script was written with AI assistance using GitHub Copilot
# with the GPT-5.3-Codex model.

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  calc_fst_modern_historic.sh [options]

Compute weighted FST between historical and modern populations from ANGSD SAF outputs.

Options:
  --results-dir PATH   Directory containing ANGSD per-population SAF files
                       (default: nf-pipelines/nf-angsd-diversity/results/angsd_pop)
  --historic NAME      Historic SAF prefix (default: CviAPal_historic)
  --modern NAME        Modern SAF prefix (default: CviCPal_modern)
  --outdir PATH        Output directory (default: results/fst_historic_vs_modern)
  --threads INT        Number of threads for realSFS -P (default: 4)
  --folded             Use folded SFS/FST (default)
  --unfolded           Use unfolded SFS/FST
  --realSFS PATH       Path to realSFS executable (default: realSFS from PATH)
  --realSFS-cmd CMD    Full command used to invoke realSFS
                       (example: "crun realSFS")
  -h, --help           Show this help and exit

Example:
  bash scripts/calc_fst_modern_historic.sh --threads 12
EOF
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

results_dir="${repo_root}/nf-pipelines/nf-angsd-diversity/results/angsd_pop"
historic_prefix="CviAPal_historic"
modern_prefix="CviCPal_modern"
outdir="${repo_root}/output/fst_historic_vs_modern"
threads=4
folded=1
realsfs_cmd="realSFS"

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

read -r -a realsfs_parts <<< "$realsfs_cmd"
if [[ "${#realsfs_parts[@]}" -eq 0 ]]; then
  echo "Error: Empty realSFS command." >&2
  exit 1
fi

if ! command -v "${realsfs_parts[0]}" >/dev/null 2>&1; then
  echo "Error: realSFS launcher not found: '${realsfs_parts[0]}'" >&2
  echo "Load your ANGSD environment/module or pass --realSFS-cmd \"crun realSFS\"" >&2
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
run_realsfs "$historic_idx" "$modern_idx" -P "$threads" "${fold_args[@]}" > "$sfs_file"

echo "[2/4] Building FST index"
run_realsfs fst index "$historic_idx" "$modern_idx" -sfs "$sfs_file" -fstout "$fst_prefix" -P "$threads" "${fold_args[@]}"

echo "[3/4] Computing genome-wide FST stats"
run_realsfs fst stats "${fst_prefix}.fst.idx" > "$stats_file"

echo "[4/4] Computing sliding-window FST (50kb windows, 10kb step)"
run_realsfs fst stats2 "${fst_prefix}.fst.idx" -win 50000 -step 10000 > "$windows_file"

# realSFS fst stats is headerless in this workflow; the second numeric value is weighted FST.
weighted_fst="$(awk 'NF>=2 {print $2; exit}' "$stats_file")"

{
  echo "pair=${historic_prefix},${modern_prefix}"
  echo "folded=${folded}"
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
