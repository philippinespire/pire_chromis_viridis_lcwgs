#!/usr/bin/env bash

# Malin Pinsky 2026
# This script was written with AI assistance using GitHub Copilot
# with the GPT-5.3-Codex model.

set -euo pipefail

# Preserve original CLI arguments for self-submission via sbatch.
original_args=("$@")

usage() {
  cat <<'EOF'
Usage:
  bash scripts/ld_prune_slurm.sh [options]

Submit and run an ngsLD SLURM job for LD pruning.

Required:
  --probs FILE         Input genotype probabilities file for ngsLD (--geno + --probs)
  --pos FILE           Position file for ngsLD
  --out FILE           Output LD table from ngsLD
  --prune-out FILE     Output file with pruned/unlinked sites from prune_graph
  --max-kb-dist INT    Maximum pairwise distance in kb (ngsLD + pruning filter)
  --min-weight FLOAT   Minimum r2 threshold for pruning filter

Optional:
  --n-ind INT          Number of individuals (auto-inferred from --probs if omitted)
  --n-sites INT        Number of sites (auto-inferred from --pos if omitted)
  --n-threads INT      Threads passed to ngsLD (default: SLURM cpus-per-task)
  --no-prune           Run ngsLD only (skip prune_graph)
  --extra-arg ARG      Extra ngsLD argument token (repeatable)
  --prune-extra-arg A  Extra prune_graph argument token (repeatable)
  --job-name NAME      SLURM job name (default: ngsLD_prune)
  --time HH:MM:SS      SLURM time (default: 12:00:00)
  --cpus INT           SLURM cpus-per-task (default: 8)
  --mem SIZE           SLURM memory (default: 32G)
  --partition NAME     SLURM partition/queue (optional)
  --log-dir DIR        SLURM log directory (default: <out_dir>/logs)
  -h, --help           Show this help and exit

Example:
  bash scripts/ld_prune_slurm.sh \
    --probs output/angsd/all.beagle.gz \
    --pos output/angsd/all.pos.gz \
    --out output/ngsld/all.ld \
    --prune-out output/ngsld/all.unlinked.pos \
    --max-kb-dist 50 \
    --min-weight 0.2 \
    --n-threads 16
EOF
}

probs=""
pos=""
n_ind=""
n_sites=""
out=""
prune_out=""
max_kb_dist=""
min_weight=""
n_threads=""
run_prune=1
extra_args=()
prune_extra_args=()

job_name="ngsLD_prune"
time_limit="12:00:00"
cpus=8
mem="32G"
partition=""
log_dir=""

# Count lines from either plain text or gzipped files.
line_count() {
  local file="$1"
  if [[ "$file" == *.gz ]]; then
    gzip -cd -- "$file" | wc -l
  else
    wc -l < "$file"
  fi
}

# Count fields from first line in either plain text or gzipped files.
first_line_field_count() {
  local file="$1"
  if [[ "$file" == *.gz ]]; then
    # Avoid SIGPIPE under `set -o pipefail` by fully consuming gzip stream.
    gzip -cd -- "$file" | sed -n '1p' | awk '{print NF}'
  else
    awk 'NR==1 {print NF; exit}' "$file"
  fi
}

# Parse command-line options.
while [[ $# -gt 0 ]]; do
  case "$1" in
    --probs)
      probs="$2"
      shift 2
      ;;
    --pos)
      pos="$2"
      shift 2
      ;;
    --n-ind)
      n_ind="$2"
      shift 2
      ;;
    --n-sites)
      n_sites="$2"
      shift 2
      ;;
    --out)
      out="$2"
      shift 2
      ;;
    --prune-out)
      prune_out="$2"
      shift 2
      ;;
    --max-kb-dist)
      max_kb_dist="$2"
      shift 2
      ;;
    --min-weight)
      min_weight="$2"
      shift 2
      ;;
    --n-threads)
      n_threads="$2"
      shift 2
      ;;
    --no-prune)
      run_prune=0
      shift
      ;;
    --extra-arg)
      extra_args+=("$2")
      shift 2
      ;;
    --prune-extra-arg)
      prune_extra_args+=("$2")
      shift 2
      ;;
    --job-name)
      job_name="$2"
      shift 2
      ;;
    --time)
      time_limit="$2"
      shift 2
      ;;
    --cpus)
      cpus="$2"
      shift 2
      ;;
    --mem)
      mem="$2"
      shift 2
      ;;
    --partition)
      partition="$2"
      shift 2
      ;;
    --log-dir)
      log_dir="$2"
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

# Validate required user inputs.
if [[ -z "$probs" || -z "$pos" || -z "$out" || -z "$max_kb_dist" || -z "$min_weight" ]]; then
  echo "Error: Missing required arguments." >&2
  usage >&2
  exit 1
fi

if [[ "$run_prune" -eq 1 && -z "$prune_out" ]]; then
  echo "Error: --prune-out is required unless --no-prune is used." >&2
  usage >&2
  exit 1
fi

if [[ ! -f "$probs" ]]; then
  echo "Error: Input probs file not found: $probs" >&2
  exit 1
fi
if [[ ! -f "$pos" ]]; then
  echo "Error: Input pos file not found: $pos" >&2
  exit 1
fi

# Infer dataset dimensions from the provided files.
beagle_ncols="$(first_line_field_count "$probs")"
if [[ -z "$beagle_ncols" || "$beagle_ncols" -lt 6 ]]; then
  echo "Error: Could not parse a valid beagle header from: $probs" >&2
  exit 1
fi

if [[ $(((beagle_ncols - 3) % 3)) -ne 0 ]]; then
  echo "Error: Beagle file appears malformed (NF=${beagle_ncols}); expected 3 + 3*n_ind columns." >&2
  exit 1
fi

inferred_n_ind="$(((beagle_ncols - 3) / 3))"
inferred_n_sites_from_pos="$(line_count "$pos" | tr -d '[:space:]')"
beagle_total_lines="$(line_count "$probs" | tr -d '[:space:]')"
inferred_n_sites_from_beagle="$((beagle_total_lines - 1))"

# Ensure beagle and position inputs describe the same set/order of sites.
if [[ "$inferred_n_sites_from_beagle" -ne "$inferred_n_sites_from_pos" ]]; then
  echo "Error: Site count mismatch between --probs and --pos." >&2
  echo "  Beagle data rows (excluding header): $inferred_n_sites_from_beagle" >&2
  echo "  Position rows: $inferred_n_sites_from_pos" >&2
  exit 1
fi

if [[ -z "$n_ind" ]]; then
  n_ind="$inferred_n_ind"
fi

if [[ -z "$n_sites" ]]; then
  n_sites="$inferred_n_sites_from_pos"
fi

if [[ "$n_ind" -ne "$inferred_n_ind" ]]; then
  echo "Error: --n-ind ($n_ind) does not match inferred n_ind from --probs ($inferred_n_ind)." >&2
  exit 1
fi

if [[ "$n_sites" -ne "$inferred_n_sites_from_pos" ]]; then
  echo "Error: --n-sites ($n_sites) does not match inferred n_sites from --pos ($inferred_n_sites_from_pos)." >&2
  exit 1
fi

# Prepare output and log directories.
out_dir="$(dirname "$out")"
mkdir -p "$out_dir"

if [[ -n "$prune_out" ]]; then
  prune_out_dir="$(dirname "$prune_out")"
  mkdir -p "$prune_out_dir"
fi

if [[ -z "$log_dir" ]]; then
  log_dir="${out_dir}/logs"
fi
mkdir -p "$log_dir"

# If not already on a SLURM node, submit this script as a batch job.
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  if ! command -v sbatch >/dev/null 2>&1; then
    echo "Error: sbatch not found in PATH." >&2
    exit 1
  fi

  sbatch_cmd=(
    sbatch
    --job-name="$job_name"
    --time="$time_limit"
    --cpus-per-task="$cpus"
    --mem="$mem"
    --output="${log_dir}/%x_%j.out"
    --error="${log_dir}/%x_%j.err"
  )

  if [[ -n "$partition" ]]; then
    sbatch_cmd+=(--partition="$partition")
  fi

  sbatch_cmd+=("$0" "${original_args[@]}")

  submit_output="$("${sbatch_cmd[@]}")"
  echo "$submit_output"
  exit 0
fi

# Runtime environment for ngsLD/prune_graph.
module load container_env ngsTools

threads="${n_threads:-${SLURM_CPUS_PER_TASK:-$cpus}}"

# Report key resolved parameters for reproducibility and debugging.
echo "Resolved dataset dimensions: n_ind=$n_ind n_sites=$n_sites"
echo "Inputs: probs=$probs pos=$pos"
echo "Outputs: ngsld=$out$(if [[ "$run_prune" -eq 1 ]]; then printf ' prune=%s' "$prune_out"; fi)"
echo "Resources: threads=$threads"

# Build and run ngsLD command.
ngsld_cmd=(
  crun ngsLD
  --geno "$probs"
  --probs
  --pos "$pos"
  --n_ind "$n_ind"
  --n_sites "$n_sites"
  --max_kb_dist "$max_kb_dist"
  --n_threads "$threads"
  --out "$out"
)

if [[ ${#extra_args[@]} -gt 0 ]]; then
  ngsld_cmd+=("${extra_args[@]}")
fi

echo "Running: ${ngsld_cmd[*]}"
"${ngsld_cmd[@]}"

echo "Done. ngsLD output: $out"

# Optionally prune linked sites using prune_graph.
if [[ "$run_prune" -eq 1 ]]; then
  prune_launcher=(prune_graph)
  if ! command -v prune_graph >/dev/null 2>&1; then
    if command -v crun >/dev/null 2>&1; then
      prune_launcher=(crun prune_graph)
    else
      echo "Error: prune_graph not found (neither direct nor via crun). Install/enable prune_graph or rerun with --no-prune." >&2
      exit 1
    fi
  fi

  max_bp_dist="$(awk -v kb="$max_kb_dist" 'BEGIN{printf "%.6f", kb*1000}')"

  # Normalize ngsLD output to a stable edge table for prune_graph:
  # site1 site2 dist r2
  # This avoids parser/column-name issues across ngsLD/prune_graph versions.
  prune_edges_tsv="${out}.prune_graph_input.tsv"
  first_dist_field="$(awk 'NR==1 {print $7; exit}' "$out")"
  if [[ "$first_dist_field" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    awk 'BEGIN{OFS="\t"; print "site1","site2","dist","r2"} {print $1,$4,$7,$8}' "$out" > "$prune_edges_tsv"
  else
    awk 'BEGIN{OFS="\t"; print "site1","site2","dist","r2"} NR>1 {print $1,$4,$7,$8}' "$out" > "$prune_edges_tsv"
  fi
  weight_filter="dist <= ${max_bp_dist} && r2 >= ${min_weight}"

  prune_cmd=(
    "${prune_launcher[@]}"
    --header
    --in "$prune_edges_tsv"
    --weight-field "r2"
    --weight-filter "$weight_filter"
    --out "$prune_out"
  )

  if [[ ${#prune_extra_args[@]} -gt 0 ]]; then
    prune_cmd+=("${prune_extra_args[@]}")
  fi

  echo "Running: ${prune_cmd[*]}"
  "${prune_cmd[@]}"

  rm -f "$prune_edges_tsv"

  echo "Done. Pruned/unlinked sites: $prune_out"
fi
