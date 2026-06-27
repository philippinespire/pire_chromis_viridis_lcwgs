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
  bash scripts/dystruct_slurm.sh [options] [-- <extra dystruct args>]

Submit and run a DYSTRUCT SLURM job from an ANGSD beagle file.
This script converts beagle genotype probabilities to an EIGENSTRAT .geno matrix
by taking the most likely genotype per sample/locus (optionally masked to 9 when
posterior confidence is below a threshold), then runs DYSTRUCT.

Required:
  --beagle FILE               Input beagle(.gz) file (e.g., Cvi.beagle.gz)
  --out-dir DIR               Output directory (e.g., output/dystruct)
  --npops INT                 Number of populations (K) passed to DYSTRUCT

Generation-time input (choose one):
  --generation-times FILE     One generation-time value per sample (required by DYSTRUCT)
  --generation-value INT      Use one constant generation-time value for all samples

Optional:
  --prefix STR                Prefix for intermediate/output files (default: beagle basename)
  --seed INT                  Random seed for DYSTRUCT (default: 1145)
  --min-posterior FLOAT       If max posterior < threshold, write missing genotype 9 (default: 0)
  --threads INT               OMP_NUM_THREADS for DYSTRUCT (default: npops)
  --dystruct-dir DIR          DYSTRUCT install directory (default: /archive/carpenterlab/pire/softwares/dystruct)
  --keep-intermediate         Keep generated .geno and sample-id files
  --dystruct-arg ARG          Extra token to pass to DYSTRUCT (repeatable)

SLURM options:
  --job-name NAME             SLURM job name (default: dystruct)
  --time HH:MM:SS             SLURM time (default: 24:00:00)
  --cpus INT                  SLURM cpus-per-task (default: 8)
  --mem SIZE                  SLURM memory (default: 32G)
  --partition NAME            SLURM partition/queue (optional)
  --log-dir DIR               SLURM log directory (default: <out-dir>/logs)
  -h, --help                  Show this help and exit

Pass-through arguments:
  Any args after "--" are forwarded directly to DYSTRUCT.
  Example: -- --epochs 100 --hold-out-fraction 0.1 --hold-out-seed 55307

Examples:
  bash scripts/dystruct_slurm.sh \
    --beagle nf-pipelines/nf-angsd-diversity-cvi-only/results/GL/Cvi.beagle.gz \
    --out-dir output/dystruct \
    --npops 3 \
    --generation-times output/dystruct/Cvi.generation_times.txt \
    --seed 1145

  bash scripts/dystruct_slurm.sh \
    --beagle nf-pipelines/nf-angsd-diversity-cvi-only/results/GL/Cvi.beagle.gz \
    --out-dir output/dystruct \
    --npops 5 \
    --generation-value 0 \
    --min-posterior 0.8 \
    -- --epochs 100 --hold-out-fraction 0.1
EOF
}

beagle=""
out_dir=""
npops=""
generation_times=""
generation_value=""
prefix=""
seed="1145"
min_posterior="0"
threads=""
dystruct_dir="/archive/carpenterlab/pire/softwares/dystruct"
keep_intermediate=0
dystruct_extra_args=()

job_name="dystruct"
time_limit="24:00:00"
cpus=8
mem="32G"
partition=""
log_dir=""

count_lines() {
  local file="$1"
  if [[ "$file" == *.gz ]]; then
    gzip -cd -- "$file" | wc -l
  else
    wc -l < "$file"
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --beagle)
      beagle="$2"
      shift 2
      ;;
    --out-dir)
      out_dir="$2"
      shift 2
      ;;
    --npops)
      npops="$2"
      shift 2
      ;;
    --generation-times)
      generation_times="$2"
      shift 2
      ;;
    --generation-value)
      generation_value="$2"
      shift 2
      ;;
    --prefix)
      prefix="$2"
      shift 2
      ;;
    --seed)
      seed="$2"
      shift 2
      ;;
    --min-posterior)
      min_posterior="$2"
      shift 2
      ;;
    --threads)
      threads="$2"
      shift 2
      ;;
    --dystruct-dir)
      dystruct_dir="$2"
      shift 2
      ;;
    --keep-intermediate)
      keep_intermediate=1
      shift
      ;;
    --dystruct-arg)
      dystruct_extra_args+=("$2")
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
    --)
      shift
      dystruct_extra_args+=("$@")
      break
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

if [[ -z "$beagle" || -z "$out_dir" || -z "$npops" ]]; then
  echo "Error: Missing required arguments." >&2
  usage >&2
  exit 1
fi

if [[ ! -f "$beagle" ]]; then
  echo "Error: beagle file not found: $beagle" >&2
  exit 1
fi

if [[ -n "$generation_times" && -n "$generation_value" ]]; then
  echo "Error: Use either --generation-times or --generation-value, not both." >&2
  exit 1
fi

if [[ -z "$generation_times" && -z "$generation_value" ]]; then
  echo "Error: Provide one of --generation-times or --generation-value." >&2
  exit 1
fi

if [[ -n "$generation_times" && ! -f "$generation_times" ]]; then
  echo "Error: generation-times file not found: $generation_times" >&2
  exit 1
fi

mkdir -p "$out_dir"
if [[ -z "$log_dir" ]]; then
  log_dir="${out_dir}/logs"
fi
mkdir -p "$log_dir"

if [[ -z "$prefix" ]]; then
  prefix="$(basename "$beagle")"
  prefix="${prefix%.gz}"
  prefix="${prefix%.beagle}"
fi

if [[ -z "$threads" ]]; then
  threads="$npops"
fi

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

dystruct_bin="${dystruct_dir}/bin/dystruct"
if [[ ! -x "$dystruct_bin" ]]; then
  echo "Error: DYSTRUCT binary not executable: $dystruct_bin" >&2
  echo "Hint: build DYSTRUCT first (cd $dystruct_dir && make)." >&2
  exit 1
fi

if [[ "$beagle" == *.gz ]]; then
  beagle_reader=(gzip -cd -- "$beagle")
else
  beagle_reader=(cat -- "$beagle")
fi

# Infer n_ind from beagle header and create sample-id list.
sample_ids_file="${out_dir}/${prefix}.sample_ids.txt"
header_info="$(${beagle_reader[@]} | awk 'NR==1 {
  ncols=NF
  if (ncols < 6) {
    print "ERR: header has too few columns"; exit 2
  }
  if (((ncols - 3) % 3) != 0) {
    print "ERR: expected 3 + 3*n_ind columns"; exit 3
  }
  nind=(ncols - 3)/3
  for (i=1; i<=nind; i++) {
    print $(3 + (i-1)*3 + 1)
  }
  print "NIND=" nind > "/dev/stderr"
  exit
}' 2>"${out_dir}/${prefix}.header_info.tmp")"

if [[ -n "$header_info" && "$header_info" == ERR:* ]]; then
  echo "Error parsing beagle header: ${header_info#ERR: }" >&2
  rm -f "${out_dir}/${prefix}.header_info.tmp"
  exit 1
fi

if [[ ! -s "${out_dir}/${prefix}.header_info.tmp" ]]; then
  echo "Error: Could not parse beagle header from $beagle" >&2
  exit 1
fi

n_ind="$(awk -F'=' '/^NIND=/{print $2; exit}' "${out_dir}/${prefix}.header_info.tmp" | tr -d '[:space:]')"
rm -f "${out_dir}/${prefix}.header_info.tmp"

if [[ -z "$n_ind" ]]; then
  echo "Error: Failed to infer number of individuals from beagle header." >&2
  exit 1
fi

printf '%s\n' "$header_info" > "$sample_ids_file"

geno_file="${out_dir}/${prefix}.geno"

echo "Converting beagle to EIGENSTRAT .geno"
echo "  Input:  $beagle"
echo "  Output: $geno_file"
echo "  n_ind:  $n_ind"
echo "  min_posterior threshold: $min_posterior"

# Convert per-locus genotype posteriors to hard genotype calls in EIGENSTRAT format.
${beagle_reader[@]} | awk -v minp="$min_posterior" '
NR==1 {
  nind=(NF-3)/3
  next
}
{
  row=""
  for (i=1; i<=nind; i++) {
    p0=$(3 + (i-1)*3 + 1)
    p1=$(3 + (i-1)*3 + 2)
    p2=$(3 + (i-1)*3 + 3)

    g=0
    maxp=p0
    if (p1 > maxp) { maxp=p1; g=1 }
    if (p2 > maxp) { maxp=p2; g=2 }
    if (maxp < minp) { g=9 }

    row=row g
  }
  print row
}
' > "$geno_file"

nloci="$(wc -l < "$geno_file" | tr -d '[:space:]')"
if [[ -z "$nloci" || "$nloci" -le 0 ]]; then
  echo "Error: Converted .geno file is empty: $geno_file" >&2
  exit 1
fi

if [[ -n "$generation_value" ]]; then
  generation_times="${out_dir}/${prefix}.generation_times.txt"
  awk -v n="$n_ind" -v g="$generation_value" 'BEGIN {for (i=1; i<=n; i++) print g}' > "$generation_times"
fi

gt_count="$(count_lines "$generation_times" | tr -d '[:space:]')"
if [[ "$gt_count" -ne "$n_ind" ]]; then
  echo "Error: Generation-times row count ($gt_count) does not match inferred n_ind ($n_ind)." >&2
  echo "  File: $generation_times" >&2
  exit 1
fi

export OMP_NUM_THREADS="$threads"

dystruct_output_prefix="${out_dir}/${prefix}.dystruct"

dystruct_cmd=(
  "$dystruct_bin"
  --input "$geno_file"
  --generation-times "$generation_times"
  --output "$dystruct_output_prefix"
  --npops "$npops"
  --nloci "$nloci"
  --seed "$seed"
)

if [[ ${#dystruct_extra_args[@]} -gt 0 ]]; then
  dystruct_cmd+=("${dystruct_extra_args[@]}")
fi

echo "Running DYSTRUCT"
echo "  OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "  npops=$npops nloci=$nloci n_ind=$n_ind"
echo "  generation_times=$generation_times"
echo "  output_prefix=$dystruct_output_prefix"
echo "  sample_ids=$sample_ids_file"
echo "Command: ${dystruct_cmd[*]}"

"${dystruct_cmd[@]}"

echo "DYSTRUCT completed."
echo "Primary outputs:"
echo "  ${dystruct_output_prefix}.theta"
echo "  ${dystruct_output_prefix}.freqs"
echo "Sample IDs used (column order from beagle):"
echo "  ${sample_ids_file}"

if [[ "$keep_intermediate" -eq 0 ]]; then
  rm -f "$geno_file"
  echo "Removed intermediate: $geno_file"
else
  echo "Kept intermediate: $geno_file"
fi
