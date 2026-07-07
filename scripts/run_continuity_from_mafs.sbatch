#!/usr/bin/env bash
#SBATCH --job-name=cvi_continuity
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Malin Pinsky 2026
# This script was written with AI assistance using GitHub Copilot
# with the GPT-5.3-Codex model.

set -euo pipefail

trap 'echo "Error: script failed at line ${LINENO} (exit=$?)" >&2' ERR

original_args=("$@")

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_continuity_from_mafs.sbatch [options]

Submit and run the continuity workflow in Slurm. This wrapper loads the needed
modules inside the batch job, then runs scripts/run_continuity_from_mafs.py.

Required:
  --samplesheet FILE        CSV with sample/bam/pop/era columns
  --mafs FILE               ANGSD .mafs or .mafs.gz file
  --output-prefix PREFIX    Output prefix for continuity outputs

Input selection:
  --bam-dir DIR             Directory containing BAM files to include
  --bam-list FILE           File with one BAM path per line

Continuity options:
  --sample-column NAME      Samplesheet sample column (default: sample)
  --bam-column NAME         Samplesheet BAM column (default: bam)
  --pop-column NAME         Samplesheet population column (default: pop)
  --era-column NAME         Samplesheet era column (default: era)
  --historic-era VALUE      Historic-era label (default: historic)
  --modern-era VALUE        Modern-era label (default: modern)
  --site-cutoff FLOAT       parse_reads_by_pop cutoff (default: 0.0)
  --min-cutoff FLOAT        coverage_filter lower percentile (default: 2.5)
  --max-cutoff FLOAT        coverage_filter upper percentile (default: 97.5)
  --num-core INT            Optimization cores (default: SLURM cpus-per-task)
  --keep-temp               Keep intermediate positions file
  --no-ref-polarization     Disable ref fallback when anc is absent
  --tmp-dir DIR             Temporary working directory parent
  --ancient-genotypes-path  Path to ancient_genotypes.py
  --samtools-cmd CMD        Samtools command passed to the Python script
                            (default in this Slurm wrapper: auto-resolved host samtools path)

Runtime options:
  --python-cmd CMD          Python command used inside the job
                            (default: /home/mpinsky/.conda/envs/AMBER/bin/python)

SLURM options:
  --job-name NAME           SLURM job name (default: cvi_continuity)
  --time HH:MM:SS           SLURM time limit (default: 96:00:00)
  --cpus INT                SLURM cpus-per-task (default: 8)
  --mem SIZE                SLURM memory (default: 32G)
  --partition NAME          SLURM partition/queue (optional)
  --log-dir DIR             SLURM log directory (default: <output-dir>/logs)
  -h, --help                Show this help and exit

Example:
  bash scripts/run_continuity_from_mafs.sbatch \
    --samplesheet nf-pipelines/nf-angsd-diversity-cvi-only/inputfiles/samplesheet.csv \
    --mafs nf-pipelines/nf-angsd-diversity-cvi-only/results/GL/Cvi.mafs.gz \
    --bam-dir nf-pipelines/nf-trim-merged-unmerged/results/data/bam \
    --output-prefix output/continuity/cvi-only
EOF
}

repo_root="/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs"
cd "$repo_root"

samplesheet=""
mafs=""
output_prefix=""
bam_dir=""
bam_list=""

sample_column="sample"
bam_column="bam"
pop_column="pop"
era_column="era"
historic_era="historic"
modern_era="modern"
site_cutoff="0.0"
min_cutoff="2.5"
max_cutoff="97.5"
num_core=""
keep_temp=0
no_ref_polarization=0
tmp_dir=""
ancient_genotypes_path="scripts/ancient_genotypes.py"
samtools_cmd="samtools"
samtools_cmd_user_set=0
python_cmd="/home/mpinsky/.conda/envs/AMBER/bin/python"
python_cmd_user_set=0

job_name="cvi_continuity"
time_limit="96:00:00"
cpus="8"
mem="32G"
partition=""
log_dir=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet)
      samplesheet="$2"
      shift 2
      ;;
    --mafs)
      mafs="$2"
      shift 2
      ;;
    --output-prefix)
      output_prefix="$2"
      shift 2
      ;;
    --bam-dir)
      bam_dir="$2"
      shift 2
      ;;
    --bam-list)
      bam_list="$2"
      shift 2
      ;;
    --sample-column)
      sample_column="$2"
      shift 2
      ;;
    --bam-column)
      bam_column="$2"
      shift 2
      ;;
    --pop-column)
      pop_column="$2"
      shift 2
      ;;
    --era-column)
      era_column="$2"
      shift 2
      ;;
    --historic-era)
      historic_era="$2"
      shift 2
      ;;
    --modern-era)
      modern_era="$2"
      shift 2
      ;;
    --site-cutoff)
      site_cutoff="$2"
      shift 2
      ;;
    --min-cutoff)
      min_cutoff="$2"
      shift 2
      ;;
    --max-cutoff)
      max_cutoff="$2"
      shift 2
      ;;
    --num-core)
      num_core="$2"
      shift 2
      ;;
    --keep-temp)
      keep_temp=1
      shift
      ;;
    --no-ref-polarization)
      no_ref_polarization=1
      shift
      ;;
    --tmp-dir)
      tmp_dir="$2"
      shift 2
      ;;
    --ancient-genotypes-path)
      ancient_genotypes_path="$2"
      shift 2
      ;;
    --samtools-cmd)
      samtools_cmd="$2"
      samtools_cmd_user_set=1
      shift 2
      ;;
    --python-cmd)
      python_cmd="$2"
      python_cmd_user_set=1
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

if [[ -z "$samplesheet" || -z "$mafs" || -z "$output_prefix" ]]; then
  echo "Error: --samplesheet, --mafs, and --output-prefix are required." >&2
  usage >&2
  exit 1
fi

if [[ ! -f "$samplesheet" ]]; then
  echo "Error: samplesheet not found: $samplesheet" >&2
  exit 1
fi

if [[ ! -f "$mafs" ]]; then
  echo "Error: mafs file not found: $mafs" >&2
  exit 1
fi

if [[ -n "$bam_dir" && ! -d "$bam_dir" ]]; then
  echo "Error: bam directory not found: $bam_dir" >&2
  exit 1
fi

if [[ -n "$bam_list" && ! -f "$bam_list" ]]; then
  echo "Error: bam list not found: $bam_list" >&2
  exit 1
fi

if [[ -z "$num_core" ]]; then
  num_core="$cpus"
fi

output_dir="$(dirname "$output_prefix")"
if [[ -z "$output_dir" || "$output_dir" == "." ]]; then
  output_dir="$repo_root"
fi
mkdir -p "$output_dir"

if [[ -z "$log_dir" ]]; then
  log_dir="$output_dir/logs"
fi
mkdir -p "$log_dir"

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  sbatch_args=(
    --job-name="$job_name"
    --time="$time_limit"
    --cpus-per-task="$cpus"
    --mem="$mem"
    --output="$log_dir/%x_%j.out"
    --error="$log_dir/%x_%j.err"
  )
  if [[ -n "$partition" ]]; then
    sbatch_args+=(--partition="$partition")
  fi
  sbatch "${sbatch_args[@]}" "$0" "${original_args[@]}"
  exit 0
fi

echo "[$(date -Is)] continuity slurm job starting"
echo "SLURM_JOB_ID=${SLURM_JOB_ID}"
echo "SLURM_JOB_NAME=${SLURM_JOB_NAME:-unknown}"
echo "HOST=$(hostname)"
echo "[$(date -Is)] checking module command availability"

# In some Slurm shells, the `module` function is not initialized by default.
if ! type module >/dev/null 2>&1; then
  if [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    set +e
    source /etc/profile.d/modules.sh
    source_status=$?
    set -e
    if [[ "$source_status" -ne 0 ]]; then
      echo "Warning: sourcing /etc/profile.d/modules.sh returned ${source_status}" >&2
    fi
  elif [[ -f /usr/share/lmod/lmod/init/bash ]]; then
    # shellcheck disable=SC1091
    set +e
    source /usr/share/lmod/lmod/init/bash
    source_status=$?
    set -e
    if [[ "$source_status" -ne 0 ]]; then
      echo "Warning: sourcing /usr/share/lmod/lmod/init/bash returned ${source_status}" >&2
    fi
  fi
fi

if ! type module >/dev/null 2>&1; then
  echo "Error: 'module' command is not available in this Slurm shell." >&2
  echo "Tip: run through a login shell or source your site module init script before load." >&2
  exit 1
fi

echo "[$(date -Is)] module command is available"

set +e
module load container_env python3 samtools
module_status=$?
set -e
if [[ "$module_status" -ne 0 ]]; then
  echo "Error: module load container_env python3 samtools failed with exit code ${module_status}" >&2
  exit "$module_status"
fi

if [[ "$samtools_cmd_user_set" -eq 0 ]]; then
  # Resolve a samtools command in the Slurm shell and pass it into Python.
  # Prefer a direct host samtools binary when available; otherwise fall back
  # to the cluster wrapper command.
 
  # If command -v succeeds, resolved_samtools gets the path. 
  # If it fails, it assigns an empty string, resolve_status is 0, and NO TRAP fires.
  resolved_samtools="$(command -v samtools 2>/dev/null)" || resolved_samtools=""
  
  if [[ -n "$resolved_samtools" ]]; then
    samtools_cmd="$resolved_samtools"
  else
    resolved_crun_samtools="$(command -v crun.samtools 2>/dev/null)" || resolved_crun_samtools=""
    
    if [[ -n "$resolved_crun_samtools" ]]; then
      samtools_cmd="$resolved_crun_samtools samtools"
    else
      echo "Warning: could not resolve samtools or crun.samtools; using fallback '$samtools_cmd'." >&2
    fi
  fi
fi

python_has_deps() {
  local cmd_text="$1"
  local -a cmd_parts
  # shellcheck disable=SC2206
  cmd_parts=($cmd_text)
  "${cmd_parts[@]}" -c 'import numpy, scipy, pandas, joblib, matplotlib' >/dev/null 2>&1
}

ensure_python_deps() {
  local cmd_text="$1"
  local -a cmd_parts
  # shellcheck disable=SC2206
  cmd_parts=($cmd_text)

  if python_has_deps "$cmd_text"; then
    return 0
  fi

  echo "Warning: selected python '$cmd_text' is missing required packages; attempting install." >&2
  "${cmd_parts[@]}" -m pip install --upgrade pip setuptools wheel >/dev/null 2>&1 || true
  "${cmd_parts[@]}" -m pip install numpy scipy pandas joblib matplotlib >/dev/null 2>&1 || true

  if ! python_has_deps "$cmd_text"; then
    echo "Error: python '$cmd_text' still missing required packages after install attempt." >&2
    echo "Tip: manually install numpy/scipy/pandas/joblib/matplotlib for this interpreter or pass --python-cmd." >&2
    return 1
  fi
  return 0
}

if [[ "$python_cmd_user_set" -eq 0 ]]; then
  need_host_python=0
  if [[ "$samtools_cmd" == *"crun.samtools"* ]]; then
    need_host_python=1
  fi

  if [[ "$need_host_python" -eq 1 && -x "$python_cmd" ]]; then
    :
  elif [[ -x /home/mpinsky/miniconda3/bin/python3 ]] && python_has_deps "/home/mpinsky/miniconda3/bin/python3"; then
    python_cmd="/home/mpinsky/miniconda3/bin/python3"
  elif python_has_deps "python3"; then
    python_cmd="python3"
  elif [[ "$need_host_python" -eq 0 ]] && command -v crun.python3 >/dev/null 2>&1 && python_has_deps "crun.python3 python3"; then
    python_cmd="crun.python3 python3"
  else
    if [[ "$need_host_python" -eq 1 ]]; then
      echo "Error: samtools command '$samtools_cmd' requires host Python, but no host Python with required packages was found." >&2
      echo "Tip: install numpy/scipy/pandas/joblib/matplotlib in '$python_cmd' or pass --python-cmd to a suitable host Python." >&2
    else
      echo "Error: could not find a Python command with required packages (numpy, scipy, pandas, joblib, matplotlib)." >&2
      echo "Tip: pass --python-cmd explicitly to a working Python environment." >&2
    fi
    exit 1
  fi
fi

if ! ensure_python_deps "$python_cmd"; then
  exit 1
fi

echo "python_cmd=$python_cmd"
echo "samtools_cmd=$samtools_cmd"
echo "PATH=$PATH"

# Ensure Python writes progress/errors promptly to Slurm logs.
export PYTHONUNBUFFERED=1

python_parts=($python_cmd)
python_call=("${python_parts[@]}" scripts/run_continuity_from_mafs.py)

python_call+=(--samplesheet "$samplesheet")
python_call+=(--mafs "$mafs")
python_call+=(--output-prefix "$output_prefix")
python_call+=(--sample-column "$sample_column")
python_call+=(--bam-column "$bam_column")
python_call+=(--pop-column "$pop_column")
python_call+=(--era-column "$era_column")
python_call+=(--historic-era "$historic_era")
python_call+=(--modern-era "$modern_era")
python_call+=(--site-cutoff "$site_cutoff")
python_call+=(--min-cutoff "$min_cutoff")
python_call+=(--max-cutoff "$max_cutoff")
python_call+=(--num-core "$num_core")
python_call+=(--ancient-genotypes-path "$ancient_genotypes_path")
python_call+=(--samtools-cmd "$samtools_cmd")

if [[ -n "$bam_dir" ]]; then
  python_call+=(--bam-dir "$bam_dir")
fi
if [[ -n "$bam_list" ]]; then
  python_call+=(--bam-list "$bam_list")
fi
if [[ -n "$tmp_dir" ]]; then
  python_call+=(--tmp-dir "$tmp_dir")
fi
if [[ "$keep_temp" -eq 1 ]]; then
  python_call+=(--keep-temp)
fi
if [[ "$no_ref_polarization" -eq 1 ]]; then
  python_call+=(--no-ref-polarization)
fi

echo "[$(date -Is)] running continuity command: ${python_call[*]}"
set +e
"${python_call[@]}"
exit_code=$?
set -e
echo "[$(date -Is)] continuity job finished with exit code ${exit_code}"
exit "$exit_code"