#!/usr/bin/env bash

# ==============================================================================
# SLURM Resource Allocation Headers
# Default parameters used when submitting this script via `sbatch`.
# Can be overridden on the CLI (e.g., `sbatch --time=48:00:00 ...`).
# ==============================================================================
#SBATCH --job-name=dystruct
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Enable strict Bash error handling:
#  -e: Exit immediately if any command returns a non-zero status.
#  -u: Treat unset variables as an error and exit immediately.
#  -o pipefail: Pipeline return status is the value of the last (rightmost) command
#               to exit with a non-zero status.
set -euo pipefail

# ==============================================================================
# Usage & Help Display
# Prints full documentation and CLI usage instructions.
# ==============================================================================
usage() {
  cat <<'EOF'
Usage:
  sbatch [slurm options] scripts/dystruct_slurm.sh [script options] [-- <extra dystruct args>]

Submit and run a DYSTRUCT job from an ANGSD beagle file.
Converts beagle genotype probabilities to an EIGENSTRAT .geno matrix
filtered by site list / .pos file, then runs DYSTRUCT.

Required:
  --beagle FILE               Input beagle(.gz) file
  --out-dir DIR               Output directory
  --npops INT                 Number of populations (K) passed to DYSTRUCT

Generation-time input (choose one):
  --generation-times FILE     One generation-time value per sample
  --generation-value INT      Use one constant generation-value for all samples

Optional Filtering:
  --sites-file, --pos-file FILE  Site list or ANGSD .pos file (chr pos, chr_pos, or chr:pos)
  --min-posterior FLOAT          If max posterior < threshold, write missing genotype 9 (default: 0)

Optional Parameters:
  --prefix STR                Base prefix for files (default: beagle basename)
  --seed INT                  Random seed for DYSTRUCT (default: 1145)
  --threads INT               OMP_NUM_THREADS for DYSTRUCT (default: npops)
  --dystruct-dir DIR          DYSTRUCT install directory
  --keep-intermediate         Keep generated .geno and sample-id files
  -h, --help                  Show this help and exit

Pass-through arguments:
  Any args after "--" are forwarded directly to DYSTRUCT.
EOF
}

# ==============================================================================
# Variable Initialization
# Set internal variable defaults before parsing command-line options.
# ==============================================================================
beagle=""
out_dir=""
npops=""
generation_times=""
generation_value=""
prefix=""
sites_file=""
seed="1145"
min_posterior="0"
threads=""
dystruct_dir="/archive/carpenterlab/pire/softwares/dystruct"
keep_intermediate=0
dystruct_extra_args=()

# ==============================================================================
# Helper Functions
# ==============================================================================

# Counts total lines in a plain text or gzipped file without uncompressing to disk.
count_lines() {
  local file="$1"
  if [[ "$file" == *.gz ]]; then
    gzip -cd -- "$file" | wc -l
  else
    wc -l < "$file"
  fi
}

# ==============================================================================
# Command-Line Argument Parsing
# Loop through provided positional flags using standard `case` matching.
# ==============================================================================
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
    --sites-file|--pos-file) # Supports both standard site-list and ANGSD .pos file flags
      sites_file="$2"
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
    --dystruct-arg) # Allows passing individual extra flags directly to DYSTRUCT
      dystruct_extra_args+=("$2")
      shift 2
      ;;
    --) # Anything following '--' is collected as pass-through parameters for DYSTRUCT
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

# ==============================================================================
# Input Validation Checks
# Ensure all mandatory parameters and files exist before running heavy jobs.
# ==============================================================================

# Check mandatory script arguments
if [[ -z "$beagle" || -z "$out_dir" || -z "$npops" ]]; then
  echo "Error: Missing required arguments." >&2
  usage >&2
  exit 1
fi

# Verify beagle file presence
if [[ ! -f "$beagle" ]]; then
  echo "Error: beagle file not found: $beagle" >&2
  exit 1
fi

# Verify site/position list file if provided
if [[ -n "$sites_file" && ! -f "$sites_file" ]]; then
  echo "Error: sites/pos file not found: $sites_file" >&2
  exit 1
fi

# Ensure generation-time inputs are mutually exclusive and at least one is provided
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

# ==============================================================================
# Output Directories & File Naming Setup
# ==============================================================================
mkdir -p "$out_dir"
log_dir="${out_dir}/logs"
mkdir -p "$log_dir"

# Automatically derive base output prefix from beagle filename if not specified
if [[ -z "$prefix" ]]; then
  prefix="$(basename "$beagle")"
  prefix="${prefix%.gz}"
  prefix="${prefix%.beagle}"
fi

# Standardized job output prefix appending population count K
run_prefix="${prefix}_K${npops}"

# Default OpenMP thread count to K (DYSTRUCT parallelization scales primarily with K)
if [[ -z "$threads" ]]; then
  threads="$npops"
fi

# Check for DYSTRUCT executable installation
dystruct_bin="${dystruct_dir}/bin/dystruct"
if [[ ! -x "$dystruct_bin" ]]; then
  echo "Error: DYSTRUCT binary not executable: $dystruct_bin" >&2
  exit 1
fi

# Select decompressed stream tool according to file extension
if [[ "$beagle" == *.gz ]]; then
  beagle_reader=(gzip -cd -- "$beagle")
else
  beagle_reader=(cat -- "$beagle")
fi

# ==============================================================================
# Beagle Header Parsing & Sample Count Inference
# Extracts sample IDs and counts total individuals ($n_{ind}$) from columns.
# Beagle layout: [marker, allele1, allele2, Ind0.g0, Ind0.g1, Ind0.g2, ...]
# ==============================================================================
sample_ids_file="${out_dir}/${run_prefix}.sample_ids.txt"
header_line=""

# Safely extract the first row without triggering pipefail issues with gzip
if [[ "$beagle" == *.gz ]]; then
  header_line="$(gzip -cd -- "$beagle" | sed -n '1p')"
else
  header_line="$(sed -n '1p' "$beagle")"
fi

# Parse column count and build sample ID list
header_info="$(printf '%s\n' "$header_line" | awk 'NR==1 {
  ncols=NF
  if (ncols < 6) { print "ERR: header has too few columns"; exit 2 }
  if (((ncols - 3) % 3) != 0) { print "ERR: expected 3 + 3*n_ind columns"; exit 3 }
  nind=(ncols - 3)/3
  for (i=1; i<=nind; i++) { print $(3 + (i-1)*3 + 1) }
  print "NIND=" nind > "/dev/stderr"
  exit
}' 2>"${out_dir}/${run_prefix}.header_info.tmp")"

if [[ -n "$header_info" && "$header_info" == ERR:* ]]; then
  echo "Error parsing beagle header: ${header_info#ERR: }" >&2
  rm -f "${out_dir}/${run_prefix}.header_info.tmp"
  exit 1
fi

# Save total individual count to variable
n_ind="$(awk -F'=' '/^NIND=/{print $2; exit}' "${out_dir}/${run_prefix}.header_info.tmp" | tr -d '[:space:]')"
rm -f "${out_dir}/${run_prefix}.header_info.tmp"

if [[ -z "$n_ind" ]]; then
  echo "Error: Failed to infer number of individuals from beagle header." >&2
  exit 1
fi

# Output extracted sample names
printf '%s\n' "$header_info" > "$sample_ids_file"

# ==============================================================================
# Beagle to EIGENSTRAT (.geno) Conversion
# Parses likelihood posteriors, selects max probability genotype (0, 1, or 2),
# filters by site/pos list, and writes dense genotype matrix format required by DYSTRUCT.
# ==============================================================================
geno_file="${out_dir}/${run_prefix}.geno"

echo "Converting beagle to EIGENSTRAT .geno"
echo "  Input:  $beagle"
echo "  Output: $geno_file"
echo "  n_ind:  $n_ind"
if [[ -n "$sites_file" ]]; then
  echo "  sites/pos filter: $sites_file"
fi

"${beagle_reader[@]}" | awk -v minp="$min_posterior" -v sites_file="$sites_file" '
# Helper function: Normalizes site IDs into standard format (chr:pos)
function canon_site_id(id, a) {
  if (match(id, /^(.*)[:_]([0-9]+)$/, a)) {
    return a[1] ":" a[2]
  }
  return id
}

BEGIN {
  use_filter = 0
  # If a site/pos file was passed, build an in-memory hash set of allowed loci
  if (sites_file != "") {
    use_filter = 1
    while ((getline line < sites_file) > 0) {
      gsub(/\r$/, "", line) # Strip Windows newline chars if present
      if (line == "" || line ~ /^#/) continue
      split(line, fields, /[ \t]+/)
      
      # Match single-column format (chr_pos or chr:pos)
      if (fields[1] != "") {
        keep[fields[1]] = 1
        keep[canon_site_id(fields[1])] = 1
      }
      # Match standard 2-column ANGSD .pos format (chr pos)
      if (fields[1] != "" && fields[2] != "") {
        keep[fields[1] ":" fields[2]] = 1
        keep[fields[1] "_" fields[2]] = 1
      }
    }
    close(sites_file)
  }
  kept_loci = 0
}

# Skip header row during genotype parsing
NR==1 {
  nind=(NF-3)/3
  next
}

{
  site_id = $1
  site_id_norm = canon_site_id(site_id)
  
  # Filter out loci not present in the allowed site set
  if (use_filter && !(site_id in keep) && !(site_id_norm in keep)) next

  row=""
  # Extract individual genotype triplet posteriors (P(RR), P(RA), P(AA))
  for (i=1; i<=nind; i++) {
    p0=$(3 + (i-1)*3 + 1) # Hom-ref
    p1=$(3 + (i-1)*3 + 2) # Het
    p2=$(3 + (i-1)*3 + 3) # Hom-alt

    g=0
    maxp=p0
    if (p1 > maxp) { maxp=p1; g=1 }
    if (p2 > maxp) { maxp=p2; g=2 }
    
    # If confidence is below --min-posterior threshold, mask as missing genotype (9)
    if (maxp < minp) { g=9 }

    row=row g
  }
  kept_loci++
  print row
}

END {
  if (kept_loci == 0) {
    if (use_filter) {
      print "Error: No loci from beagle matched entries in site list / .pos file." > "/dev/stderr"
    } else {
      print "Error: No loci found in beagle input." > "/dev/stderr"
    }
    exit 4
  }
}
' > "$geno_file"

# Calculate total retained loci count ($nloci$) required by DYSTRUCT arguments
nloci="$(wc -l < "$geno_file" | tr -d '[:space:]')"
if [[ -z "$nloci" || "$nloci" -le 0 ]]; then
  echo "Error: Converted .geno file is empty." >&2
  exit 1
fi

# ==============================================================================
# Generation-Times Construction & Verification
# ==============================================================================
# Generate constant generation value file if --generation-value was passed
if [[ -n "$generation_value" ]]; then
  generation_times="${out_dir}/${run_prefix}.generation_times.txt"
  awk -v n="$n_ind" -v g="$generation_value" 'BEGIN {for (i=1; i<=n; i++) print g}' > "$generation_times"
fi

# Verify row count in generation-times matches inferred individual count
gt_count="$(count_lines "$generation_times" | tr -d '[:space:]')"
if [[ "$gt_count" -ne "$n_ind" ]]; then
  echo "Error: Generation-times row count ($gt_count) does not match inferred n_ind ($n_ind)." >&2
  exit 1
fi

# ==============================================================================
# DYSTRUCT Execution
# Set OpenMP threads and pass formatted parameters to the binary executable.
# ==============================================================================
export OMP_NUM_THREADS="$threads"
dystruct_output_prefix="${out_dir}/${run_prefix}.dystruct"

# Construct command array
dystruct_cmd=(
  "$dystruct_bin"
  --input "$geno_file"
  --generation-times "$generation_times"
  --output "$dystruct_output_prefix"
  --npops "$npops"
  --nloci "$nloci"
  --seed "$seed"
)

# Append any pass-through arguments provided after '--'
if [[ ${#dystruct_extra_args[@]} -gt 0 ]]; then
  dystruct_cmd+=("${dystruct_extra_args[@]}")
fi

echo "Running DYSTRUCT"
echo "  OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "  npops=$npops nloci=$nloci n_ind=$n_ind"
echo "  Command: ${dystruct_cmd[*]}"

# Execute DYSTRUCT binary
"${dystruct_cmd[@]}"

# ==============================================================================
# Intermediate File Cleanup
# Remove converted temporary .geno file unless --keep-intermediate was specified.
# ==============================================================================
if [[ "$keep_intermediate" -eq 0 ]]; then
  rm -f "$geno_file"
  echo "Removed intermediate matrix: $geno_file"
else
  echo "Kept intermediate matrix: $geno_file"
fi