#!/usr/bin/env bash
# Rebuilds nf-trim-merged-unmerged symlinks from raw FASTQs, falling back to repaired
# FASTQs when a pair contains at least one corrupted raw file.

set -euo pipefail

repo_root=/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs
raw_dir=/archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/fq_raw
repaired_dir=/archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/fq_fp1_clmp_fp2_fqscrn_rprd
symlink_dir="$repo_root/nf-pipelines/nf-trim-merged-unmerged/data/symlinks"

dry_run=1   # default: dry run mode. use --apply to actually create/update symlinks.

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --apply)
      dry_run=0
      shift
      ;;
    --help|-h)
      printf 'Usage: %s [OPTIONS]\n' "$(basename "$0")"
      printf '\nOptions:\n'
      printf '  --apply   Apply changes (create/update symlinks); default is dry-run preview\n'
      printf '  --help    Show this help message\n'
      exit 0
      ;;
    *)
      printf 'ERROR: Unknown option: %s\n' "$1" >&2
      printf 'Use --help for usage information.\n' >&2
      exit 1
      ;;
  esac
done

# Build or load the set of corrupt FASTQ files
declare -A corrupt_files

build_corrupt_set() {
  local report_file corrupt_report_path
  
  # Look for an existing corrupt report file (with various date patterns)
  for pattern in "$repo_root"/2026-*_corrupt_fastq_report.txt "$repo_root"/*_corrupt_fastq_report.txt; do
    if [[ -f "$pattern" ]]; then
      corrupt_report_path="$pattern"
      break
    fi
  done
  
  if [[ -n "$corrupt_report_path" ]]; then
    # Parse existing report: lines like "CORRUPT: /path/to/symlink_name.fastq.gz"
    # Extract symlink base names and convert to sample IDs for matching raw files
    while IFS= read -r line; do
      if [[ "$line" =~ CORRUPT:\ .*/([Cc]vi[AC][Pp]al[0-9]{3}_Ex[0-9]+_L[0-9]+_R[12])\.fastq\.gz$ ]]; then
        symlink_base=${BASH_REMATCH[1]}
        # Mark any raw file matching this sample as corrupt
        # We'll check against the actual raw file paths when needed
        corrupt_files["$symlink_base"]=1
      fi
    done < "$corrupt_report_path"
    printf 'INFO: Loaded corruption data from %s (%d corrupt symlinks)\n' "$corrupt_report_path" "${#corrupt_files[@]}" >&2
  else
    # Generate report by testing all raw FASTQs
    report_file="$repo_root/corrupt_fastq_report_$(date +%Y-%m-%d).txt"
    printf 'INFO: Generating corruption report to %s\n' "$report_file" >&2
    while IFS= read -r fastq; do
      if ! gzip -t "$fastq" >/dev/null 2>&1; then
        printf 'CORRUPT: %s\n' "$fastq" >> "$report_file"
        corrupt_files["$fastq"]=1
      fi
    done < <(find "$raw_dir" -maxdepth 1 -type f -name '*.fq.gz' | sort)
  fi
}

is_corrupt_fastq() {
  local fastq="$1"
  local base grp num ex readn symlink_name
  
  # Check if the full path is in the array (for newly generated reports)
  if [[ -v corrupt_files["$fastq"] ]]; then
    return 0
  fi
  
  # For symlink-based reports, extract the sample ID and check if that's marked corrupt
  base=$(basename "$fastq")
  if [[ "$base" =~ ^Cvi-([AC]Pal)_([0-9]{3})-Ex([0-9]+)-[^.]+-lcwgs-[0-9]+-[0-9]+\.([12])\.fq\.gz$ ]]; then
    grp=${BASH_REMATCH[1]}
    num=${BASH_REMATCH[2]}
    ex=${BASH_REMATCH[3]}
    readn=${BASH_REMATCH[4]}
    symlink_name="Cvi${grp}${num}_Ex${ex}_L4_R${readn}"
    if [[ -v corrupt_files["$symlink_name"] ]]; then
      return 0
    fi
  fi
  
  return 1
}

find_repaired_target() {
  local link_path="$1"
  local readn="$2"
  local base grp num ex pattern

  base=$(basename "$link_path")
  if [[ "$base" =~ ^Cvi([AC]Pal)([0-9]{3})_Ex([0-9]+)_L[0-9]+_R([12])\.fastq\.gz$ ]]; then
    grp=${BASH_REMATCH[1]}
    num=${BASH_REMATCH[2]}
    ex=${BASH_REMATCH[3]}
    pattern="Cvi-${grp}_${num}-Ex${ex}-*.clmp.fp2_repr.R${readn}.fq.gz"
    find "$repaired_dir" -maxdepth 1 -type f -name "$pattern" | sort | head -n 1
  fi
}

link_target() {
  local link_path="$1"
  local target="$2"

  if [[ -z "$target" ]]; then
    printf 'SKIP_NO_TARGET\t%s\n' "$link_path" >&2
    return 0
  fi

  if [[ "$dry_run" -eq 1 ]]; then
    printf 'DRY_RUN_UPDATE\t%s\t->\t%s\n' "$link_path" "$target"
  else
    ln -sfn "$target" "$link_path"
    printf 'UPDATED\t%s\t->\t%s\n' "$link_path" "$target"
  fi
}

if [[ ! -d "$symlink_dir" ]]; then
  printf 'SKIP_NO_SYMLINK_DIR\t%s\n' "$symlink_dir" >&2
  exit 0
fi

# Initialize corruption data set
build_corrupt_set

link_name_for_key() {
  local key="$1"
  local readn="$2"
  local grp num ex

  if [[ "$key" =~ ^Cvi-([AC]Pal)_([0-9]{3})-Ex([0-9]+)$ ]]; then
    grp=${BASH_REMATCH[1]}
    num=${BASH_REMATCH[2]}
    ex=${BASH_REMATCH[3]}
    printf '%s/Cvi%s%s_Ex%s_L4_R%s.fastq.gz' "$symlink_dir" "$grp" "$num" "$ex" "$readn"
  fi
}

declare -A raw_r1_paths
declare -A raw_r2_paths
declare -A sample_seen

while IFS= read -r raw_file; do
  base=$(basename "$raw_file")
  if [[ "$base" =~ ^Cvi-([AC]Pal)_([0-9]{3})-Ex([0-9]+)-[^.]+-lcwgs-[0-9]+-[0-9]+\.([12])\.fq\.gz$ ]]; then
    grp=${BASH_REMATCH[1]}
    num=${BASH_REMATCH[2]}
    ex=${BASH_REMATCH[3]}
    readn=${BASH_REMATCH[4]}
    key="Cvi-${grp}_${num}-Ex${ex}"
    sample_seen["$key"]=1
    if [[ "$readn" == "1" ]]; then
      raw_r1_paths["$key"]="$raw_file"
    else
      raw_r2_paths["$key"]="$raw_file"
    fi
  fi
done < <(find "$raw_dir" -maxdepth 1 -type f -name '*.fq.gz' | sort)

for key in "${!sample_seen[@]}"; do
  raw_r1=${raw_r1_paths[$key]:-}
  raw_r2=${raw_r2_paths[$key]:-}

  r1_link=$(link_name_for_key "$key" 1)
  r2_link=$(link_name_for_key "$key" 2)

  if [[ -z "$raw_r1" || -z "$raw_r2" || -z "$r1_link" || -z "$r2_link" ]]; then
    printf 'PAIR_SKIP_MISSING_RAW\t%s\n' "$key" >&2
    continue
  fi

  r1_corrupt=0
  r2_corrupt=0
  if is_corrupt_fastq "$raw_r1"; then
    r1_corrupt=1
  fi
  if is_corrupt_fastq "$raw_r2"; then
    r2_corrupt=1
  fi
  # Decide pair targets atomically so both mates stay in the same source directory.
  if [[ "$r1_corrupt" -eq 0 && "$r2_corrupt" -eq 0 ]]; then
    link_target "$r1_link" "$raw_r1"
    link_target "$r2_link" "$raw_r2"
    continue
  fi

  repaired_r1=$(find_repaired_target "$r1_link" 1)
  repaired_r2=$(find_repaired_target "$r2_link" 2)

  if [[ -n "$repaired_r1" && -n "$repaired_r2" ]]; then
    link_target "$r1_link" "$repaired_r1"
    link_target "$r2_link" "$repaired_r2"
  else
    [[ -z "$repaired_r1" && "$r1_corrupt" -eq 1 ]] && printf 'PAIR_SKIP_NO_REPAIRED_R1\t%s\n' "$r1_link" >&2
    [[ -z "$repaired_r2" && "$r2_corrupt" -eq 1 ]] && printf 'PAIR_SKIP_NO_REPAIRED_R2\t%s\n' "$r2_link" >&2
    printf 'PAIR_SKIP_SPLIT_SOURCE\t%s\n' "$key" >&2
    continue
  fi
done
