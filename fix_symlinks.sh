#!/usr/bin/env bash
# reads 2026-05-12_corrupt_fastq_report.txt and either deletes or repoints symlinks for corrupt FASTQs

set -euo pipefail

repo_root=/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs
tsv="$repo_root/2026-05-13_corrupt_symlink_target_similarity_check.tsv"
repaired_dir=/archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/fq_fp1_clmp_fp2_fqscrn_rprd

dry_run=0   # set to 0 to apply changes

awk -F '\t' 'NR>1 {print $1 "\t" $3 "\t" $4}' "$tsv" |
while IFS=$'\t' read -r rel sim_pattern match_field; do
  link_path="$repo_root/$rel"

  # If no similar repaired file exists, delete the symlink.
  if [[ "$match_field" == "NO_SIMILAR_FILE" ]]; then
    if [[ -L "$link_path" ]]; then
      if [[ "$dry_run" -eq 1 ]]; then
        printf 'DRY_RUN_DELETE\t%s\n' "$link_path"
      else
        rm -f "$link_path"
        printf 'DELETED\t%s\n' "$link_path"
      fi
    else
      printf 'SKIP_NOT_SYMLINK\t%s\n' "$link_path" >&2
    fi
    continue
  fi

  # Otherwise repoint symlink to repaired FASTQ.
  fq_pattern=$(printf '%s' "$sim_pattern" | sed -E 's/\.R([12])_\*$/\.R\1.fq.gz/')
  target=$(find "$repaired_dir" -maxdepth 1 -type f -name "$fq_pattern" | sort | head -n 1)

  if [[ -z "$target" ]]; then
    printf 'SKIP_NO_TARGET\t%s\t(pattern %s)\n' "$link_path" "$fq_pattern" >&2
    continue
  fi

  if [[ "$dry_run" -eq 1 ]]; then
    printf 'DRY_RUN_UPDATE\t%s\t->\t%s\n' "$link_path" "$target"
  else
    ln -sfn "$target" "$link_path"
    printf 'UPDATED\t%s\t->\t%s\n' "$link_path" "$target"
  fi
done
