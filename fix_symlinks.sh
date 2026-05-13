#!/usr/bin/env bash
# reads 2026-05-12_corrupt_fastq_report.txt and either deletes or repoints symlinks for corrupt FASTQs
# workflow summary:
# 1) use the similarity TSV to delete or repoint individual corrupt symlinks
# 2) verify each R1/R2 pair resolves to the same source directory
# 3) if a pair is split across directories, repoint both mates to repaired files

set -euo pipefail

repo_root=/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs
tsv="$repo_root/2026-05-13_corrupt_symlink_target_similarity_check.tsv"
repaired_dir=/archive/carpenterlab/pire/pire_chromis_viridis_lcwgs/2nd_sequencing_run/fq_fp1_clmp_fp2_fqscrn_rprd

dry_run=0   # set to 0 to apply changes

# Map a symlink basename to the best matching repaired FASTQ in repaired_dir.
find_repaired_target() {
  local link_path="$1"
  local base grp num ex readn pattern

  base=$(basename "$link_path")
  if [[ "$base" =~ ^Cvi([AC]Pal)([0-9]{3})_Ex([0-9]+)_L[0-9]+_R([12])\.fastq\.gz$ ]]; then
    grp=${BASH_REMATCH[1]}
    num=${BASH_REMATCH[2]}
    ex=${BASH_REMATCH[3]}
    readn=${BASH_REMATCH[4]}
    pattern="Cvi-${grp}_${num}-Ex${ex}-*.R${readn}_*"
    find "$repaired_dir" -maxdepth 1 -type f -name "$pattern" | sort | head -n 1
  fi
}

awk -F '\t' 'NR>1 {print $1 "\t" $3 "\t" $4}' "$tsv" |
while IFS=$'\t' read -r rel sim_pattern match_field; do
  # TSV stores symlink paths relative to repo_root.
  link_path="$repo_root/$rel"

  # If no similar repaired file exists, delete the symlink.
  if [[ "$match_field" == "NO_SIMILAR_FILE" ]]; then
    # Only remove real symlinks so regular files are never deleted here.
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

# Ensure paired-end links from the same individual resolve to the same source directory.
# symlink_dir is inferred from the first data row in the TSV.
symlink_dir=$(awk -F '\t' 'NR==2 {print $1}' "$tsv" | xargs dirname)
symlink_dir="$repo_root/$symlink_dir"

if [[ -d "$symlink_dir" ]]; then
  declare -A r1_links
  declare -A r2_links

  # Index symlinks by sample ID prefix so mates can be compared together.
  while IFS= read -r link; do
    base=$(basename "$link")
    if [[ "$base" =~ ^(Cvi[AC]Pal[0-9]{3}_Ex[0-9]+_L[0-9]+)_R([12])\.fastq\.gz$ ]]; then
      id=${BASH_REMATCH[1]}
      readn=${BASH_REMATCH[2]}
      if [[ "$readn" == "1" ]]; then
        r1_links["$id"]="$link"
      else
        r2_links["$id"]="$link"
      fi
    fi
  done < <(find "$symlink_dir" -maxdepth 1 -type l -name 'Cvi*Pal*_Ex*_L*_R*.fastq.gz' | sort)

  for id in "${!r1_links[@]}"; do
    r1_link=${r1_links[$id]}
    r2_link=${r2_links[$id]:-}

    if [[ -z "$r2_link" ]]; then
      printf 'PAIR_MISSING_R2\t%s\n' "$r1_link" >&2
      continue
    fi

    r1_target=$(readlink -f "$r1_link" 2>/dev/null || true)
    r2_target=$(readlink -f "$r2_link" 2>/dev/null || true)

    if [[ -z "$r1_target" || -z "$r2_target" ]]; then
      printf 'PAIR_UNRESOLVED\t%s\t%s\n' "$r1_link" "$r2_link" >&2
      continue
    fi

    if [[ "$(dirname "$r1_target")" != "$(dirname "$r2_target")" ]]; then
      # When mates come from different directories, move both to repaired targets.
      new_r1=$(find_repaired_target "$r1_link")
      new_r2=$(find_repaired_target "$r2_link")

      if [[ -z "$new_r1" || -z "$new_r2" ]]; then
        printf 'PAIR_SKIP_NO_REPAIRED_TARGET\t%s\t%s\n' "$r1_link" "$r2_link" >&2
        continue
      fi

      if [[ "$dry_run" -eq 1 ]]; then
        printf 'DRY_RUN_PAIR1_UPDATE\t%s\t->\t%s\n' "$r1_link" "$new_r1"
        printf 'DRY_RUN_PAIR2_UPDATE\t%s\t->\t%s\n' "$r2_link" "$new_r2"
      else
        ln -sfn "$new_r1" "$r1_link"
        ln -sfn "$new_r2" "$r2_link"
        printf 'PAIR_UPDATED\t%s\t->\t%s\n' "$r1_link" "$new_r1"
        printf 'PAIR_UPDATED\t%s\t->\t%s\n' "$r2_link" "$new_r2"
      fi
    fi
  done
else
  printf 'SKIP_NO_SYMLINK_DIR\t%s\n' "$symlink_dir" >&2
fi
