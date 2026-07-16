#!/usr/bin/env bash

################################################################################
# List full paths of all BAM files in a directory, one per line
#
# USAGE:
#   bash list_bam_paths.sh <bam_dir> <output_file>
#
# ARGUMENTS:
#   $1: Path to directory containing BAM files (required)
#   $2: Path to output text file (required)
#
# OUTPUT:
#   A plain text file with one absolute BAM file path per line, sorted.
#
################################################################################

set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <bam_dir> <output_file>" >&2
    exit 1
fi

bam_dir="$1"
out_file="$2"

if [[ ! -d "${bam_dir}" ]]; then
    echo "ERROR: BAM directory not found: ${bam_dir}" >&2
    exit 1
fi

mkdir -p "$(dirname "${out_file}")"

find "${bam_dir}" -maxdepth 1 -type f -name '*.bam' | sort > "${out_file}"

echo "Wrote $(wc -l < "${out_file}") BAM paths to: ${out_file}"
