#!/usr/bin/env bash

################################################################################
# Create a samplesheet.csv from BAM files in a directory
#
# DESCRIPTION:
#   Scans a directory for all .bam files and generates a CSV samplesheet with
#   one row per BAM file. Extracts sample, population, and era information from
#   the BAM filename.
#
# ARGUMENTS:
#   $1: Path to directory containing BAM files
#       (default: /archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-merged-unmerged/results/data/bam)
#   $2: Path to output CSV file (samplesheet)
#       (default: {script_dir}/../inputfiles/samplesheet.csv)
#
# FILENAME PATTERN:
#   Expected format: CviXPalNNNNN.*.bam
#   - Cvi: species prefix
#   - X: either A (historic) or C (modern)
#   - Pal: location
#   - NNNNN: numeric sample identifier
#
# OUTPUT COLUMNS:
#   sample: extracted sample ID (e.g., CviAPal001)
#   bam: absolute path to BAM file
#   pop: population code (CviAPal or CviCPal)
#   era: time period (historic for APal, modern for CPal)
#   region: always "Bali"
#
################################################################################

set -euo pipefail

# Default locations for this repository.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "${script_dir}/.." && pwd)"

bam_dir_default="/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-merged-unmerged/results/data/bam"
out_csv_default="${project_dir}/inputfiles/samplesheet.csv"

bam_dir="${1:-${bam_dir_default}}"
out_csv="${2:-${out_csv_default}}"

if [[ ! -d "${bam_dir}" ]]; then
    echo "ERROR: BAM directory not found: ${bam_dir}" >&2
    exit 1
fi

mkdir -p "$(dirname "${out_csv}")"

tmp_csv="$(mktemp)"
trap 'rm -f "${tmp_csv}"' EXIT

printf 'sample,bam,pop,era,region\n' > "${tmp_csv}"

while IFS= read -r bam_path; do
    bam_file="$(basename "${bam_path}")"

    if [[ "${bam_file}" =~ ^(Cvi(APal|CPal)[0-9]+).*\.bam$ ]]; then
        sample="${BASH_REMATCH[1]}"
        pop="Cvi${BASH_REMATCH[2]}"
    else
        echo "WARNING: Skipping unrecognized BAM filename pattern: ${bam_file}" >&2
        continue
    fi

    if [[ "${pop}" == "CviAPal" ]]; then
        era="historic"
    elif [[ "${pop}" == "CviCPal" ]]; then
        era="modern"
    else
        echo "WARNING: Skipping unsupported population in ${bam_file}" >&2
        continue
    fi

    printf '%s,%s,%s,%s,%s\n' "${sample}" "${bam_path}" "${pop}" "${era}" "Bali" >> "${tmp_csv}"
done < <(find "${bam_dir}" -maxdepth 1 -type f -name '*.bam' | sort)

mv "${tmp_csv}" "${out_csv}"

echo "Wrote samplesheet: ${out_csv}"