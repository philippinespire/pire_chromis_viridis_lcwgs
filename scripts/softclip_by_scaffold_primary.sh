#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  softclip_by_scaffold_primary.sh <input.bam> [output_dir]

Description:
  Summarize soft clipping per scaffold using mapped primary alignments only.
  Excludes unmapped (0x4), secondary (0x100), QC-fail (0x200), duplicate (0x400), and supplementary (0x800) reads.

Output file name:
  softclip_by_scaffold_<input_bam_basename>.tsv

Examples:
  softclip_by_scaffold_primary.sh /path/to/sample.realn.bam
  softclip_by_scaffold_primary.sh /path/to/sample.realn.bam /path/to/output_dir
EOF
}

module load container_env
module load samtools/1.19

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage >&2
  exit 1
fi

bam_path="$1"
out_dir="${2:-.}"

if [[ ! -f "$bam_path" ]]; then
  echo "Error: BAM file not found: $bam_path" >&2
  exit 1
fi

if ! command -v crun.samtools >/dev/null 2>&1; then
  echo "Error: samtools is not available in PATH." >&2
  exit 1
fi

if ! command -v gawk >/dev/null 2>&1; then
  echo "Error: gawk is not available in PATH." >&2
  exit 1
fi

mkdir -p "$out_dir"

bam_base="$(basename "$bam_path")"
out_file="$out_dir/softclip_by_scaffold_${bam_base}.tsv"

# Exclude: 0x4 (unmapped), 0x100 (secondary), 0x200 (QC-fail), 0x400 (duplicate), 0x800 (supplementary)
{
  echo -e "scaffold\treads\treads_with_softclip\tfrac_reads_with_softclip\tquery_bases\tsoftclip_bases\tfrac_query_bases_softclipped"
  crun.samtools samtools view -F 3844 "$bam_path" | gawk '
{
  scaf=$3
  cig=$6
  if (scaf=="*" || cig=="*") next

  reads[scaf]++
  soft=0
  qbases=0

  while (match(cig, /([0-9]+)([MIDNSHP=X])/, m)) {
    len=m[1]+0
    op=m[2]
    if (op=="S") soft += len
    if (op ~ /[MIS=X]/) qbases += len
    cig=substr(cig, RSTART+RLENGTH)
  }

  if (soft > 0) reads_soft[scaf]++
  soft_bases[scaf] += soft
  query_bases[scaf] += qbases
}
END {
  for (s in reads) {
    fr = (reads[s] > 0) ? reads_soft[s] / reads[s] : 0
    fb = (query_bases[s] > 0) ? soft_bases[s] / query_bases[s] : 0
    print s, reads[s], reads_soft[s] + 0, fr, query_bases[s] + 0, soft_bases[s] + 0, fb
  }
}
' | sort -k1,1
} > "$out_file"

echo "Wrote: $out_file"
