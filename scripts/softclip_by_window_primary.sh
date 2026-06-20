#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  softclip_by_window_primary.sh <input.bam> <window_size_bp> [output_dir]

Description:
  Summarize soft clipping in non-overlapping windows per scaffold using mapped primary alignments only.
  Excludes unmapped (0x4), secondary (0x100), QC-fail (0x200), duplicate (0x400), and supplementary (0x800) reads.
  Each alignment is assigned to the window containing its leftmost mapped position (POS).

Output file name:
  softclip_by_scaffold_window_<window_size_bp>_<input_bam_basename>.tsv

Examples:
  softclip_by_window_primary.sh /path/to/sample.realn.bam 100000
  softclip_by_window_primary.sh /path/to/sample.realn.bam 100000 /path/to/output_dir
EOF
}

module load container_env
module load samtools/1.19

if [[ ${1:-} == "-h" || ${1:-} == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 2 || $# -gt 3 ]]; then
  usage >&2
  exit 1
fi

bam_path="$1"
window_size="$2"
out_dir="${3:-.}"

if [[ ! -f "$bam_path" ]]; then
  echo "Error: BAM file not found: $bam_path" >&2
  exit 1
fi

if [[ ! "$window_size" =~ ^[1-9][0-9]*$ ]]; then
  echo "Error: window_size_bp must be a positive integer." >&2
  exit 1
fi

if command -v crun.samtools >/dev/null 2>&1; then
  samtools_cmd=(crun.samtools samtools)
elif command -v samtools >/dev/null 2>&1; then
  samtools_cmd=(samtools)
else
  echo "Error: samtools is not available in PATH." >&2
  exit 1
fi

if ! command -v gawk >/dev/null 2>&1; then
  echo "Error: gawk is not available in PATH." >&2
  exit 1
fi

bam_base="$(basename "$bam_path")"
mkdir -p "$out_dir"
out_file="$out_dir/softclip_by_scaffold_window_${window_size}_${bam_base}.tsv"

tmp_scaf="$(mktemp)"
trap 'rm -f "$tmp_scaf"' EXIT

"${samtools_cmd[@]}" view -H "$bam_path" | gawk '
BEGIN { OFS="\t"; ord=0 }
$1=="@SQ" {
  sn=""
  ln=0
  for (i=2; i<=NF; i++) {
    if ($i ~ /^SN:/) sn=substr($i,4)
    else if ($i ~ /^LN:/) ln=substr($i,4)+0
  }
  if (sn!="" && ln>0) {
    ord++
    print sn, ln, ord
  }
}
' > "$tmp_scaf"

if [[ ! -s "$tmp_scaf" ]]; then
  echo "Error: no scaffold definitions found in BAM header." >&2
  exit 1
fi

# Exclude: 0x4 (unmapped), 0x100 (secondary), 0x200 (QC-fail), 0x400 (duplicate), 0x800 (supplementary)
"${samtools_cmd[@]}" view -F 3844 "$bam_path" | gawk -v OFS='\t' -v W="$window_size" '
function add_cov_segment(scaf, seg_start, seg_end,    ws, we, key, ov_start, ov_end) {
  if (seg_end < 1 || seg_start > scaf_len[scaf]) return
  if (seg_start < 1) seg_start = 1
  if (seg_end > scaf_len[scaf]) seg_end = scaf_len[scaf]

  ws = int((seg_start - 1) / W) * W + 1
  while (ws <= seg_end) {
    we = ws + W - 1
    if (we > scaf_len[scaf]) we = scaf_len[scaf]
    ov_start = (seg_start > ws ? seg_start : ws)
    ov_end = (seg_end < we ? seg_end : we)
    if (ov_end >= ov_start) {
      key = scaf SUBSEP ws
      cov_bases[key] += (ov_end - ov_start + 1)
    }
    ws += W
  }
}
BEGIN {
  print "scaffold","window_start","window_end","reads","reads_with_softclip","frac_reads_with_softclip","query_bases","softclip_bases","frac_query_bases_softclipped","avg_depth"
}
NR==FNR {
  scaf=$1
  len=$2+0
  ord=$3+0
  scaf_len[scaf]=len
  scaf_by_order[ord]=scaf
  max_order=(ord>max_order ? ord : max_order)
  next
}
{
  scaf=$3
  pos=$4+0
  cig=$6
  cig_work=cig
  if (!(scaf in scaf_len)) next
  if (scaf=="*" || cig=="*" || pos<=0) next

  win_start = int((pos-1)/W)*W + 1
  if (win_start > scaf_len[scaf]) next

  key=scaf SUBSEP win_start
  reads[key]++

  soft=0
  qbases=0
  ref_pos=pos
  while (match(cig_work, /([0-9]+)([MIDNSHP=X])/, m)) {
    len=m[1]+0
    op=m[2]
    if (op=="S") soft += len
    if (op ~ /[MIS=X]/) qbases += len
    if (op ~ /[M=X]/) {
      add_cov_segment(scaf, ref_pos, ref_pos + len - 1)
      ref_pos += len
    } else if (op ~ /[DN]/) {
      ref_pos += len
    }
    cig_work=substr(cig_work, RSTART+RLENGTH)
  }

  if (soft > 0) reads_soft[key]++
  soft_bases[key] += soft
  query_bases[key] += qbases
}
END {
  for (i=1; i<=max_order; i++) {
    scaf=scaf_by_order[i]
    len=scaf_len[scaf]
    for (ws=1; ws<=len; ws+=W) {
      we=ws+W-1
      if (we>len) we=len
      key=scaf SUBSEP ws

      r=reads[key]+0
      rs=reads_soft[key]+0
      qb=query_bases[key]+0
      sb=soft_bases[key]+0
      cb=cov_bases[key]+0
      wl=we-ws+1
      fr=(r>0 ? rs/r : 0)
      fb=(qb>0 ? sb/qb : 0)
      ad=(wl>0 ? cb/wl : 0)

      print scaf, ws, we, r, rs, fr, qb, sb, fb, ad
    }
  }
}
' "$tmp_scaf" - > "$out_file"

echo "Wrote: $out_file"
