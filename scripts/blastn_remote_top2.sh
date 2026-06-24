#!/usr/bin/env bash
set -euo pipefail

# Run remote BLASTn against NCBI nt and keep top 2 hits per query,
# then write a compact summary table.
#
# Usage:
#   scripts/blastn_remote_top2.sh [query_fasta] [out_prefix]
#
# Defaults:
#   query_fasta = mitoz/MitoZ_output.fasta
#   out_prefix  = mitoz/MitoZ_output_vs_nt

QUERY_FASTA="${1:-mitoz/MitoZ_output.fasta}"
OUT_PREFIX="${2:-mitoz/MitoZ_output_vs_nt}"
SOFTWARE_DIR="${BLAST_SOFTWARE_DIR:-/archive/carpenterlab/pire/mpinsky/software}"

ALL_TSV="${OUT_PREFIX}_all.tsv"
TOP2_TSV="${OUT_PREFIX}_top2.tsv"
COMPACT_TSV="${OUT_PREFIX}_top2_compact.tsv"

if [[ ! -f "$QUERY_FASTA" ]]; then
  echo "ERROR: Query FASTA not found: $QUERY_FASTA" >&2
  exit 1
fi

if [[ ! -d "$SOFTWARE_DIR" ]]; then
  echo "ERROR: BLAST software directory not found: $SOFTWARE_DIR" >&2
  exit 1
fi

BLASTN_BIN="$(find "$SOFTWARE_DIR" -type f -name blastn -perm -u+x 2>/dev/null | head -n 1)"
if [[ -z "$BLASTN_BIN" ]]; then
  echo "ERROR: No executable blastn found under $SOFTWARE_DIR" >&2
  exit 1
fi

echo "Running remote BLASTn on $QUERY_FASTA ..."
"$BLASTN_BIN" \
  -query "$QUERY_FASTA" \
  -db nt \
  -remote \
  -task blastn \
  -evalue 1e-20 \
  -max_target_seqs 50 \
  -max_hsps 1 \
  -outfmt "6 qseqid sacc pident length qcovs evalue bitscore staxids sscinames stitle" \
  -out "$ALL_TSV"

echo "Selecting top 2 hits per query ..."
awk 'BEGIN{FS=OFS="\t"} {if(++n[$1] <= 2) print}' "$ALL_TSV" > "$TOP2_TSV"

echo "Writing compact summary table ..."
awk 'BEGIN{
  FS=OFS="\t";
  print "qseqid","rank","sacc","species","pident","qcovs","evalue","bitscore","title"
}
{
  r[$1]++;
  title=$10;
  if (length(title) > 90) title=substr(title,1,87)"...";
  print $1,r[$1],$2,$9,$3,$5,$6,$7,title
}' "$TOP2_TSV" > "$COMPACT_TSV"

echo "Done."
echo "All hits:        $ALL_TSV"
echo "Top 2 per query: $TOP2_TSV"
echo "Compact table:   $COMPACT_TSV"
