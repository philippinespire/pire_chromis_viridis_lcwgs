#!/usr/bin/env bash
set -euo pipefail

# Malin Pinsky 2026
# This script was written with AI assistance using GitHub Copilot
# with the GPT-5.3-Codex model.

# Run remote BLASTn against NCBI nt and keep top 2 hits per query,
# then write a compact summary table.
#
# Usage:
#   scripts/blastn_remote_top2.sh [query_fasta] [out_prefix]
#
# Defaults:
#   query_fasta = output/MitoZ_output.fasta
#   out_prefix  = output/MitoZ_output_vs_nt

QUERY_FASTA="${1:-output/MitoZ_output.fasta}"
OUT_PREFIX="${2:-output/MitoZ_output_vs_nt}"
SOFTWARE_DIR="${BLAST_SOFTWARE_DIR:-/archive/carpenterlab/pire/mpinsky/software}"
BLAST_TASK="${BLAST_TASK:-blastn}" # megablast would be cheaper
MAX_RETRIES="${BLAST_MAX_RETRIES:-3}"
RETRY_SLEEP_SECONDS="${BLAST_RETRY_SLEEP_SECONDS:-5}"
PARALLEL_JOBS="${BLAST_PARALLEL_JOBS:-3}"
SUBMISSION_DELAY_SECONDS="${BLAST_SUBMISSION_DELAY_SECONDS:-10}"

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

if ! [[ "$PARALLEL_JOBS" =~ ^[0-9]+$ ]] || [[ "$PARALLEL_JOBS" -lt 1 ]]; then
  echo "ERROR: BLAST_PARALLEL_JOBS must be a positive integer (got: $PARALLEL_JOBS)" >&2
  exit 1
fi

if ! [[ "$SUBMISSION_DELAY_SECONDS" =~ ^[0-9]+$ ]]; then
  echo "ERROR: BLAST_SUBMISSION_DELAY_SECONDS must be a non-negative integer (got: $SUBMISSION_DELAY_SECONDS)" >&2
  exit 1
fi

echo "Running remote BLASTn on $QUERY_FASTA ..."
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

# Split multi-fasta into one-query files to reduce NCBI remote CPU-limit failures.
awk -v outdir="$TMP_DIR" '
  /^>/ {
    n++
    file=sprintf("%s/query_%06d.fasta", outdir, n)
  }
  { print >> file }
' "$QUERY_FASTA"

if ! ls "$TMP_DIR"/query_*.fasta >/dev/null 2>&1; then
  echo "ERROR: No FASTA entries found in $QUERY_FASTA" >&2
  exit 1
fi

query_runner="$TMP_DIR/run_one_query.sh"
cat > "$query_runner" <<'RUNNER'
#!/usr/bin/env bash
set -euo pipefail

query_fasta="$1"
query_tsv="${query_fasta%.fasta}.tsv"
fail_marker="${query_fasta%.fasta}.failed"
sequence_id="$(awk 'NR==1{sub(/^>/, "", $0); print $1; exit}' "$query_fasta")"

for ((attempt=1; attempt<=MAX_RETRIES; attempt++)); do
  ts="$(date '+%Y-%m-%d %H:%M:%S')"
  echo "[$ts] Submitting sequence ${sequence_id} (attempt $attempt/$MAX_RETRIES)"
  if "$BLASTN_BIN" \
    -query "$query_fasta" \
    -db nt \
    -remote \
    -task "$BLAST_TASK" \
    -evalue 1e-20 \
    -max_target_seqs 50 \
    -max_hsps 1 \
    -outfmt "6 qseqid sacc pident length qcovs evalue bitscore staxids sscinames stitle" \
    -out "$query_tsv"; then
    exit 0
  fi

  echo "WARNING: BLAST failed for $(basename "$query_fasta") attempt $attempt/$MAX_RETRIES" >&2
  sleep "$RETRY_SLEEP_SECONDS"
done

echo "$query_fasta" > "$fail_marker"
echo "WARNING: Giving up on $(basename "$query_fasta") after $MAX_RETRIES attempts" >&2
exit 0
RUNNER
chmod +x "$query_runner"

export BLASTN_BIN BLAST_TASK MAX_RETRIES RETRY_SLEEP_SECONDS
echo "Launching BLAST requests with PARALLEL_JOBS=$PARALLEL_JOBS and SUBMISSION_DELAY_SECONDS=$SUBMISSION_DELAY_SECONDS ..."
for query_fasta in "$TMP_DIR"/query_*.fasta; do
  while [[ "$(jobs -rp | wc -l)" -ge "$PARALLEL_JOBS" ]]; do
    sleep 1
  done

  "$query_runner" "$query_fasta" &
  sleep "$SUBMISSION_DELAY_SECONDS"
done
wait

: > "$ALL_TSV"
for query_fasta in "$TMP_DIR"/query_*.fasta; do
  query_tsv="${query_fasta%.fasta}.tsv"
  if [[ -s "$query_tsv" ]]; then
    cat "$query_tsv" >> "$ALL_TSV"
  fi
done

failed_queries=$(find "$TMP_DIR" -maxdepth 1 -type f -name 'query_*.failed' | wc -l)

if [[ ! -s "$ALL_TSV" ]]; then
  echo "ERROR: No BLAST hits were collected. Check remote BLAST errors in stdout/stderr." >&2
  exit 1
fi

if [[ "$failed_queries" -gt 0 ]]; then
  echo "WARNING: $failed_queries query(ies) failed after retries and were skipped." >&2
fi

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
