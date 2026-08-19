#!/bin/bash
#SBATCH --job-name=compare_loci
#SBATCH --output=logs/compare_loci_%j.out
#SBATCH --error=logs/compare_loci_%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4

# ==============================================================================
# Script Name   : run_loci_comparison.sbatch
# Description   : Extracts site coordinates from historic and modern ANGSD SAF
#                 index files (.saf.idx) in parallel, then streams a two-pointer 
#                 comparison to quantify shared vs. era-unique genomic loci.
#                 Be sure to set the input *.saf.idx files to compare on lines 28-32
# Author        : Malin Pinsky with Gemini (Google AI)
# Usage         : sbatch run_loci_comparison.sbatch
# Dependencies  : realSFS (ANGSD), Python 3
# ==============================================================================

set -e

mkdir -p temp

#echo "[$(date)] Extracting historic and modern loci in parallel..."
module load container_env angsd python3

# Extract historic loci in background
#crun.angsd realSFS print nf-pipelines/nf-angsd-diversity-generode/results/angsd_pop/CviAPal_historic.saf.idx | cut -f1,2 > temp/hist.pos &
#PID_HIST=$!

# Extract modern loci in background
#crun.angsd realSFS print nf-pipelines/nf-angsd-diversity-generode/results/angsd_pop/CviCPal_modern.saf.idx | cut -f1,2 > temp/mod.pos &
#PID_MOD=$!

# Wait for both extraction tasks to finish
#wait $PID_HIST $PID_MOD

echo "[$(date)] Extractions complete. Streaming comparison..."

# Stream and compare files using embedded Python
crun.python3 python3 - temp/hist.pos temp/mod.pos << 'EOF'
import sys

hist_path, mod_path = sys.argv[1], sys.argv[2]

chr_order = {}
def get_chr_idx(chrom):
    if chrom not in chr_order:
        chr_order[chrom] = len(chr_order)
    return chr_order[chrom]

hist_unique = mod_unique = shared = 0

with open(hist_path, 'r') as fh, open(mod_path, 'r') as fm:
    def parse(line):
        if not line: return None
        p = line.strip().split()
        return (p[0], int(p[1])) if len(p) >= 2 else None

    h = parse(fh.readline())
    m = parse(fm.readline())

    while h and m:
        h_cidx, m_cidx = get_chr_idx(h[0]), get_chr_idx(m[0])
        
        if h_cidx < m_cidx or (h_cidx == m_cidx and h[1] < m[1]):
            hist_unique += 1
            h = parse(fh.readline())
        elif m_cidx < h_cidx or (h_cidx == m_cidx and m[1] < h[1]):
            mod_unique += 1
            m = parse(fm.readline())
        else:
            shared += 1
            h = parse(fh.readline())
            m = parse(fm.readline())

    while h:
        hist_unique += 1
        h = parse(fh.readline())
    while m:
        mod_unique += 1
        m = parse(fm.readline())

print("\n================ LOCI COMPARISON RESULTS ================")
print(f"Historic-only loci: {hist_unique:,}")
print(f"Modern-only loci:   {mod_unique:,}")
print(f"Shared loci:        {shared:,}")
print("=========================================================\n")
EOF

echo "[$(date)] Job finished successfully."