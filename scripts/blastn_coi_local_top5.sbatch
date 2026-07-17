#!/usr/bin/env bash
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=40G
#SBATCH --time=02:00:00
#SBATCH --job-name=blast_coi
#SBATCH --output=logs/blastn_coi_%j.out

# Malin Pinsky 2026
# with assistance from Google Gemini
# Script to run BLASTn on consolidated master COI multi-FASTA files
# (eg, output by scripts/cut_COI_from_mia.sbatch) against 
# a local BLAST database, keeping the top 2 hits per individual query.

# Set in-file names
#infile1="all_samples_10x_0.9_coi.fasta" # for COI barcode
#infile2="all_samples_3x_0.67_coi.fasta"
infile1="all_samples_10x_0.9_coi_full.fasta" # for full COI
infile2="all_samples_3x_0.67_coi_full.fasta"

# Define paths
FASTA_DIR="/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/output/mia"
BLAST_DB="/archive/carpenterlab/pire/pire_blast/nt/nt"
OUT_DIR="/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/output/mia"
MASTER_OUT="${OUT_DIR}/all_samples_coi_blast_results.txt"


# Adds the Wahab compiled BLAST binaries to the PATH
export PATH=/archive/carpenterlab/pire/softwares/ncbi-blast-2.17.0+/bin:$PATH

# Ensure output directory exists
mkdir -p "$OUT_DIR"

# Write a clean, tab-separated header to the master output file
echo -e "Individual_ID\tQuery_ID\tSubject_Accession\tPercent_Identity\tAlignment_Length\tE_Value\tBit_Score\tSpecies_Title" > "$MASTER_OUT"

echo "Starting BLAST queries for master multi-FASTA files..."
echo "--------------------------------------------------"

# Loop directly through the two master files
for file_name in $infile1 $infile2; do
    fasta_path="${FASTA_DIR}/${file_name}"
    
    if [[ -f "$fasta_path" ]]; then
        echo "[$(date +'%Y-%m-%d %H:%M:%S')]: BLASTing ${file_name}..."

        
        # Run blastn:
        #  -num_threads: utilizes the 8 cores requested in the Slurm header
        #  -max_target_seqs 5: limits results to the top 5 matching reference subjects PER individual query
        #  -max_hsps 1: keeps only the single best alignment segment per subject
        #  awk: takes column 1 (e.g., CviAPal031_10x_0.9_COI), splits it by the underscore,
        #       and prints the first chunk (CviAPal031) as a new leading column.
        blastn \
            -query "$fasta_path" \
            -db "$BLAST_DB" \
            -num_threads "$SLURM_CPUS_PER_TASK" \
            -max_target_seqs 5 \
            -max_hsps 1 \
            -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" | \
            awk 'BEGIN {FS=OFS="\t"} {split($1, arr, "_"); print arr[1], $0}' >> "$MASTER_OUT"
    else
        echo "WARNING: Expected master file not found at $fasta_path"
    fi
done

echo "--------------------------------------------------"
echo "BLAST complete! Combined master table written to: $MASTER_OUT"