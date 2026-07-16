#!/usr/bin/env bash
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=40G
#SBATCH --time=04:00:00
#SBATCH --job-name=blast_mitogenomes
#SBATCH --output=/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/logs/blast_%j.out

# Malin Pinsky 2026
# with assistance from Google Gemini
# Script to run BLASTn on a set of filtered consensus fasta files against 
# a local BLAST database, keeping the top 2 hits per query and writing a compact summary table.
# assumes that the files are of the form (see line 36)
#       <Individual_ID>_Ex1_L4.maln.F.mia_consensus.10x_0.9_filtered.fasta
# check the paths below if reusing this script
# run as
#       sbatch scripts/blastn_local_top2.sbatch

# Adds our Wahab compiled BLAST binaries to the PATH
export PATH=/archive/carpenterlab/pire/softwares/ncbi-blast-2.17.0+/bin:$PATH

# Define paths
FASTA_DIR="/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/output/mia"
BLAST_DB="/archive/carpenterlab/pire/pire_blast/nt"
OUT_DIR="/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/output/mia"
MASTER_OUT="${OUT_DIR}/mia_CviAPal_blast_results.txt"

# Ensure output directory exists
mkdir -p "$OUT_DIR"

# Write a clean, tab-separated header to your master output file
echo -e "Individual_ID\tQuery_ID\tSubject_Accession\tPercent_Identity\tAlignment_Length\tE_Value\tBit_Score\tSpecies_Title" > "$MASTER_OUT"

echo "Starting BLAST database queries..."
echo "--------------------------------------------------"

# Loop through all 10x filtered consensus fasta files
for fasta_path in "${FASTA_DIR}"/*_Ex1_L4.maln.F.mia_consensus.10x_0.9_filtered.fasta; do
    if [[ -f "$fasta_path" ]]; then
        # Extract the base filename
        filename=$(basename "$fasta_path")
        
        # Extract the clean Individual ID (e.g., CviAPal031) by stripping everything after and including "_Ex1"
        ind_id="${filename%%_Ex1*}"
        
        echo "Processing: ${ind_id}"

        # Run blastn:
        #  -num_threads: utilizes the 8 cores requested in the Slurm header
        #  -max_target_seqs 2: limits results to the top 2 matching reference subjects
        #  -max_hsps 1: keeps only the single best alignment segment per subject (prevents messy duplication)
        #  -outfmt 6 (with custom fields): outputs query name, subject accession, identity %, alignment length, e-value, bit score, and species description
        #  awk: dynamically injects the Individual ID as a new first column
        blastn \
            -query "$fasta_path" \
            -db "$BLAST_DB" \
            -num_threads "$SLURM_CPUS_PER_TASK" \
            -max_target_seqs 2 \
            -max_hsps 1 \
            -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" | \
            awk -v id="$ind_id" 'BEGIN {FS=OFS="\t"} {print id, $0}' >> "$MASTER_OUT"
    fi
done

echo "--------------------------------------------------"
echo "BLAST complete! Master output written to: $MASTER_OUT"