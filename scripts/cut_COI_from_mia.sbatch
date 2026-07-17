#!/bin/bash
#SBATCH --job-name=cut_COI
#SBATCH --output=logs/cut_coi_from_mia_%j.out
#SBATCH --error=logs/cut_coi_from_mia_%j.err
#SBATCH --ntasks=1
#SBATCH --time=1:00:00

# Malin Pinsky, 2026 with the assistance of Google Gemini

### Parameters
# set the target COI sequence for BLAST extraction
# target="data/chromis_viridis_COI_barcode_ref.fasta" # for COI barcode
target="data/chromis_viridis_COI_full.fasta" # for full COI gene

# set the output files
#out10x="output/mia/all_samples_10x_0.9_coi.fasta" # for barcode
#out3x="output/mia/all_samples_3x_0.67_coi.fasta" # for barcode
out10x="output/mia/all_samples_10x_0.9_coi_full.fasta" # for full gene
out3x="output/mia/all_samples_3x_0.67_coi_full.fasta" # for full gene

# Ensure local BLAST binaries are in the path
export PATH=/archive/carpenterlab/pire/softwares/ncbi-blast-2.17.0+/bin:$PATH

# to hold temporary COIs from individuals
mkdir -p temp/extracted_coi

# Loop through the MIA consensus files and isolate the COI gene
for file in output/mia/*.maln.F.mia_consensus.*_filtered.fasta; do
    if [[ -f "$file" ]]; then
        filename=$(basename "$file")
        
        # Extract clean individual ID (e.g., CviAPal031)
        sample="${filename%%_Ex1*}"
        
        # Detect the MIA fasta filtering strategy
        if [[ "$filename" == *"10x_0.9"* ]]; then
            strategy="10x_0.9"
        elif [[ "$filename" == *"3x_0.67"* ]]; then
            strategy="3x_0.67"
        else
            strategy="unknown"
        fi
        
        echo "Processing: $sample ($strategy)"
        
        # Define output path with strategy included to prevent collisions
        out_file="temp/extracted_coi/${sample}_${strategy}_coi.fasta"
        
        # Write unique FASTA header (e.g., >CviAPal031_10x_0.9_COI)
        echo ">${sample}_${strategy}_COI" > "$out_file"
        
        # Run BLAST extraction
        blastn -query "$target" \
               -subject "$file" \
               -max_target_seqs 1 \
               -max_hsps 1 \
               -outfmt "6 sseq" | tr -d '-' >> "$out_file"
    fi
done

#Generate separate master files for each strategy
cat temp/extracted_coi/*_10x_0.9_coi.fasta > $out10x
cat temp/extracted_coi/*_3x_0.67_coi.fasta > $out3x