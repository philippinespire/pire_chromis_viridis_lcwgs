#!/bin/bash
# ==============================================================================
# SCRIPT NAME: analyze_soft_clips.sh
# Malin Pinsky, 2026, with assistance from Gemini
#
# DESCRIPTION: Calculates the total count and fraction of soft-clipped bases 
#              per individual (BAM file) and per contig (aggregated across all 
#              individuals). Generates two distribution histograms (.png).
#
# DEPENDENCIES:
#   - samtools  (For streaming and filtering BAM alignments)
#   - awk       (For high-performance text and math parsing)
#   - R         (for plotting)
#
# USAGE:
#   1. Interactive / Local execution:
#      ./analyze_soft_clips.sbatch /path/to/bams/*.bam
#
#   2. Cluster Submission (SLURM):
#      sbatch analyze_soft_clips.sbatch /path/to/bams/*.bam
#
# ARGUMENTS:
#   <bam_file1> [bam_file2 ...] : Space-separated list of target BAM files.
#                                 Accepts wildcards (e.g., *.bam).
#
# OUTPUTS GENERATED:
#   In a directory specified by OUT_DIR (default: output/softclip_analysis):
#   - per_individual_stats.txt            : Tab-delimited metrics per sample.
#   - per_contig_stats.txt                : Tab-delimited metrics per contig.
#   - histogram_soft_clipped_per_individual.png : Sample-wide fraction distribution.
#   - histogram_soft_clipped_per_contig.png     : Contig-wide fraction distribution.
#
# NOTE: To avoid inflated counts, this script filters out unmapped, secondary, 
#       and supplementary alignments (-F 2308).
# ==============================================================================

#SBATCH --job-name=soft_clip
#SBATCH --output=logs/soft_clip_%j.out
#SBATCH --error=logs/soft_clip_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2          # Samtools and awk will share resources safely
#SBATCH --mem=16G                  # Adjust based on BAM index sizes/depth
#SBATCH --time=04:00:00            # Adjust walltime based on cohort size

set -euo pipefail

# --- Usage & Dependency Check ---
# When running via sbatch, pass your BAMs as arguments, e.g., sbatch analyze_soft_clips.sh *.bam
if [ "$#" -lt 1 ]; then
    echo "Error: No BAM files provided." >&2
    echo "Usage: sbatch $0 <bam_file1> [bam_file2 ...]" >&2
    exit 1
fi

# --- Initialization ---
#OUT_DIR="output/softclip_analysis" # first run on in-house genome
OUT_DIR="output/softclip_analysis-iridian" # second run on iridian genome
mkdir -p "$OUT_DIR"

OUTPUT_IND="$OUT_DIR/per_individual_stats.txt"
OUTPUT_CTG="$OUT_DIR/per_contig_stats.txt"
TMP_DIR="$OUT_DIR/tmp_soft_clip_counts_${SLURM_JOB_ID:-local}"

echo -e "Sample\tSoftClippedBases\tTotalBases\tFraction" > "$OUTPUT_IND"
mkdir -p "$TMP_DIR"

# --- Step 1: Process Each BAM File ---
module load container_env samtools  # Ensure samtools is available for BAM processing

for bam in "$@"; do
    if [ ! -f "$bam" ]; then
        echo "Warning: File '$bam' not found. Skipping..." >&2
        continue
    fi
    
    sample=$(basename "$bam" .bam)
    echo "Processing sample: $sample..."

    # -F 2308 excludes unmapped (4), secondary (256), and supplementary (2048) alignments
    crun samtools view -F 2308 "$bam" | awk -v sample="$sample" '
    {
        cigar = $6; contig = $3; len = length($10);
        soft = 0;
        
        # Fast regex matching to extract and sum all soft-clipped (S) lengths per read
        while (match(cigar, /[0-9]+S/)) {
            soft += substr(cigar, RSTART, RLENGTH-1);
            cigar = substr(cigar, RSTART+RLENGTH);
        }
        
        ind_soft += soft; 
        ind_len += len;
        ctg_soft[contig] += soft; 
        ctg_len[contig] += len;
    }
    END {
        # Output individual totals
        print sample, ind_soft, ind_len, (ind_len > 0 ? ind_soft/ind_len : 0) >> "'"$OUTPUT_IND"'"
        
        # Stream contig totals to a temporary file
        for (c in ctg_soft) {
            print c, ctg_soft[c], ctg_len[c] >> "'"$TMP_DIR"'/contigs.tmp"
        }
    }' OFS="\t"
done

# --- Step 2: Aggregate Contig Data Across All Individuals ---
echo "Aggregating contig metrics across individuals..."
echo -e "Contig\tSoftClippedBases\tTotalBases\tFraction" > "$OUTPUT_CTG"

if [ -f "$TMP_DIR/contigs.tmp" ]; then
    awk '
    {
        ctg_soft[$1] += $2;
        ctg_len[$1] += $3;
    }
    END {
        for (c in ctg_soft) {
            print c, ctg_soft[c], ctg_len[c], (ctg_len[c] > 0 ? ctg_soft[c]/ctg_len[c] : 0);
        }
    }' OFS="\t" "$TMP_DIR/contigs.tmp" >> "$OUTPUT_CTG"
fi

# Clean up temporary storage
rm -rf "$TMP_DIR"

# --- Step 3: Plot Histograms using R ---
echo "Generating plots with R..."
module load container_env R
crun Rscript - "$OUT_DIR" << 'EOF'

# Force R to use Cairo instead of X11 for headless server plotting
options(bitmapType = 'cairo')

# Catch the directory path passed from Bash
args <- commandArgs(trailingOnly = TRUE)
dir_path <- args[1]

# create a plotting function to avoid repetition
plot_histogram <- function(filename, col_name, color_hex, output_png, title, xlab, ylab, bins) {
  file_path <- file.path(dir_path, filename)
  if (!file.exists(file_path)) return(NULL)
  
  # Read dataset
  data <- read.table(file_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  if (nrow(data) == 0) return(NULL)
  
  # Initialize PNG graphic device
  png(file.path(dir_path, output_png), width=7, height=4.5, units="in", res=300, type="cairo")  

  # Set clean margins
  par(mar=c(4.5, 4.5, 2.5, 1.5))
  
  # Pre-calculate histogram to layer background grids cleanly underneath the bars
  h <- hist(data[[col_name]], breaks=bins, plot=FALSE)
  
  # Draw an empty plot framework first
  plot(h, col=NA, border=NA, main=title, xlab=xlab, ylab=ylab, las=1, 
       cex.main=1.1, cex.lab=0.9, cex.axis=0.8)
  
  # Add soft horizontal grid lines
  grid(nx=NA, ny=NULL, lty="dashed", col="gray85")
  
  # Overlay the styled histogram bars
  plot(h, col=color_hex, border="black", add=TRUE)
  
  dev.off()
}

# Plot 1: Per Individual (Blue theme)
plot_histogram(
  filename = "per_individual_stats.txt", 
  col_name = "Fraction", 
  color_hex = "#4C72B0",
  output_png = "histogram_soft_clipped_per_individual.png",
  title = "Soft-Clipped Base Fraction Per Individual",
  xlab = "Fraction of Soft-Clipped Bases",
  ylab = "Count of Individuals",
  bins = 15
)

# Plot 2: Per Contig (Orange theme)
plot_histogram(
  filename = "per_contig_stats.txt", 
  col_name = "Fraction", 
  color_hex = "#DD8452",
  output_png = "histogram_soft_clipped_per_contig.png",
  title = "Soft-Clipped Base Fraction Per Contig (All Samples)",
  xlab = "Fraction of Soft-Clipped Bases",
  ylab = "Count of Contigs",
  bins = 20
)

print("Plotting complete!")
EOF

echo "Analysis finished successfully."