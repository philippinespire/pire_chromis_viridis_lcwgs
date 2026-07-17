#!/bin/bash
#SBATCH --job-name=depth_vs_reads
#SBATCH --output=logs/depth_vs_reads_%j.out
#SBATCH --error=logs/depth_vs_reads_%j.err
#SBATCH --ntasks=1
#SBATCH --time=10:00:00

# Usage: sbatch scripts/plot_depth_vs_reads.sh <path_to_fastq_dir> <path_to_depth_dir> [output_data.txt] [output_plot.pdf]
# The fastq files must be of the form {indiv_id}_*.fastq.gz
# Depth files must be of the form {indiv_id}.*.bam.dpstats.txt with one number per file
#
# Malin Pinsky, 2026 with the assistance of Google Gemini


FASTQ_DIR=$1
DEPTH_DIR=$2
OUTPUT_DATA=${3:-"output/depth_vs_reads.txt"}
OUTPUT_PLOT=${4:-"output/depth_vs_reads.pdf"}

if [ -z "$FASTQ_DIR" ] || [ -z "$DEPTH_DIR" ]; then
	echo "Usage: $0 <path_to_fastq_dir> <path_to_depth_dir> [output_data.txt] [output_plot.pdf]"
    echo ""
    echo "Optional arguments:"
    echo "  [output_data.txt]  Path to save the summary table (default: plot_data.txt)"
    echo "  [output_plot.pdf]  Path to save the generated plot (default: reads_vs_depth.pdf)"
    exit 1
fi

echo -e "Individual\tTotalReads\tMappedDepth" > "$OUTPUT_DATA"

echo "Extracting data and counting reads. This may take some time depending on fastq sizes..."

# Iterate over each depth file
for depth_file in "$DEPTH_DIR"/*.bam.dpstats.txt; do
    # 1. Extract Individual ID (e.g., CviAPal001 from CviAPal001.Q25.L107.bam.dpstats.txt)
    filename=$(basename "$depth_file")
    indiv_id="${filename%%.*}" 
    
    # 2. Extract Depth (assumes the file contains a single number)
    depth=$(cat "$depth_file" | tr -d '[:space:]')
    
    # 3. Count Total Reads
    # Matches all fastq.gz files for this individual (all lanes, R1 & R2)
    fq_pattern="${FASTQ_DIR}/${indiv_id}_*.fastq.gz"
    
    # Check if files exist to avoid errors
    if ls $fq_pattern 1> /dev/null 2>&1; then
        # zcat all matching files, count lines, and divide by 4 for read count
        total_lines=$(zcat $fq_pattern | wc -l)
        total_reads=$((total_lines / 4))
    else
        total_reads=0
        echo "Warning: No fastq files found for $indiv_id"
    fi
    
    echo -e "${indiv_id}\t${total_reads}\t${depth}" >> "$OUTPUT_DATA"
done

echo "Data extraction complete. Generating plot..."

# 4. Generate Plot using Rscript
module load container_env R
crun Rscript - <<EOF
# Load ggplot2 quietly
suppressPackageStartupMessages(library(ggplot2))

# Read the generated data
df <- read.table("$OUTPUT_DATA", header=TRUE, sep="\t")

# Ensure numeric types
df\$TotalReads <- as.numeric(df\$TotalReads)
df\$MappedDepth <- as.numeric(df\$MappedDepth)

# Create the scatter plot
p <- ggplot(df, aes(x=TotalReads, y=MappedDepth)) +
    geom_point(alpha=0.7, color="#0072B2", size=3) +
    geom_smooth(method="lm", color="black", linetype="dashed", se=FALSE) +
    scale_x_continuous(labels = scales::comma) +
    theme_minimal(base_size=14) +
    labs(title="Mapped Depth vs. Total Reads per Individual",
         subtitle=paste("N =", nrow(df), "individuals"),
         x="Total Reads",
         y="Mapped Depth") +
    theme(
        plot.title = element_text(face="bold"),
        panel.grid.minor = element_blank()
    )

# Save the plot
ggsave("$OUTPUT_PLOT", plot=p, width=8, height=6, dpi=300)
EOF

echo "Done! Plot saved as $OUTPUT_PLOT and data saved as $OUTPUT_DATA"