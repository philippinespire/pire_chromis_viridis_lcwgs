#!/usr/bin/env bash
#SBATCH -n 1
#SBATCH --time 4-00:00:00
#SBATCH --job-name=MIA
#SBATCH --output=/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/logs/run_MIA_%j.out

### MIA (Mapping Iterative Assembler): for mitogenome reconstruction
# Call it with one argument: the path to the .fastq.gz file to run, eg
# sbatch run_MIA.sbatch /path/to/fastq/my.fastq.gz
# If the fastq file is paired-end and ends with _R1.fastq.gz, it will 
# look for the second read in the same directory with the same name but ending in _R2.fastq.gz

# set the path to mia and ma
export PATH=/archive/carpenterlab/pire/softwares/mapping-iterative-assembler/bin:$PATH

### Read arguments
SAMPLE_PATH=$1

### Set up user parameters
MIA_OUTDIR=/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/output/mia/
ANCIENT_DNA_MATRIX=/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/data/ancient.submat.solexa.pe.txt
MIA_REFERENCE_MITOGENOME=/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/data/chromis_viridis_mitogenome_MT199208.1.fasta
MIA_FASTA_PY=/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/scripts/create_fasta_from_mia.py

### Set up other parameters from arguments so far
# find file name by stripping off the directory
FILE_NAME="${SAMPLE_PATH##*/}"

# Determine sample name and check for R2 dynamically
if [[ "$FILE_NAME" == *"_R1"* ]]; then
    # Paired-end path: Strip "_R1" and everything after it for the clean sample name
    SAMPLE_NAME="${FILE_NAME%%_R1*}"
    R2_PATH="${SAMPLE_PATH/_R1/_R2}"
    
    # Check if the R2 file actually exists on disk
    if [[ -f "$R2_PATH" ]]; then
        echo "LOG: Paired-end reads detected. Streaming R1 and R2."
        READS_TO_STREAM="$SAMPLE_PATH $R2_PATH"
    else
        echo "WARNING: '_R1' found, but matching R2 file not found at: $R2_PATH"
        echo "LOG: Falling back to Single-End mode (R1 only)."
        READS_TO_STREAM="$SAMPLE_PATH"
    fi
else
    # Single-end path: Strip extensions (after the first dot) for the clean sample name
    SAMPLE_NAME="${FILE_NAME%%.*}"
    echo "LOG: Single-end read detected. Streaming R1 only."
    READS_TO_STREAM="$SAMPLE_PATH"
fi

# make output directory if it doesn't exist
mkdir -p ${MIA_OUTDIR}

# Running MIA -- Shotgun data, not too diverged from reference (if capture, remove duplicates; if diverged, use wormhole MIA)
	# slow approach with gzip
	# try to replace with zcat and stdin
	#gzip -d $PROCESSED_READS/${SAMPLE_NAME}

	# run mia from the fastq file
	zcat $READS_TO_STREAM | mia -r $MIA_REFERENCE_MITOGENOME -f /dev/stdin -c -C -U -s ${ANCIENT_DNA_MATRIX} -i -F -k 14 -m ${MIA_OUTDIR}/${SAMPLE_NAME}.maln

	# slow approach with gzip
	# try to replace with zcat and stdin
	# gzip $PROCESSED_READS/${SAMPLE_NAME}_trimmed_merged.fastq

	# output the column summary format
	ma -M ${MIA_OUTDIR}/${SAMPLE_NAME}.maln.* -f 3 > ${MIA_OUTDIR}/${SAMPLE_NAME}.maln.F.mia_stats.txt

	# output the line summary format
	ma -M ${MIA_OUTDIR}/${SAMPLE_NAME}.maln.* -f 2 > ${MIA_OUTDIR}/${SAMPLE_NAME}.maln.F.mia_coverage_per_site.txt

	# output the fasta format
	ma -M ${MIA_OUTDIR}/${SAMPLE_NAME}.maln.* -f 5 > ${MIA_OUTDIR}/${SAMPLE_NAME}.maln.F.mia_consensus.fasta

	# output the 13-column format for the python script
	ma -M ${MIA_OUTDIR}/${SAMPLE_NAME}.maln.* -f 41 > ${MIA_OUTDIR}/${SAMPLE_NAME}.maln.F.inputfornext.txt


module load container_env python3

# Generating consensus mitogenomes
	# -c = minimum depth of coverage required for a base to be called
	# -p = proportion of read agreement required for a base to be called
	python3 ${MIA_FASTA_PY} -c 3 -p 0.67 -I ${MIA_OUTDIR}/${SAMPLE_NAME}_3x_0.67 -m ${MIA_OUTDIR}/${SAMPLE_NAME}.maln.F.inputfornext.txt -o ${MIA_OUTDIR}/${SAMPLE_NAME}.maln.F.mia_consensus.3x_0.67_filtered.fasta

	python3 ${MIA_FASTA_PY} -c 10 -p 0.9 -I ${MIA_OUTDIR}/${SAMPLE_NAME}_10x_0.9 -m ${MIA_OUTDIR}/${SAMPLE_NAME}.maln.F.inputfornext.txt -o ${MIA_OUTDIR}/${SAMPLE_NAME}.maln.F.mia_consensus.10x_0.9_filtered.fasta


#Change name of alignment to something more informative
#This is more easy to do infile
#sed -i "1s/.*/>$SAMPLE_NAME/" ${MIA_OUTDIR}/${SAMPLE_NAME}.maln.F.mia_consensus.10x_0.9_filtered.fasta