// Marianne Dehasque and Malin Pinsky, 2026
// Pipeline optimization, dynamic channel routing (DSL2), and logic structuring 
// were developed with the assistance of GPT5.3 and Google Gemini (July 2026).

nextflow.enable.dsl=2
include { BWA_MERGED as BWA_MERGED_PASS1; BWA_MERGED as BWA_MERGED_PASS2 } from './mapping_modules.nf'
include { BWA_UNMERGED as BWA_UNMERGED_PASS1; BWA_UNMERGED as BWA_UNMERGED_PASS2 } from './mapping_modules.nf'

// --- Default Parameters ---
params.samplesheet = "${projectDir}/inputfiles/samplesheet.csv" // Preferred way to provide sample metadata, including sample IDs and eras (modern or historical).
params.samples_file = "${projectDir}/inputfiles/fastq_filenames.txt" // Legacy way to provide sample metadata. One-column text file with sample IDs. All modern by default.
params.indir        = "${projectDir}/data/symlinks" // Directory where raw FASTQ files are expected.
params.outdir       = "${projectDir}/results" // Directory where all output files will be written.
params.reference    = "${projectDir}/data/reference/reference.ssl.Cvi20k_rename.fasta" // Path to reference genome FASTA file.
params.bed_file     = "${projectDir}/data/reference/reference.ssl.Cvi20k_rename.repma.bed" // input bed file. Only used if params.run_repeatmasking is set to false.
params.reference_prefix         = params.reference.tokenize('/').last().replaceAll(/\.(fa|fasta|fna)$/, '') // Extracts base name of reference.
params.historical_era           = "historical"
params.historical_mapper        = "mem" // Default starting point for historical read mapping ("aln" or "mem")
params.run_repeatmasking        = true // Whether to run RepeatModeler and RepeatMasker on the reference genome (true/false). If false, it needs an input bed file of repeats and CpG sites for downstream ANGSD analyses.
params.run_historical_fastqc    = true
params.run_historical_mapdamage = true
params.run_historical_amber     = true
params.run_modern_amber         = true
params.split_script = "${projectDir}/scripts/split_reads.sh"
params.rmdup_script = "${projectDir}/scripts/samremovedup.py"
params.amber_script = "/archive/carpenterlab/pire/softwares/AMBERv2/AMBER" 
params.bwa_threads  = 4
params.bam_q        = 25 // Mapping quality threshold. 
params.trimlength   = 85 // Default fallback trim length if no historical reads map successfully

def resolve_reads = { String sample_id ->
    def r1 = file("${params.indir}/${sample_id}_R1.fastq.gz")
    def r2 = file("${params.indir}/${sample_id}_R2.fastq.gz")

    if( !r1.exists() ) error "R1 file not found for ${sample_id}: ${r1}"
    if( !r2.exists() ) error "R2 file not found for ${sample_id}: ${r2}"

    return [r1, r2]
}

def samplesheet_file = file(params.samplesheet)
def legacy_samples_file = file(params.samples_file)

// --- Input Channels ---
// Sample metadata is taken from a CSV when available.
// Fallback: legacy one-column sample list, which defaults all samples to modern.
sample_metadata_ch = samplesheet_file.exists() ?
    Channel
        .fromPath(samplesheet_file, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample.trim(), row.era?.trim()?.toLowerCase() ?: 'modern') }
        .filter { sample_id, era -> sample_id }
    :
    Channel
        .fromPath(legacy_samples_file, checkIfExists: true)
        .splitText()
        .map { it.trim() }
        .filter { it.length() > 0 }
        .map { sample_id -> tuple(sample_id, 'modern') }

modern_samples_ch = sample_metadata_ch
    .filter { sample_id, era -> era != params.historical_era }
    .map { sample_id, era -> tuple(sample_id, resolve_reads(sample_id)) }

historical_samples_ch = sample_metadata_ch
    .filter { sample_id, era -> era == params.historical_era }
    .map { sample_id, era ->
        def reads = resolve_reads(sample_id)
        tuple(sample_id, reads[0], reads[1])
    }

// --- Processes ---

process REPEAT_MODELER {
    tag "De novo repeat discovery on ${ref_fasta.baseName}"
    label 'process_high' // Modeler usually requires significantly more RAM/CPU
    publishDir "${params.outdir}/data/reference/modeler", mode: 'copy'

    input:
    path ref_fasta

    output:
    path "consensi.fa", emit: model_library
    path "families.stk", emit: stockholm_library, optional: true

    script:
    """
    # Step A: Build the sequence database for modeling
    BuildDatabase -name species_db ${ref_fasta}

    # Step B: Run de novo classification modeling
    RepeatModeler -database species_db -pa ${task.cpus}

    # Step C: Rescue the output files if they got left in the subfolder
    if [ ! -f consensi.fa ]; then
        echo "Outputs not found in root. Searching inside RM subdirectories..."
        cp RM_*/consensi.fa . 2>/dev/null || true
        cp RM_*/families.stk . 2>/dev/null || true
    fi
    """
}

process EXTRACT_CPG {
    tag "Extracting CpG sites from ${ref_fasta.baseName}"
    label 'process_low'
    module 'container_env:python3'
    publishDir "${params.outdir}/data/reference", mode: 'copy'

    input:
    path ref_fasta

    output:
    path "cpg_sites.bed", emit: cpg_bed

    script:
    """
    crun python3 -c '
    import sys

    ref_fasta = "${ref_fasta}"
    out_bed = "cpg_sites.bed"

    with open(ref_fasta, "r") as f_in, open(out_bed, "w") as f_out:
        chr_name = ""
        seq = []

        def process_seq(header, sequence):
            full_seq = "".join(sequence).upper()
            pos = 0
            while True:
                pos = full_seq.find("CG", pos)
                if pos == -1:
                    break
                f_out.write(f"{header}\\t{pos}\\t{pos+2}\\n")
                pos += 1

        for line in f_in:
            line = line.strip()
            if line.startswith(">"):
                if chr_name:
                    process_seq(chr_name, seq)
                chr_name = line.split()[0][1:]
                seq = []
            else:
                seq.append(line)

        if chr_name:
            process_seq(chr_name, seq)
    '
    """
}

process REPEAT_MASKER {
    tag "Masking ${ref_fasta.baseName}"
    label 'process_medium'
    publishDir "${params.outdir}/data/reference", mode: 'copy'

    input:
    path ref_fasta
    path repeat_library // This accepts the 'consensi.fa' file from REPEAT_MODELER
    path cpg_bed
    
    output:
    path "${ref_fasta.baseName}.combined_mask.bed", emit: mask_bed
    path "${ref_fasta.baseName}.cleaned.regions",     emit: regions

    script:
    """
    # Step A: Scan assembly against the custom de novo library
    RepeatMasker \
      -lib ${repeat_library} \
      -gff \
      -pa ${task.cpus} \
      ${ref_fasta}
    
    # Step B: Standardize structural output coordinates to a 0-based BED layout
    # Grabs chromosome, start (1-based adjusted to 0), and end coordinates
    if ls *.out.gff 1> /dev/null 2>&1; then
        awk 'OFS="\\t" {if (\$1 !~ /^#/) print \$1, \$4-1, \$5}' *.out.gff > repeats.bed
    else
        echo "Warning: No GFF file produced (zero repeats found). Generating empty track."
        touch repeats.bed
    fi

    # Step C: Combine repeats BED with Python-generated CpG BED
    cat repeats.bed ${cpg_bed} | sort -k1,1 -k2,2n > combined_sorted.bed
    
    # Fast native custom merge line replacing 'bedtools merge' to avoid container dependency errors
    awk 'OFS="\\t" { \
        if (NR==1) {chr=\$1; start=\$2; end=\$3; next} \
        if (\$1==chr && \$2<=end) { if (\$3>end) end=\$3 } \
        else { print chr, start, end; chr=\$1; start=\$2; end=\$3 } \
    } END { print chr, start, end }' combined_sorted.bed > ${ref_fasta.baseName}.combined_mask.bed

    # Generate reference filter tracking indexes for downstream parsing (ANGSD/BCFtools)
    cut -f1 ${ref_fasta.baseName}.combined_mask.bed | uniq | awk '{print \$0 ":"}' > ${ref_fasta.baseName}.cleaned.regions
    """
}


process PREP_REFERENCE_REPEAT {
    publishDir "${params.outdir}/data/reference", mode: 'copy'

    input:
    path ref
    path repeat_bed

    output:
    path "${params.reference_prefix}.repma.angsd.txt", emit: sites
    path "${params.reference_prefix}.repma.angsd.txt.idx", emit: sites_idx
    path "${params.reference_prefix}.repma.angsd.txt.bin", emit: sites_bin
    path "${params.reference_prefix}.regions", emit: regions
    path "${params.reference_prefix}.chrs", emit: chrs

    script:
    """
    awk '{print \$1"\\t"\$2+1"\\t"\$3}' ${repeat_bed} > ${params.reference_prefix}.repma.angsd.txt
    angsd sites index ${params.reference_prefix}.repma.angsd.txt
    cut -f1 ${params.reference_prefix}.repma.angsd.txt | awk '!seen[\$0]++' | awk '{print \$0 ":"}' > ${params.reference_prefix}.regions
    cut -f1 ${params.reference_prefix}.repma.angsd.txt | sort | uniq > ${params.reference_prefix}.chrs
    """
}

process FASTP_MERGED {
    // takes paired-end raw sequencing reads, cleans them up (quality filtering, adapter trimming, poly-G tail removal), 
    // and attempts to merge overlapping forward and reverse reads into single, longer reads.
     tag "$sample_id"
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: "*.{html,json}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.fastq.gz"), emit: merged
    tuple val(sample_id), path("${sample_id}_R1_unmerged.fq.gz"), path("${sample_id}_R2_unmerged.fq.gz"), emit: unmerged
    path "*.json"
    path "*.html"

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} -p -c --trim_poly_g --merge --merged_out=${sample_id}_trimmed_merged.fastq.gz -o ${sample_id}_R1_unmerged.fq.gz -O ${sample_id}_R2_unmerged.fq.gz -h ${sample_id}_fastp_report.html -j ${sample_id}_fastp_report.json -R "${sample_id}" -w ${task.cpus} -l 30 --overlap_diff_limit 1 --overlap_len_require 11
    """
}

process SPLIT_MERGED {
    // processes a FASTQ file and splits any read that is >trim_len right down the middle, 
    // turning a single long read into two "pseudo-paired" reads. 
    // Any read that is already shorter than the target length is left completely untouched.
    tag "$sample_id"
    publishDir "${params.outdir}/data/fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(merged_fq), val(trim_len)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.L${trim_len}.fastq.gz")

    script:
    """
    bash ${params.split_script} ${merged_fq} ${sample_id}_trimmed_merged.L${trim_len}.fastq.gz ${trim_len}
    """
}

process QC_MERGED {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(merged_fq)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    # Specify which java to use, so that it can find libfreetype.so.6
    echo "[DEBUG QC_MERGED] CONDA_PREFIX = \$CONDA_PREFIX"
    echo "[DEBUG QC_MERGED] Java path = \$CONDA_PREFIX/bin/java"
    java -version 2>&1 | head -n 1 | sed 's/^/[DEBUG QC_MERGED] Java version: /'
    fastqc --java "\$CONDA_PREFIX/bin/java" -o . -t ${task.cpus} --extract ${merged_fq}
    """
}

process SEQTK_TRIM {
    // takes unmerged paired-end reads (the ones that didn't overlap enough to be combined) 
    // and  truncates them down to a maximum length (trim_len) using seqtk    
    tag "$sample_id"
    publishDir "${params.outdir}/data/fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2), val(trim_len)

    output:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.L${trim_len}.fastq.gz"), path("${sample_id}_R2_trimmed.L${trim_len}.fastq.gz")

    script:
    """
    seqtk trimfq -L ${trim_len} ${r1} | gzip > ${sample_id}_R1_trimmed.L${trim_len}.fastq.gz
    seqtk trimfq -L ${trim_len} ${r2} | gzip > ${sample_id}_R2_trimmed.L${trim_len}.fastq.gz
    """
}

process QC_UNMERGED {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    # Specify which java to use, so that it can find libfreetype.so.6
    echo "[DEBUG QC_UNMERGED] CONDA_PREFIX = \$CONDA_PREFIX"
    echo "[DEBUG QC_UNMERGED] Java path = \$CONDA_PREFIX/bin/java"
    java -version 2>&1 | head -n 1 | sed 's/^/[DEBUG QC_UNMERGED] Java version: /'
    fastqc --java "\$CONDA_PREFIX/bin/java" -o . -t ${task.cpus} --extract ${r1}
    fastqc --java "\$CONDA_PREFIX/bin/java" -o . -t ${task.cpus} --extract ${r2}
    """
}

process CALC_HISTORICAL_TRIMLEN_BAM {
    publishDir "${params.outdir}/stats", mode: 'copy'
    tag "Calculating mapped length"

    input:
    path bams

    output:
    stdout emit: trim_len
    path "historical_trimlength.txt", emit: trim_file

    script:
    """
    # Loop through all BAMs and stream the mapped reads into a single awk process
    for bam in ${bams}; do
        samtools view -F 2308 "\$bam"
    done | awk '
      BEGIN { bases=0; reads=0 }
      {
        # Column 10 is the sequence. Ignore if it is a missing "*"
        if (\$10 != "*") {
          bases += length(\$10)
          reads++
        }
      }
      END {
        if (reads > 0) {
          # Calculate average and round to nearest integer
          printf "%d\\n", bases/reads + 0.5
        } else {
          # Fallback to default if no mapped reads exist
          print ${params.trimlength}
        }
      }' | tee historical_trimlength.txt
    """
}

process MARKDUP_MERGED {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(merged_bam)
    val trimlength

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.L${trimlength}.sorted.rmdup.bam")

    script:
    """
    samtools view -@ ${task.cpus} -h ${merged_bam} | python3 ${params.rmdup_script} | samtools view -b -o ${sample_id}_trimmed_merged.L${trimlength}.sorted.rmdup.bam
    """
}

process MARKDUP_UNMERGED {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(unmerged_bam)
    val trimlength

    output:
    tuple val(sample_id), path("${sample_id}.L${trimlength}.sorted.rmdup.bam")

    script:
    """
    samtools collate -o ${sample_id}.sorted.namecollate.bam ${unmerged_bam}
    samtools fixmate -m ${sample_id}.sorted.namecollate.bam ${sample_id}.sorted.fixmate.bam
    samtools sort -o ${sample_id}.sorted.fixmate.positionsort.bam ${sample_id}.sorted.fixmate.bam
    samtools markdup -r ${sample_id}.sorted.fixmate.positionsort.bam ${sample_id}.L${trimlength}.sorted.rmdup.bam

    # Cleanup
    rm ${sample_id}.sorted.namecollate.bam ${sample_id}.sorted.fixmate.bam ${sample_id}.sorted.fixmate.positionsort.bam
    """
}

process MERGE_BAMS {
    tag "$sample_name"

    input:
    tuple val(sample_name), path(bams)
    val trimlength

    output:
    tuple val(sample_name), path("${sample_name}.merged.L${trimlength}.bam"), path("${sample_name}.merged.L${trimlength}.bam.bai")

    script:
    """
    samtools merge -@ ${task.cpus} ${sample_name}.merged.L${trimlength}.bam ${bams}
    samtools index ${sample_name}.merged.L${trimlength}.bam
    """
}

process INDEL_REALN {
    tag "$sample_name"
    
    input:
    tuple val(sample_name), path(bam), path(bai)
    path ref_bundle
    val trimlength

    output:
    tuple val(sample_name), path("${sample_name}.merged.L${trimlength}.realn.bam")

    script:
    // ref_bundle[0] is the fasta file
    """
    # Target Creator
    java -jar /usr/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator -R ${ref_bundle[0]} -I ${bam} -o ${sample_name}.realn_targets.list -nt ${task.cpus}

    # Indel Realigner
    java -jar /usr/GenomeAnalysisTK.jar \
        -T IndelRealigner -R ${ref_bundle[0]} -I ${bam} -targetIntervals ${sample_name}.realn_targets.list -o ${sample_name}.merged.L${trimlength}.realn.bam
    """
}

process INDEX_REALIGNED {
    tag "$sample_name"
    publishDir "${params.outdir}/data/bam", mode: 'copy'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path(bam), path("${bam}.bai")

    script:
    """
    samtools index ${bam}
    """
}

process BAM_QC {
    tag "$sample_name"
    publishDir "${params.outdir}/depth", mode: 'copy'

    input:
    tuple val(sample_name), path(bam), path(bai)
    val trimlength
    path bed_file

    output:
    path "${sample_name}.Q25.L${trimlength}.bam.dpstats.txt"

    script:
    """
    # calculate depth
    samtools depth -a -Q 30 -q 25 -b ${bed_file} ${bam} > ${sample_name}.Q25.bam.dp
    awk '{sum+=\$3} END { print sum/NR }' ${sample_name}.Q25.bam.dp > ${sample_name}.Q25.L${trimlength}.bam.dpstats.txt
    rm ${sample_name}.Q25.bam.dp
    """
}

process MAPDAMAGE {
    tag "$sample_name"
    // Send BAM and BAI files to the rescaled bam folder
    publishDir "${params.outdir}/data/bam_rescaled", mode: 'copy', pattern: "*.bam*"
    // Send the PDF and TXT stats files to the mapdamage folder
    publishDir "${params.outdir}/mapdamage", mode: 'copy', pattern: "*.{pdf,txt}"    
    
    module 'container_env:mapdamage2'

    input:
    tuple val(sample_name), path(bam), path(bai)
    path ref_bundle

    output:
    tuple val(sample_name), path("${sample_name}.rescaled.bam"), path("${sample_name}.rescaled.bam.bai"), emit: rescaled_indexed
    path "*.pdf", emit: plots, optional: true
    path "*.txt", emit: stats, optional: true

    script:
    def input_prefix = bam.name.replaceAll(/\.bam$/, '')
    """
    # Run mapDamage for damage assessment and rescaling    
    crun mapDamage -i ${bam} -r ${ref_bundle[0]} -d mapd_output --rescale --merge-reference-sequences

    # Copy rescaled BAM to standard output name and index it
    if [ -f "mapd_output/${input_prefix}.rescaled.bam" ]; then
        cp mapd_output/${input_prefix}.rescaled.bam ${sample_name}.rescaled.bam
        samtools index ${sample_name}.rescaled.bam
    fi
    
    # Loop through PDFs, prepend the sample name, and copy to root
    for pdf in mapd_output/*.pdf; do
        if [ -f "\$pdf" ]; then
            cp "\$pdf" "${sample_name}_\$(basename "\$pdf")"
        fi
    done

    # Loop through TXTs, prepend the sample name, and copy to root
    for txt in mapd_output/*.txt; do
        if [ -f "\$txt" ]; then
            cp "\$txt" "${sample_name}_\$(basename "\$txt")"
        fi
    done
    """
}

process AMBER_PREP {
    tag "$sample_name"
    
    input:
    tuple val(sample_name), path(bam), path(bai)
    path ref_bundle

    output:
    tuple val(sample_name), path("amber_input.txt"), path("${sample_name}.tags.MQ25.bam")

    script:
    """
    # Adding MD tags and filtering for MQ25
    samtools calmd -b ${bam} ${ref_bundle[0]} > ${sample_name}.tags.bam
    samtools view -bq 25 ${sample_name}.tags.bam > ${sample_name}.tags.MQ25.bam

    # Create input file for AMBER
    echo -e "${sample_name}\\t${sample_name}.tags.MQ25.bam" > amber_input.txt
    """
}

process AMBER {
    tag "$sample_name"
    publishDir "${params.outdir}/amber", mode: 'copy'

    input:
    tuple val(sample_name), path(amber_input), path(bam_file)

    output:
    tuple val(sample_name), path("${sample_name}.amber_MQ25.pdf"), path("${sample_name}.amber_MQ25.txt")

    script:
    """
    python3 ${params.amber_script} --bamfiles ${amber_input} --output ${sample_name}.amber_MQ25
    """
}

workflow MODERN_PIPELINE {
    // this takes raw reads as input
    take:
    modern_samples
    trim_length_ch
    mapper_choice_ch
    ref_bundle_ch
    bed_ch

    main:
    // 1. Run the initial FASTP_MERGED
    modern_fastp = FASTP_MERGED(modern_samples)

    // 2. Handle Merged Reads (Split them if they are longer than the trim length)
    modern_merged_split = SPLIT_MERGED(modern_fastp.merged.combine(trim_length_ch))

    // 3. Handle Unmerged Reads (Truncate them directly with SEQTK)
    modern_unmerged_trimmed = SEQTK_TRIM(modern_fastp.unmerged.combine(trim_length_ch))

    // 4. QC steps
    QC_MERGED(modern_merged_split)
    QC_UNMERGED(modern_unmerged_trimmed)

    // 5. Combine with the mapper choice
    modern_merged_for_mapping = modern_merged_split.combine(mapper_choice_ch)
    modern_unmerged_for_mapping = modern_unmerged_trimmed.combine(mapper_choice_ch)
    
    // 6. Mapping
    modern_bwa_merged = BWA_MERGED_PASS1(modern_merged_for_mapping, params.bam_q, trim_length_ch, ref_bundle_ch)
    modern_bwa_unmerged = BWA_UNMERGED_PASS1(modern_unmerged_for_mapping, params.bam_q, trim_length_ch, ref_bundle_ch)

    // 7. Mark duplicates
    modern_markdup_merged = MARKDUP_MERGED(modern_bwa_merged, trim_length_ch)
    modern_markdup_unmerged = MARKDUP_UNMERGED(modern_bwa_unmerged, trim_length_ch)

    // 8. Merge BAMs
    modern_markdup_merged
        .mix(modern_markdup_unmerged)
        .map { id, bam ->
            def sample_name = id.split('_')[0]
            [sample_name, bam]
        }
        .groupTuple()
        .set { modern_bams_to_merge }
    modern_merge_bams = MERGE_BAMS(modern_bams_to_merge, trim_length_ch)

    // 9. Indel realignment and indexing
    modern_realn = INDEL_REALN(modern_merge_bams, ref_bundle_ch, trim_length_ch)
    modern_indexed = INDEX_REALIGNED(modern_realn)
    BAM_QC(modern_indexed, trim_length_ch, bed_ch)

    // 10. Optional AMBER analysis
    if (params.run_modern_amber) {
        AMBER_PREP(modern_indexed, ref_bundle_ch)
        AMBER(AMBER_PREP.out)
    }

    emit:
    indexed = modern_indexed
}

workflow HISTORICAL_PIPELINE {
    // This takes mapped reads as input. The QC and mapping is done in the main workflow.
    take:
    historical_bwa_merged
    historical_bwa_unmerged
    ref_bundle_ch
    trim_length_ch
    bed_ch

    main:
    // 1. Mark duplicates
    historical_markdup_merged = MARKDUP_MERGED(historical_bwa_merged, trim_length_ch)
    historical_markdup_unmerged = MARKDUP_UNMERGED(historical_bwa_unmerged, trim_length_ch)

    // 2. Merge BAMs
    historical_markdup_merged
        .mix(historical_markdup_unmerged)
        .map { id, bam -> [id.split('_')[0], bam] }
        .groupTuple()
        .set { historical_bams_to_merge }
        
    historical_merge_bams = MERGE_BAMS(historical_bams_to_merge, trim_length_ch)

    // 3. Indel realignment and indexing
    historical_realn = INDEL_REALN(historical_merge_bams, ref_bundle_ch, trim_length_ch)
    historical_indexed = INDEX_REALIGNED(historical_realn)
    BAM_QC(historical_indexed, trim_length_ch, bed_ch)

    // 4. Optional MapDamage rescaling
    historical_rescaled = Channel.empty()
    if (params.run_historical_mapdamage) {
        historical_mapdamage = MAPDAMAGE(historical_indexed, ref_bundle_ch)
        historical_rescaled = historical_mapdamage.rescaled_indexed
    }

    // 5. Optional AMBER analysis
    if (params.run_historical_amber) {
        AMBER_PREP(historical_indexed, ref_bundle_ch)
        AMBER(AMBER_PREP.out)
    }

    emit:
    indexed = historical_indexed
    rescaled = historical_rescaled
}

workflow {
    // Package the reference files into a channel
    ref_bundle_ch = Channel.fromPath( "${params.reference}*" )
        .mix( Channel.fromPath( params.reference.replaceAll(/\.(fa|fasta|fna)$/, '.dict'), checkIfExists: false ) )
        .collect()
        .map { files -> 
            // Safely find the actual sequence file
            def fasta = files.find { it.name =~ /\.(fa|fasta|fna)$/ }
            // Get all other index/dict files, remove duplicates, and sort them
            def others = files.findAll { it.name != fasta.name }.unique { it.name }.sort { it.name }
            // Return a reconstructed list where the fasta is index [0]
            return [fasta] + others
        }

    // channnels for reference files and for the initial mapper to use for historical reads
    fasta_ref_ch = Channel.fromPath(params.reference, checkIfExists: true).first()
    ref_fai_ch  = Channel.fromPath("${params.reference}.fai", checkIfExists: true).first()
    ref_dict_ch = Channel.fromPath(params.reference.replaceAll(/\.(fa|fasta|fna)$/, '.dict'), checkIfExists: true).first()
    pass1_mapper_ch = Channel.value(params.historical_mapper)

    // 1. Conduct optional repeat modeling, repeat-masking and CpG site filtering on the reference genome
    bed_ch = Channel.empty()
    if (params.run_repeatmasking) {
        modeler_output = REPEAT_MODELER(fasta_ref_ch) // resource-intensive step, may require high RAM and CPU
        cpg_output     = EXTRACT_CPG(fasta_ref_ch)        
        masking_output = REPEAT_MASKER(fasta_ref_ch, modeler_output.model_library, cpg_output.cpg_bed)
        bed_ch = masking_output.mask_bed 
        regions_ch  = masking_output.regions
    } else {
        // Ensure the manual BED file exists before proceeding
        bed_file_obj = file(params.bed_file)
        if ( !bed_file_obj.exists() ) {
            error "params.run_repeatmasking is false, but the required manual BED file was not found at: ${params.bed_file}"
        }

        // If repeatmasking is skipped, use the provided bed file for downstream analyses
        bed_ch = Channel.fromPath(params.bed_file, checkIfExists: true).first()
    }
    // then prepare the reference for downstream ANGSD analyses
    PREP_REFERENCE_REPEAT(fasta_ref_ch, bed_ch)

    // 2. FASTP on Historical Reads
    historical_fastp_input = historical_samples_ch.map { id, r1, r2 -> tuple(id, [r1, r2]) }    
    historical_fastp = FASTP_MERGED(historical_fastp_input)
    
    if (params.run_historical_fastqc) {
        QC_MERGED(historical_fastp.merged)
        QC_UNMERGED(historical_fastp.unmerged)
    }

    // 3. Initial historical mapping based on user-chosen mapper
    hist_merged_pass1_in = historical_fastp.merged.combine(pass1_mapper_ch)
    hist_unmerged_pass1_in = historical_fastp.unmerged.combine(pass1_mapper_ch)
    
    hist_bwa_merged_pass1 = BWA_MERGED_PASS1(hist_merged_pass1_in, params.bam_q, params.trimlength, ref_bundle_ch)
    hist_bwa_unmerged_pass1 = BWA_UNMERGED_PASS1(hist_unmerged_pass1_in, params.bam_q, params.trimlength, ref_bundle_ch)

    pass1_bams_ch = hist_bwa_merged_pass1.map { id, bam -> bam }
        .mix(hist_bwa_unmerged_pass1.map { id, bam -> bam })
        .collect()

    // 4. Calculate average mapped length and choose optimal mapper
    hist_trim_output = CALC_HISTORICAL_TRIMLEN_BAM(pass1_bams_ch)
    trim_length_ch = hist_trim_output.trim_len.map { it.trim() as Integer }

    mapper_ch = trim_length_ch.map { trim_len -> 
        def chosen = trim_len <= 80 ? "aln" : "mem"
        log.info """
        ===================================================================
         PASS 1 MAPPING COMPLETE:
         -> User-Specified Mapper             : bwa ${params.historical_mapper}
         -> Calculated Mapped Read Length     : ${trim_len} bp
         -> Optimal Chosen Mapper             : bwa ${chosen}
        ===================================================================
        """.stripIndent()
        return chosen
    }

    // 5. Conditional remapping of historical reads (Pass 2 Routing Gates)
    hist_merged_pass2_in = historical_fastp.merged
        .combine(mapper_ch)
        .filter { id, fq, chosen -> chosen != params.historical_mapper }
        
    hist_unmerged_pass2_in = historical_fastp.unmerged
        .combine(mapper_ch)
        .filter { id, r1, r2, chosen -> chosen != params.historical_mapper }

    hist_bwa_merged_pass2 = BWA_MERGED_PASS2(hist_merged_pass2_in, params.bam_q, trim_length_ch, ref_bundle_ch)
    hist_bwa_unmerged_pass2 = BWA_UNMERGED_PASS2(hist_unmerged_pass2_in, params.bam_q, trim_length_ch, ref_bundle_ch)

    // 6. Merge logic: Route correct BAMs downstream based on condition matching
    final_hist_bwa_merged = hist_bwa_merged_pass1
        .combine(mapper_ch)
        .filter { id, bam, chosen -> chosen == params.historical_mapper }
        .map { id, bam, chosen -> tuple(id, bam) }
        .mix(hist_bwa_merged_pass2)

    final_hist_bwa_unmerged = hist_bwa_unmerged_pass1
        .combine(mapper_ch)
        .filter { id, bam, chosen -> chosen == params.historical_mapper }
        .map { id, bam, chosen -> tuple(id, bam) }
        .mix(hist_bwa_unmerged_pass2)
    
    // 7. Sub-pipelines execution to produce indexed bams and optional rescaled bams for historical samples
    modern_pipeline = MODERN_PIPELINE(modern_samples_ch, trim_length_ch, mapper_ch, ref_bundle_ch, bed_ch)
    historical_pipeline = HISTORICAL_PIPELINE(final_hist_bwa_merged, final_hist_bwa_unmerged, ref_bundle_ch, trim_length_ch, bed_ch)
}