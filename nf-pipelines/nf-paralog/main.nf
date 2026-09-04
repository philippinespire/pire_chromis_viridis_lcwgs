// Malin Pinsky, August 2026
// FASTP trims, removes duplidates, and maps modern reads without quality filters
// Then runs ANGSD HWE and depth tests and ngsParalogs to identify putative paralogs.
// Developed from nf-trim-generode with the assistance of Google Gemini.

nextflow.enable.dsl=2

// --- Default Parameters ---
params.samplesheet = "${projectDir}/inputfiles/samplesheet.csv" // Preferred way to provide sample metadata, including sample IDs and eras (modern or historical).
params.indir        = "${projectDir}/data/symlinks" // Directory where raw FASTQ files are expected.
params.outdir       = "${projectDir}/results" // Directory where all output files will be written.
params.reference    = "/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/data/GCA_051013605.1_ASM5101360v1_genomic_20kb.fna" // Path to reference genome FASTA file.
params.reference_prefix = params.reference.tokenize('/').last().replaceAll(/\.(fa|fasta|fna)$/, '') // Extracts base name of reference.
params.modern_era   = "modern" // label in the "era" column of the samplesheet that identifies modern individuals
params.rmdup_script = "${projectDir}/scripts/samremovedup.py"
params.ngsparalog_bin = "/archive/carpenterlab/pire/softwares/ngsParalog/ngsParalog" // path to ngsParalog binary

def resolve_reads = { String sample_id ->
    def r1 = file("${params.indir}/${sample_id}_R1.fastq.gz")
    def r2 = file("${params.indir}/${sample_id}_R2.fastq.gz")

    if( !r1.exists() ) error "R1 file not found for ${sample_id}: ${r1}"
    if( !r2.exists() ) error "R2 file not found for ${sample_id}: ${r2}"

    return [r1, r2]
}


// --- Input Channels ---
// Sample metadata is taken from a CSV when available.
samplesheet_file = file(params.samplesheet)

sample_metadata_ch = Channel
    .fromPath(samplesheet_file, checkIfExists: true)
    .splitCsv(header: true)
    .map { row -> tuple(row.sample.trim(), row.era?.trim()?.toLowerCase() ?: 'modern') }
    .filter { sample_id, era -> sample_id }

modern_samples_ch = sample_metadata_ch
    .filter { sample_id, era -> era == params.modern_era }
    .map { sample_id, era -> tuple(sample_id, resolve_reads(sample_id)) }

// --- Processes ---

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


process MARKDUP_MERGED {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(merged_bam)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.sorted.rmdup.bam")

    script:
    """
    samtools view -@ ${task.cpus} -h ${merged_bam} | python3 ${params.rmdup_script} | samtools view -b -o ${sample_id}_trimmed_merged.sorted.rmdup.bam
    """
}

process MARKDUP_UNMERGED {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(unmerged_bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.rmdup.bam")

    script:
    """
    samtools collate -o ${sample_id}.sorted.namecollate.bam ${unmerged_bam}
    samtools fixmate -m ${sample_id}.sorted.namecollate.bam ${sample_id}.sorted.fixmate.bam
    samtools sort -o ${sample_id}.sorted.fixmate.positionsort.bam ${sample_id}.sorted.fixmate.bam
    samtools markdup -r ${sample_id}.sorted.fixmate.positionsort.bam ${sample_id}.sorted.rmdup.bam

    # Cleanup
    rm ${sample_id}.sorted.namecollate.bam ${sample_id}.sorted.fixmate.bam ${sample_id}.sorted.fixmate.positionsort.bam
    """
}

process MERGE_BAMS {
    tag "$sample_name"

    input:
    tuple val(sample_name), path(bams)

    output:
    tuple val(sample_name), path("${sample_name}.merged.bam"), path("${sample_name}.merged.bam.bai")

    script:
    """
    samtools merge -@ ${task.cpus} ${sample_name}.merged.bam ${bams}
    samtools index ${sample_name}.merged.bam
    """
}

process INDEL_REALN {
    tag "$sample_name"
    
    input:
    tuple val(sample_name), path(bam), path(bai)
    path ref_bundle

    output:
    tuple val(sample_name), path("${sample_name}.merged.realn.bam")

    script:
    // ref_bundle[0] is the fasta file
    """
    # Target Creator
    java -jar /usr/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator -R ${ref_bundle[0]} -I ${bam} -o ${sample_name}.realn_targets.list -nt ${task.cpus}

    # Indel Realigner
    java -jar /usr/GenomeAnalysisTK.jar \
        -T IndelRealigner -R ${ref_bundle[0]} -I ${bam} -targetIntervals ${sample_name}.realn_targets.list -o ${sample_name}.merged.realn.bam
    """
}

process INDEX_REALIGNED {
    tag "$sample_name"
    publishDir "${params.outdir}/large_data/bam", mode: 'copy'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path(bam), path("${bam}.bai")

    script:
    """
    samtools index ${bam}
    """
}

process BWA_MERGED {
    // map the merged reads with no quality filter
    tag "$sample_id"

    input:
    tuple val(sample_id), path(merged_fq)
    path reference_bundle // a bundle of reference files, including the fasta and index files. fasta is first.

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.sorted.bam")

    script:
    def fields = sample_id.split('_')
	def name = fields[0]
    def lib  = fields.size() > 1 ? fields[1] : 'lib1'
    def rg   = fields.size() > 2 ? fields[2] : sample_id
    
    """
    bwa mem -R "@RG\\tID:${rg}\\tSM:${name}\\tPL:ILLUMINA\\tLB:${name}_${lib}\\tPU:${rg}" \
        -t ${task.cpus} ${reference_bundle[0]} ${merged_fq} \
        | samtools view -q 0 -F 1028 -@ ${task.cpus} -bSh - \
        | samtools sort -m 4G -o ${sample_id}_trimmed_merged.sorted.bam -T ${sample_id}.sorting -@ ${task.cpus} -
    """

}

process BWA_UNMERGED {
    // map the unmerged reads with no quality filter
    tag "$sample_id"

    input:
    tuple val(sample_id), path(r1), path(r2)
    path reference_bundle // a bundle of reference files, including the fasta and index files. fasta is first.

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
    def fields = sample_id.split('_')
	def name = fields[0]
    def lib  = fields.size() > 1 ? fields[1] : 'lib1'
    def rg   = fields.size() > 2 ? fields[2] : sample_id

	"""
    bwa mem -R "@RG\\tID:${rg}\\tSM:${name}\\tPL:ILLUMINA\\tLB:${name}_${lib}\\tPU:${rg}" \
        -t ${task.cpus} ${reference_bundle[0]} ${r1} ${r2} \
        | samtools view -q 0 -F 1028 -@ ${task.cpus} -bSh - \
        | samtools sort -m 4G -o ${sample_id}.sorted.bam -T ${sample_id}.sorting -@ ${task.cpus} -
    """
}

process ANGSD_HWE_DEPTH {
    tag "ANGSD_Analysis"
	// Large output: Raw ANGSD calculations routed to large_data/angsd
    publishDir "${params.outdir}/large_data/angsd", mode: 'copy', pattern: "cvi_angsd.*"
    // Lightweight BED file routed to top-level angsd folder
    publishDir "${params.outdir}/angsd", mode: 'copy', pattern: "hwe_excess_het.bed"

    input:
    path bams
    path bais
    path ref_bundle

    output:
    path "cvi_angsd.*"
    path "hwe_excess_het.bed"

    script:
    """
	ls *.bam > bam.filelist

    # 1. Total genome size from reference index (.fai)
    genome_size=\$(awk '{s += \$2} END {print s}' ${ref_bundle[0]}.fai)

    # 2. Total mapped reads across all BAMs (instant execution via idxstats)
    total_reads=0
    for b in *.bam; do
        r=\$(samtools idxstats \$b | awk '{s += \$3} END {print s}')
        total_reads=\$((total_reads + r))
    done

    # 3. Estimate mean total population depth (assuming ~130bp average read length)
    mean_depth=\$(awk -v tr="\$total_reads" -v gs="\$genome_size" 'BEGIN { print (tr * 130) / gs }')

    # 4. Set maxDepth as 5x mean total depth (with a safety floor of 500 for low-coverage runs)
    max_depth=\$(awk -v md="\$mean_depth" 'BEGIN { m = int(md * 5); print (m > 500 ? m : 500) }')

    echo "Estimated Mean Population Depth: \${mean_depth}x"
    echo "Setting ANGSD -setMaxDepth to: \${max_depth}"

    angsd -bam bam.filelist -ref ${ref_bundle[0]} -out cvi_angsd \
        -minMapQ 25 \
        -minQ 30 \
        -remove_bads 1 \
        -uniqueOnly 1 \
        -baq 1 \
        -C 50 \
        -minInd ${min_ind} \
        -setMaxDepth \$max_depth \
        -skipTriallelic 0 \
        -doHWE 1 -doCounts 1 -doDepth 1 -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -P ${task.cpus}

    # Extract sites with excess heterozygosity (p-value < 1e-4) into 1-based BED
    zcat cvi_angsd.hwe.gz | awk 'NR>1 && \$9 < 1e-4 {print \$1 "\t" \$2 "\t" \$2}' > hwe_excess_het.bed
    """
}

process GENERATE_CONTIG_BEDS {
    tag "Generate_Contig_BEDs"

    input:
    path ref_bundle

    output:
    path "*.bed", emit: contig_beds

    script:
    """
    # Extract each scaffold/contig from the FASTA index (.fai) into a full-length BED (0 to length)
	if [ ! -f *.fai ]; then
        samtools faidx ${ref_bundle[0]}
    fi
    fai_file=\$(ls *.fai)
	awk 'OFS="\t" { print \$1, 0, \$2 > (\$1 ".bed"); close(\$1 ".bed") }' \$fai_file    
	"""
}

process NGSPARALOG_CALCLR {
    tag "$contig_bed.baseName"

    input:
    path contig_bed
    path bams
    path bais
    path ref_bundle

    output:
    path "${contig_bed.baseName}.lr.txt"

    script:
	def num_bams = bams instanceof List ? bams.size() : 1
    def min_ind  = Math.floor(num_bams / 2) as int    
    """
    ls *.bam > bam.filelist

    samtools mpileup -b bam.filelist -f ${ref_bundle[0]} -l ${contig_bed} -q 0 -Q 0 | \
        ${params.ngsparalog_bin} calcLR -infile - -outfile - -minQ 20 -minind ${min_ind} > ${contig_bed.baseName}.lr.txt
    """
}

process COMBINE_PARALOGS {
    tag "Combine_Paralogs"
	// Large output: Full likelihood-ratio table sent to large_data/paralogs
    publishDir "${params.outdir}/large_data/paralogs", mode: 'copy', pattern: "cvi_ngsparalog.lr.txt"
    // Lightweight filtered BED file sent to paralogs folder
    publishDir "${params.outdir}/paralogs", mode: 'copy', pattern: "ngsparalog_sites.bed"
    
    input:
    path lr_files

    output:
    path "cvi_ngsparalog.lr.txt"
    path "ngsparalog_sites.bed"

    script:
    """
	cat *.lr.txt > cvi_ngsparalog.lr.txt
    awk '\$3 > 10 {print \$1 "\t" \$2 "\t" \$2}' cvi_ngsparalog.lr.txt > ngsparalog_sites.bed    """
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

    // 1. Run the initial FASTP_MERGED
    modern_fastp = FASTP_MERGED(modern_samples_ch)
    
    // 2. Mapping
    modern_bwa_merged = BWA_MERGED(modern_fastp.merged, ref_bundle_ch)
    modern_bwa_unmerged = BWA_UNMERGED(modern_fastp.unmerged, ref_bundle_ch)

    // 3. Mark duplicates
    modern_markdup_merged = MARKDUP_MERGED(modern_bwa_merged)
    modern_markdup_unmerged = MARKDUP_UNMERGED(modern_bwa_unmerged)

    // 4. Merge BAMs
    modern_markdup_merged
        .mix(modern_markdup_unmerged)
        .map { id, bam ->
            def sample_name = id.split('_')[0]
            [sample_name, bam]
        }
        .groupTuple()
        .set { modern_bams_to_merge }
    modern_merge_bams = MERGE_BAMS(modern_bams_to_merge)

    // 5. Indel realignment and indexing
    modern_realn = INDEL_REALN(modern_merge_bams, ref_bundle_ch)
    modern_indexed = INDEX_REALIGNED(modern_realn)
    
    // Collect all final BAMs and BAIs into single list emissions
    all_bams_ch = modern_indexed.map { sample, bam, bai -> bam }.collect()
    all_bais_ch = modern_indexed.map { sample, bam, bai -> bai }.collect()

    // 6. Run ANGSD HWE & Depth filtering
    angsd_out = ANGSD_HWE_DEPTH(all_bams_ch, all_bais_ch, ref_bundle_ch)

    // 7. Generate candidate BED file by scaffold/contig
	contig_beds_ch = GENERATE_CONTIG_BEDS(ref_bundle_ch).contig_beds.flatten()

    // 8. Run ngsParalog across full scaffolds concurrently
    lr_chunks = NGSPARALOG_CALCLR(contig_beds_ch, all_bams_ch, all_bais_ch, ref_bundle_ch)

    // 9. Combine into final output files
    COMBINE_PARALOGS(lr_chunks.collect())

}