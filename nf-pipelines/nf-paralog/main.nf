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
params.min_ind_ratio = 0.5 // Require coverage in at least 50% of samples
params.high_depth_quantile = 0.995 // target high depth percentile cutoff
params.lr_quantile = 0.999  // Target percentile cutoff (e.g., 0.999 = top 0.1% highest LR sites)
params.rmdup_script = "${projectDir}/scripts/samremovedup.py"
params.ngsparalog_bin = "/archive/carpenterlab/pire/softwares/ngsParalog/ngsParalog" // path to ngsParalog binary
params.run_duphmm    = true          // Toggle dupHMM execution
params.duphmm_script = "/archive/carpenterlab/pire/softwares/ngsParalog/dupHMM.R"     // Path to dupHMM.R (or place inside pipeline bin/)
params.duphmm_emit   = 0              // 0 = LR only, 1 = both LR and Coverage

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
    debug true  // Streams standard output directly into .nextflow.log and console for capturing the mean depth calc

	// Large output: Raw ANGSD calculations routed to large_data/angsd
    publishDir "${params.outdir}/large_data/angsd", mode: 'copy', pattern: "cvi_angsd.*"
    // Lightweight BED file routed to top-level angsd folder
    publishDir "${params.outdir}/angsd", mode: 'copy', pattern: "*.bed"

    input:
    path bams
    path bais
    path ref_bundle

    output:
    path "cvi_angsd.*"
    path "hwe_excess_het.bed"
    path "high_depth_regions.bed"

    script:
    def num_bams = bams instanceof List ? bams.size() : 1
    def min_ind  = Math.max(1, Math.floor(num_bams * params.min_ind_ratio) as int)

    // 1. Logs min_ind to .nextflow.log
    log.info "[ANGSD_HWE_DEPTH] Total BAMs: ${num_bams} | Calculated minInd: ${min_ind}"

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

    # 3. Empirically calculate average read length by sampling 1,000 reads per BAM
    avg_read_len=\$(for b in *.bam; do samtools view "\$b" | head -n 1000; done | awk '{sum += length(\$10); count++} END {if (count > 0) print sum/count; else print 150}')

    # 4. Estimate mean total population depth using empirical read length
    mean_depth=\$(awk -v tr="\$total_reads" -v rl="\$avg_read_len" -v gs="\$genome_size" 'BEGIN { print (tr * rl) / gs }')

    # 5. Set maxDepth as 5x mean total depth (with a safety floor of 500 for low-coverage runs)
    # max_depth is an input filter that prevents extreme coverage spikes from distorting ANGSD calculations
    max_depth=\$(awk -v md="\$mean_depth" 'BEGIN { m = int(md * 5); print (m > 500 ? m : 500) }')

    echo "Calculated Average Read Length: \${avg_read_len} bp"
    echo "Estimated Mean Population Depth: \${mean_depth}x"
    echo "Setting ANGSD -setMaxDepth to: \${max_depth}"

    # Run ANGSD to generate cvi_angsd.depthGlobal histogram
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

    # Extract sites with excess heterozygosity (p-value < 1e-3) and F<0 into 1-based BED
    zcat cvi_angsd.hwe.gz | awk 'NR>1 && \$7 < 0 && \$9 < 1e-3 {print \$1 "\t" \$2 - 1 "\t" \$2}' > hwe_excess_het.bed

    # Calculate Total Depth Cutoff at Target Quantile (e.g., 99.5th percentile) from ANGSD Histogram
    # DEPTH_CUTOFF is a diagnostic threshold used post-analysis to extract and output those specific high-depth regions into a BED file
    DEPTH_CUTOFF=\$(Rscript -e '
        counts <- scan("cvi_angsd.depthGlobal", quiet = TRUE)
        depths <- 0:(length(counts) - 1)
        cdf <- cumsum(counts) / sum(counts)
        cutoff <- depths[min(which(cdf >= ${params.high_depth_quantile}))]
        cat(cutoff)
    ')

    echo "High depth cutoff (${params.high_depth_quantile} quantile): \$DEPTH_CUTOFF total depth"

    # Stream genome depths and write contiguous regions exceeding the quantile cutoff
    samtools depth -q 25 -Q 30 -f bam.filelist | awk -v cutoff="\$DEPTH_CUTOFF" '
    {
        sum = 0
        for (i = 3; i <= NF; i++) sum += \$i
        if (sum >= cutoff) {
            chrom = \$1
            start = \$2 - 1
            end   = \$2

            if (chrom == prev_chrom && start == prev_end) {
                prev_end = end
            } else {
                if (prev_chrom != "") print prev_chrom "\t" prev_start "\t" prev_end
                prev_chrom = chrom
                prev_start = start
                prev_end   = end
            }
        }
    }
    END {
        if (prev_chrom != "") print prev_chrom "\t" prev_start "\t" prev_end
    }' > high_depth_regions.bed
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
    # Ensure reference index exists
    if [ ! -f *.fai ]; then
        samtools faidx ${ref_bundle[0]}
    fi
    fai_file=\$(ls *.fai)

    # Slice scaffolds into 30 Mb BED windows
    awk -v chunk=30000000 'OFS="\t" {
        chrom = \$1
        len = \$2
        start = 0
        while (start < len) {
            end = start + chunk
            if (end > len) end = len
            filename = sprintf("%s_%d_%d.bed", chrom, start, end)
            print chrom, start, end > filename
            close(filename)
            start = end
        }
    }' \$fai_file
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
    def min_ind  = Math.max(1, Math.floor(num_bams * params.min_ind_ratio) as int)
    """
    set -e -o pipefail # fail immediately if any command fails in the pipeline
    # Increase stack size to avoid potential segmentation faults in ngsParalog
    ulimit -s unlimited

    ls *.bam > bam.filelist

    # Uses mpileup -d 1000 to truncate extreme coverage spikes at repetitive loci before data reaches ngsParalog
    # to avoid excessive memory usage and potential crashes
    # Redirects ngsParalog standard error (2> ngsParalog_calcLR.log) so verbose site warnings 
    # write to a file rather than overflowing the Bash process buffer.
    samtools mpileup -d 1000 -b bam.filelist -f ${ref_bundle[0]} -l ${contig_bed} -q 0 -Q 0 | \
        ${params.ngsparalog_bin} calcLR \
            -infile - \
            -outfile ${contig_bed.baseName}.lr.txt \
            -minQ 20 \
            -minind ${min_ind} \
            -allow_overwrite 1 2> ngsParalog_calcLR.log
    """
}

process COMBINE_PARALOGS {
    tag "Combine_Paralogs"
	// Combined likelihood ratio table and BED file routed to large_data/paralogs
    publishDir "${params.outdir}/large_data/paralogs", mode: 'copy'
    
    input:
    path lr_files

    output:
    path "cvi_ngsparalog.lr.txt", emit: lr_txt
    path "ngsparalog_sites.bed", emit: bed
    path "ngsparalog_threshold.txt", emit: threshold

    script:
    """
    # 1. Concatenate all chunk outputs and sort numerically by chromosome and position
    cat *.lr.txt | sort -k1,1 -k2,2n > cvi_ngsparalog.lr.txt
    
    # 2. Compute dynamic LR cutoff at the targeted quantile via stride sampling
    SAMPLE_COUNT=\$(awk 'NR % 1000 == 0' cvi_ngsparalog.lr.txt | wc -l)
    if [ "\$SAMPLE_COUNT" -gt 0 ]; then
        SAMPLER="awk 'NR % 1000 == 0 {print \$3}' cvi_ngsparalog.lr.txt"
    else
        SAMPLER="awk '{print \$3}' cvi_ngsparalog.lr.txt"
    fi

    THRESHOLD=\$(eval \$SAMPLER | \
        Rscript -e '
            x <- scan("stdin", quiet = TRUE)
            cutoff <- quantile(x, probs = ${params.lr_quantile}, na.rm = TRUE)
            cat(sprintf("%.4f", cutoff))
        ')

    echo "Calculated LR threshold (${params.lr_quantile} quantile): \$THRESHOLD" > ngsparalog_threshold.txt

    # 3. Filter using calculated dynamic threshold and merge contiguous coordinates
    awk -v threshold="\$THRESHOLD" '
    \$3 > threshold {
        chrom = \$1
        start = \$2 - 1
        end   = \$2

        if (chrom == prev_chrom && start == prev_end) {
            prev_end = end
        } else {
            if (prev_chrom != "") {
                print prev_chrom "\\t" prev_start "\\t" prev_end
            }
            prev_chrom = chrom
            prev_start = start
            prev_end   = end
        }
    }
    END {
        if (prev_chrom != "") {
            print prev_chrom "\\t" prev_start "\\t" prev_end
        }
    }' cvi_ngsparalog.lr.txt > ngsparalog_sites.bed    
    """
}

process CALCULATE_AVG_DEPTH {
    tag "Avg_Depth"
    publishDir "${params.outdir}/large_data/paralogs", mode: 'copy'

    input:
    path lr_file
    path bams
    path bais
    path ref_bundle

    output:
    path "cvi_avg_depth.tsv", emit: avg_depth

    script:
    def num_bams = bams instanceof List ? bams.size() : 1
    """
    ls *.bam > bam.filelist

    # 1. Create a temporary BED file of targeted LR sites
    awk '{print \$1 "\t" \$2 - 1 "\t" \$2}' ${lr_file} > lr_sites.bed

    # 2. Run mpileup restricted strictly to those sites
    samtools mpileup -d 1000 -b bam.filelist -f ${ref_bundle[0]} -l lr_sites.bed -q 0 -Q 0 | \
        awk -v n_samples=${num_bams} '{print \$1 "\t" \$2 "\t" \$4 / n_samples}' > depth_raw.tsv

    # 3. Map depths back to match lr_file line-for-line (default 0 depth for uncovered sites)
    awk 'NR==FNR { depth[\$1 "\t" \$2] = \$3; next } { d = ((\$1 "\t" \$2) in depth) ? depth[\$1 "\t" \$2] : 0; print \$1 "\t" \$2 "\t" d }' depth_raw.tsv ${lr_file} > cvi_avg_depth.tsv

    rm -f lr_sites.bed depth_raw.tsv
    """
}

// Train HMM parameters globally across all genome sites
process DUPHMM_TRAIN {
    tag "dupHMM_Train"

    input:
    path lr_file
    path cov_file

    output:
    path "genome_trained.par", emit: param_file

    script:
    def cov_arg = (params.duphmm_emit == 1 && cov_file) ? "--covfile=${cov_file}" : ""
    """
    if [ "${params.duphmm_emit}" -eq 1 ]; then
        LR_LINES=\$(wc -l < ${lr_file})
        COV_LINES=\$(wc -l < ${cov_file})

        if [ "\$LR_LINES" -ne "\$COV_LINES" ]; then
            echo "ERROR: Site mismatch detected between LR (\$LR_LINES sites) and Depth (\$COV_LINES sites) files!" >&2
            exit 1
        fi
    fi

    Rscript ${params.duphmm_script} \
        --lrfile=${lr_file} \
        ${cov_arg} \
        --emit=${params.duphmm_emit} \
        --paramOnly=1 \
        --outfile=genome_trained
    """
}

process SPLIT_SCAFFOLDS {
    tag "Split_Scaffolds"

    input:
    path lr_file
    path cov_file

    output:
    path "*.lr.txt", emit: lr_files
    path "*.cov.txt", emit: cov_files, optional: true

    script:
    """
    awk '{file = \$1 ".lr.txt"; print >> file; close(file)}' ${lr_file}
    if [ "${params.duphmm_emit}" -eq 1 ]; then
        awk '{file = \$1 ".cov.txt"; print >> file; close(file)}' ${cov_file}
    fi
    """
}

// Infer duplication states scaffold-by-scaffold using trained parameters
process DUPHMM_INFER {
    tag "$scaff"

    input:
    tuple val(scaff), path(scaff_lr), path(scaff_cov)
    path param_file

    output:
    path "${scaff}_duphmm.rf", emit: rf, optional: true

    script:
    def cov_arg = params.duphmm_emit == 1 ? "--covfile=${scaff_cov}" : ""
    """
    if [ \$(wc -l < ${scaff_lr}) -gt 20 ]; then
        Rscript ${params.duphmm_script} \
            --lrfile=${scaff_lr} \
            ${cov_arg} \
            --emit=${params.duphmm_emit} \
            --paramfile=${param_file} \
            --outfile="${scaff}_duphmm"
    fi
    """
}

process COMBINE_DUPHMM {
    tag "Combine_dupHMM"
    publishDir "${params.outdir}/paralogs", mode: 'copy'

    input:
    path rf_files

    output:
    path "ngsparalog_duphmm_regions.bed"

    script:
    """
    touch ngsparalog_duphmm_regions.bed
    if ls *_duphmm.rf 1> /dev/null 2>&1; then
        cat *_duphmm.rf | awk -v OFS="\\t" '(\$2 + 0 == \$2 && \$2 != "") {print \$1, \$2 - 1, \$3, \$4}' > ngsparalog_duphmm_regions.bed
    fi
    """
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

    if (params.run_duphmm) {
        // Calculate per-site average depth from mpileup if it will be used (emit=1)   
        if (params.duphmm_emit == 1) {
            CALCULATE_AVG_DEPTH(COMBINE_PARALOGS.out.lr_txt, all_bams_ch, all_bais_ch, ref_bundle_ch)
            ch_cov = CALCULATE_AVG_DEPTH.out.avg_depth
        } else {
            ch_cov = Channel.of([])
        }

        // Train global HMM parameters
        DUPHMM_TRAIN(COMBINE_PARALOGS.out.lr_txt, ch_cov)

        // Split whole genome files into scaffold-level channels
        SPLIT_SCAFFOLDS(COMBINE_PARALOGS.out.lr_txt, ch_cov)

        // Match LR and Depth chunks by scaffold name into a unified tuple channel
        ch_lr = SPLIT_SCAFFOLDS.out.lr_files.flatten().map { f -> tuple(f.name.replace('.lr.txt', ''), f) }

        if (params.duphmm_emit == 1) {
            ch_cov_split = SPLIT_SCAFFOLDS.out.cov_files.flatten().map { f -> tuple(f.name.replace('.cov.txt', ''), f) }
            ch_duphmm_inputs = ch_lr.join(ch_cov_split)
        } else {
            ch_duphmm_inputs = ch_lr.map { scaff, lr -> tuple(scaff, lr, []) }
        }

        // Infer duplication states concurrently per scaffold
        DUPHMM_INFER(ch_duphmm_inputs, DUPHMM_TRAIN.out.param_file)

        // Collect all scaffold output files into the final BED
        COMBINE_DUPHMM(DUPHMM_INFER.out.rf.collect().ifEmpty([]))
    }
}