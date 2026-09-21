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
params.duphmm_emit   = 1              // 0 = LR only, 1 = both LR and Coverage

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

// find HWE and depth outliers using ANGSD, also get ready to calculate allele balance at heterozygote sites and proportion heterozygotes at each site
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
    path "cvi_angsd.*",            emit: angsd_files
    path "cvi_angsd.mafs.gz",      emit: mafs_gz
    path "cvi_angsd.pos.gz",       emit: pos_gz
    path "cvi_angsd.hwe.gz",       emit: hwe_gz
    path "hwe_excess_het.bed",     emit: hwe_bed
    path "high_depth_regions.bed", emit: high_depth_bed
    path "cvi_angsd.geno.gz",       emit: geno_gz
    path "cvi_angsd.counts.gz",     emit: counts_gz

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
    echo "Setting ANGSD -setMaxDepth input filter for depthGlobal histogram to \${max_depth} to filter out coverage spikes"

    # Run ANGSD to generate cvi_angsd.depthGlobal histogram, .pos.gz, .mafs.gz, .geno.gz, .counts.gz, and .hwe.gz files
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
        -doHWE 1 -doCounts 1 -dumpCounts 3 -doDepth 1 -doGeno 8 -doPost 1 -GL 1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -P ${task.cpus}

    # Extract sites with excess heterozygosity (p-value < 1e-3) and F<0 into 1-based BED
    zcat cvi_angsd.hwe.gz | awk 'NR>1 && \$7 < 0 && \$9 < 1e-3 {print \$1 "\t" \$2 - 1 "\t" \$2}' > hwe_excess_het.bed

    # Calculate Total Depth Cutoff at Target Quantile (e.g., 99.5th percentile) from ANGSD Histogram
    # DEPTH_CUTOFF is a diagnostic threshold used post-analysis to extract and output those specific high-depth regions into a BED file
    DEPTH_CUTOFF=\$(awk -v q="${params.high_depth_quantile}" '{
        total = 0;
        for (i = 1; i <= NF; i++) total += \$i;
        target = total * q;
        cum = 0;
        for (i = 1; i <= NF; i++) {
            cum += \$i;
            if (cum >= target) {
                print (i - 1);
                exit;
            }
        }
    }' cvi_angsd.depthGlobal)

    echo "High depth cutoff (${params.high_depth_quantile} quantile) for flagging problematic high-depth regions: \$DEPTH_CUTOFF total depth"

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
    publishDir "${params.outdir}/paralogs", mode: 'copy', saveAs: { fn -> fn == 'cvi_ngsparalog.lr.txt' ? null : fn }
    publishDir "${params.outdir}/large_data/paralogs", mode: 'copy', pattern: "cvi_ngsparalog.lr.txt"

    input:
    path lr_files

    output:
    path "cvi_ngsparalog.lr.txt",    emit: lr_txt
    path "ngsparalog_threshold.txt", emit: threshold_txt
    path "ngsparalog_sites.bed",     emit: paralog_bed

    script:
    """
    # 1. Concatenate all chunk outputs and sort numerically by chromosome and position
    cat *.lr.txt | sort -k1,1 -k2,2n > cvi_ngsparalog.lr.txt

    # 2. Compute dynamic LR threshold (99.9th percentile) using AWK & sort -g (No R dependency)
    LINE_COUNT=\$(wc -l < cvi_ngsparalog.lr.txt)

    if [ "\$LINE_COUNT" -gt 100000 ]; then
        THRESHOLD=\$(awk 'NR % 100 == 0 {print \$3}' cvi_ngsparalog.lr.txt | sort -g | awk -v q="${params.lr_quantile}" '
            { a[NR] = \$1 }
            END {
                if (NR == 0) { print 0; exit }
                idx = int(NR * q);
                if (idx < 1) idx = 1;
                printf "%.4f", a[idx]
            }')
    else
        THRESHOLD=\$(awk '{print \$3}' cvi_ngsparalog.lr.txt | sort -g | awk -v q="${params.lr_quantile}" '
            { a[NR] = \$1 }
            END {
                if (NR == 0) { print 0; exit }
                idx = int(NR * q);
                if (idx < 1) idx = 1;
                printf "%.4f", a[idx]
            }')
    fi

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
                print prev_chrom "\t" prev_start "\t" prev_end
            }
            prev_chrom = chrom
            prev_start = start
            prev_end   = end
        }
    }
    END {
        if (prev_chrom != "") {
            print prev_chrom "\t" prev_start "\t" prev_end
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
    def cov_arg = (params.duphmm_emit == 1 && cov_file) ? "--covfile=${cov_file}" : "" // Only pass --covfile if emit 1 is active and a coverage file exists
    def rcmd = task.ext.rscript ?: 'Rscript' // allow nextflow.config to override Rscript path if needed (eg, with crun Rscript)
    """
    if [ "${params.duphmm_emit}" -eq 1 ]; then
        LR_LINES=\$(wc -l < ${lr_file})
        COV_LINES=\$(wc -l < ${cov_file})

        if [ "\$LR_LINES" -ne "\$COV_LINES" ]; then
            echo "ERROR: Site mismatch detected between LR (\$LR_LINES sites) and Depth (\$COV_LINES sites) files!" >&2
            exit 1
        fi
    fi

    ${rcmd} ${params.duphmm_script} \
        --lrfile=${lr_file} \
        ${cov_arg} \
        --emit=${params.duphmm_emit} \
        --paramOnly=1 \
        --outfile=genome_trained
    """
}

// split the big lr_file from ngsParalog calcLR apart by scaffold/contig
process SPLIT_SCAFFOLDS {
    tag "Split_Scaffolds"

    input:
    path lr_file
    path cov_file

    output:
    path "split_lr/*.lr.txt", emit: lr_files
    path "split_cov/*.cov.txt", emit: cov_files, optional: true

    script:
    // To guard against limits on how many files can be open at once, 
    // this function will close files after every write if the number of scaffolds exceeds a threshold (250).
    // To speed up writing, however, it will keep files open if it can (useful for many-GB files)
    """
    mkdir -p split_lr split_cov

    split_stream() {
        local infile="\$1"
        local ext="\$2"
        local outdir="\$3"

        awk -v ext="\$ext" -v outdir="\$outdir" '
        BEGIN { max_open = 250 }
        {
            file = outdir "/" \$1 ext

            if (!(file in seen)) {
                seen[file] = 1
                scaff_count++
                if (scaff_count > max_open) {
                    mode = "close_every"
                }
            }

            if (mode == "close_every") {
                print >> file
                close(file)
            } else {
                print > file
            }
        }' "\$infile"
    }

    split_stream "${lr_file}" ".lr.txt" "split_lr"

    if [ "${params.duphmm_emit}" -eq 1 ] && [ -s "${cov_file}" ]; then
        split_stream "${cov_file}" ".cov.txt" "split_cov"
    fi

    true
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
    def cov_arg = params.duphmm_emit == 1 ? "--covfile=${scaff_cov}" : "" // Only pass --covfile if emit 1 is active and a coverage file exists
    def rcmd = task.ext.rscript ?: 'Rscript' // allow nextflow.config to override Rscript path if needed (eg, with crun Rscript)
    """
    if [ \$(wc -l < ${scaff_lr}) -gt 20 ]; then
        ${rcmd} ${params.duphmm_script} \
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

// Extract allele balance and het proportion from ANGSD outputs
// Weights allele counts by individual posterior heterozygosity probabilities to soft-call heterozygous sites without requiring hard genotype cutoffs.
process PARSE_ALLELE_BALANCE {
    tag "Parse_Allele_Balance"
    publishDir "${params.outdir}/large_data/stats", mode: 'copy'

    input:
    path geno_gz
    path counts_gz

    output:
    path "site_allele_balance.tsv", emit: site_stats

    script:
    """
    python3 - << 'EOF' > site_allele_balance.tsv
    import gzip
    from itertools import chain

    print("chr\tpos\tallele_balance\tprop_het")

    with gzip.open("${geno_gz}", "rt") as f_geno, gzip.open("${counts_gz}", "rt") as f_counts:
        # 1. Check f_counts for header
        first_c = f_counts.readline()
        if not first_c:
            exit(0)

        if not first_c.strip().split()[0].lstrip("-").isdigit():
            f_counts_iter = f_counts
        else:
            f_counts_iter = chain([first_c], f_counts)

        # 2. Check f_geno for header
        first_g = f_geno.readline()
        if not first_g:
            exit(0)

        if not first_g.strip().split()[1].lstrip("-").isdigit():
            f_geno_iter = f_geno
        else:
            f_geno_iter = chain([first_g], f_geno)

        # 3. Process matched rows inside the open file context
        for line_g, line_c in zip(f_geno_iter, f_counts_iter):
            parts_g = line_g.strip().split()
            parts_c = line_c.strip().split()

            if not parts_g or not parts_c:
                continue

            chrom = parts_g[0]
            pos   = parts_g[1]

            try:
                posteriors = [float(x) for x in parts_g[2:]]
                counts     = [int(x) for x in parts_c]
            except ValueError:
                continue

            num_samples = len(posteriors) // 3
            if num_samples == 0:
                continue

            p_het_sum = 0.0
            weighted_minor_depth = 0.0
            weighted_total_depth = 0.0

            for i in range(num_samples):
                p_het = posteriors[i * 3 + 1]
                p_het_sum += p_het

                if len(counts) >= (i + 1) * 2:
                    maj_cnt = counts[i * 2]
                    min_cnt = counts[i * 2 + 1]
                    tot_cnt = maj_cnt + min_cnt

                    if tot_cnt > 0:
                        weighted_minor_depth += p_het * min_cnt
                        weighted_total_depth += p_het * tot_cnt

            prop_het = p_het_sum / float(num_samples)

            if weighted_total_depth > 0 and prop_het > 0:
                allele_balance = weighted_minor_depth / weighted_total_depth
                print(f"{chrom}\t{pos}\t{allele_balance:.4f}\t{prop_het:.4f}")
    EOF
    """
}

// Plot ngsParalog Likelihood Ratio (LR)
process PLOT_PARALOG_LR {
    tag "Plot_Paralog_LR"
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path lr_file
    path threshold_file
    path duphmm_bed
    path ref_bundle

    output:
    path "ngsparalog_lr_manhattan.png"

    script:
    def rcmd = task.ext.rscript ?: 'Rscript'
    """
    awk '{print \$1 "\t" \$2 "\t" \$5}' ${lr_file} | awk '\$2 ~ /^[0-9]+\$/' > lr_parsed.tsv

    ${rcmd} -e '
    options(bitmapType = "cairo")
    library(ggplot2)

    # 1. Standardize chromosome order & offsets from reference index (.fai)
    fai_file  <- list.files(pattern="\\\\.fai\$")[1]
    fai       <- read.table(fai_file, header=FALSE, stringsAsFactors=FALSE)
    chr_order <- fai\$V1
    chr_lens  <- setNames(as.numeric(fai\$V2), fai\$V1)
    chr_offsets <- setNames(c(0, head(cumsum(chr_lens), -1)), fai\$V1)

    # 2. Read threshold and parsed LR data
    thresh_line <- readLines("${threshold_file}")[1]
    threshold   <- as.numeric(sub(".*: ", "", thresh_line))

    df <- read.table("lr_parsed.tsv", header=FALSE, col.names=c("chr", "pos", "lr"))
    df <- df[df\$chr %in% chr_order, ]
    df\$chr <- factor(df\$chr, levels=chr_order)

    # 3. Downsample low-LR background while preserving ALL high-LR points
    target_points <- 1000000
    df_high <- df[df\$lr >= threshold, ]
    df_low  <- df[df\$lr < threshold, ]

    max_low <- target_points - nrow(df_high)
    if (max_low > 0 && nrow(df_low) > max_low) {
      df_low <- df_low[sample(nrow(df_low), max_low), ]
    }

    df <- rbind(df_high, df_low)
    df\$cum_pos <- df\$pos + chr_offsets[as.character(df\$chr)]

    # 4. Parse dupHMM regions if available
    dh_file <- "${duphmm_bed}"
    has_duphmm <- file.exists(dh_file) && file.info(dh_file)\$size > 0

    p <- ggplot(df, aes(x=cum_pos, y=lr, color=chr))

    if (has_duphmm) {
      dh <- read.table(dh_file, header=FALSE, stringsAsFactors=FALSE)[, 1:3]
      colnames(dh) <- c("chr", "start", "end")
      dh <- dh[dh\$chr %in% chr_order, ]
      if (nrow(dh) > 0) {
        dh\$xmin <- dh\$start + chr_offsets[as.character(dh\$chr)]
        dh\$xmax <- dh\$end   + chr_offsets[as.character(dh\$chr)]
        p <- p + geom_rect(
          data=dh,
          aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf),
          fill="#E41A1C", alpha=0.3, inherit.aes=FALSE
        )
      }
    }

    p <- p +
      geom_point(alpha=0.5, size=0.4) +
      geom_hline(yintercept=threshold, color="red", linetype="dashed", linewidth=0.6) +
      scale_color_manual(values=rep(c("#2B5C8F", "#D95F02"), length.out=length(chr_order))) +
      labs(
        title="ngsParalog Likelihood Ratio (LR) Manhattan Plot",
        subtitle=if(has_duphmm) "Red shaded regions = dupHMM excluded calls | Red dashed line = LR cutoff" else paste("Red dashed line = cutoff threshold (", round(threshold, 2), ")", sep=""),
        x="Genome Position", 
        y="LR Statistic"
      ) +
      theme_minimal() +
      theme(legend.position="none", panel.grid.minor=element_blank())

    ggsave("ngsparalog_lr_manhattan.png", plot=p, width=12, height=5, dpi=300, type="cairo")
    '
    rm -f lr_parsed.tsv
    """
}

// Plot ANGSD Read Depth
process PLOT_ANGSD_DEPTH {
    tag "Plot_ANGSD_Depth"
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path pos_gz
    path ref_bundle

    output:
    path "angsd_depth_manhattan.png"

    script:
    def rcmd = task.ext.rscript ?: 'Rscript'
    """
    ${rcmd} -e '
    options(bitmapType = "cairo")
    library(ggplot2)

    fai_file  <- list.files(pattern="\\\\.fai\$")[1]
    fai       <- read.table(fai_file, header=FALSE, stringsAsFactors=FALSE)
    chr_order <- fai\$V1
    chr_lens  <- setNames(as.numeric(fai\$V2), fai\$V1)
    chr_offsets <- setNames(c(0, head(cumsum(chr_lens), -1)), fai\$V1)

    df <- read.table(gzfile("${pos_gz}"), header=TRUE)
    colnames(df) <- c("chr", "pos", "depth")
    df <- df[df\$chr %in% chr_order, ]
    df\$chr <- factor(df\$chr, levels=chr_order)

    if (nrow(df) > 1000000) df <- df[sample(nrow(df), 1000000), ]
    df\$cum_pos <- df\$pos + chr_offsets[as.character(df\$chr)]

    p <- ggplot(df, aes(x=cum_pos, y=depth, color=chr)) +
      geom_point(alpha=0.5, size=0.4) +
      scale_color_manual(values=rep(c("#1B9E77", "#7570B3"), length.out=length(chr_order))) +
      labs(title="ANGSD Genome Read Depth Manhattan Plot", x="Genome Position", y="Total Read Depth") +
      theme_minimal() +
      theme(legend.position="none", panel.grid.minor=element_blank())

    ggsave("angsd_depth_manhattan.png", plot=p, width=12, height=5, dpi=300, type="cairo")
    '
    """
}

// Plot ANGSD HWE p-values and Inbreeding Coefficient (F)
process PLOT_ANGSD_HWE {
    tag "Plot_ANGSD_HWE"
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path hwe_gz
    path ref_bundle

    output:
    path "angsd_hwe_manhattan.png"

    script:
    def rcmd = task.ext.rscript ?: 'Rscript'
    """
    fai_file=\$(ls *.fai)
    zcat ${hwe_gz} | awk 'NR>1 {print \$1 "\t" \$2 "\t" \$7 "\t" \$9}' > hwe_parsed.tsv

    ${rcmd} -e '
    options(bitmapType = "cairo")
    library(ggplot2)

    fai_file  <- list.files(pattern="\\\\.fai\$")[1]
    fai       <- read.table(fai_file, header=FALSE, stringsAsFactors=FALSE)
    chr_order <- fai\$V1
    chr_lens  <- setNames(as.numeric(fai\$V2), fai\$V1)
    chr_offsets <- setNames(c(0, head(cumsum(chr_lens), -1)), fai\$V1)

    df <- read.table("hwe_parsed.tsv", header=FALSE, col.names=c("chr", "pos", "F", "pval"))
    df <- df[df\$chr %in% chr_order, ]
    df\$chr <- factor(df\$chr, levels=chr_order)

    if (nrow(df) > 1000000) df <- df[sample(nrow(df), 1000000), ]

    df\$logp <- -log10(ifelse(df\$pval <= 0, 1e-16, df\$pval))
    df\$cum_pos <- df\$pos + chr_offsets[as.character(df\$chr)]

    df_logp <- df[, c("chr", "cum_pos", "logp")]
    colnames(df_logp) <- c("chr", "cum_pos", "value")
    df_logp\$metric <- "-log10(p-value)"

    df_f <- df[, c("chr", "cum_pos", "F")]
    colnames(df_f) <- c("chr", "cum_pos", "value")
    df_f\$metric <- "Inbreeding Coeff (F)"

    df_long <- rbind(df_logp, df_f)
    df_long\$metric <- factor(df_long\$metric, levels=c("-log10(p-value)", "Inbreeding Coeff (F)"))

    p <- ggplot(df_long, aes(x=cum_pos, y=value, color=chr)) +
      geom_point(alpha=0.5, size=0.4) +
      scale_color_manual(values=rep(c("#E66101", "#5E3C99"), length.out=length(chr_order))) +
      facet_wrap(~ metric, ncol=1, scales="free_y", strip.position="left") +
      labs(title="ANGSD HWE Analysis", x="Genome Position", y=NULL) +
      theme_minimal() +
      theme(
        legend.position="none",
        panel.grid.minor=element_blank(),
        strip.placement="outside",
        strip.text=element_text(size=10, face="bold")
      )

    ggsave("angsd_hwe_manhattan.png", plot=p, width=12, height=8, dpi=300, type="cairo")
    '
    rm -f hwe_parsed.tsv
    """
}

// Compare excluded regions across filtering methods (Depth, HWE, Paralog/dupHMM)
process PLOT_FILTER_OVERLAPS {
    tag "Plot_Filter_Overlaps"
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path depth_bed
    path hwe_bed
    path paralog_bed
    val paralog_label
    path ref_bundle

    output:
    path "paralog_filter_overlaps.png"
    path "filter_region_length_histograms.png"
    path "filter_overlaps_manhattan.png"

    script:
    def rcmd = task.ext.rscript ?: 'Rscript'
    """
    cat << 'EOF' > plot_overlaps.R
    options(bitmapType = "cairo")
    library(ggplot2)
    library(grid)

    read_bed <- function(f, label) {
    if (!file.exists(f) || file.info(f)\$size == 0) {
        return(data.frame(chr=character(), start=integer(), end=integer(), width=numeric(), filter=character(), stringsAsFactors=FALSE))
    }
    df <- read.table(f, header=FALSE, stringsAsFactors=FALSE)[, 1:3]
    colnames(df) <- c("chr", "start", "end")
    df\$width <- df\$end - df\$start
    df\$filter <- label
    return(df)
    }

    b_depth   <- read_bed("${depth_bed}", "High Depth")
    b_hwe     <- read_bed("${hwe_bed}", "HWE Excess Het")
    p_label   <- "${paralog_label}"
    b_paralog <- read_bed("${paralog_bed}", p_label)

    # =========================================================================
    # PLOT 1: Summary Bar Plots (Existing Overlap Breakdown)
    # =========================================================================
    all_chrs <- unique(c(b_depth\$chr, b_hwe\$chr, b_paralog\$chr))
    combo_counts <- numeric(8)
    names(combo_counts) <- c("None", "Paralog", "HWE", "HWE+Paralog", "Depth", "Depth+Paralog", "Depth+HWE", "All_3")

    in_regions <- function(pos, bed_df) {
    if (nrow(bed_df) == 0) return(FALSE)
    any(pos >= bed_df\$start & pos <= bed_df\$end)
    }

    for (chr_name in all_chrs) {
    sub_d <- b_depth[b_depth\$chr == chr_name, ]
    sub_h <- b_hwe[b_hwe\$chr == chr_name, ]
    sub_p <- b_paralog[b_paralog\$chr == chr_name, ]

    pts <- sort(unique(c(sub_d\$start, sub_d\$end, sub_h\$start, sub_h\$end, sub_p\$start, sub_p\$end)))
    if (length(pts) < 2) next

    for (i in 1:(length(pts) - 1)) {
        s_i <- pts[i]
        e_i <- pts[i+1]
        len_i <- e_i - s_i
        if (len_i <= 0) next
        mid_i <- (s_i + e_i) / 2

        is_d <- in_regions(mid_i, sub_d)
        is_h <- in_regions(mid_i, sub_h)
        is_p <- in_regions(mid_i, sub_p)

        idx <- 1 + (if (is_d) 4 else 0) + (if (is_h) 2 else 0) + (if (is_p) 1 else 0)
        combo_counts[idx] <- combo_counts[idx] + len_i
    }
    }

    combo_mb <- combo_counts / 1e6
    total_depth   <- sum(combo_mb[c("Depth", "Depth+Paralog", "Depth+HWE", "All_3")])
    total_hwe     <- sum(combo_mb[c("HWE", "HWE+Paralog", "Depth+HWE", "All_3")])
    total_paralog <- sum(combo_mb[c("Paralog", "HWE+Paralog", "Depth+Paralog", "All_3")])

    df_indiv <- data.frame(
    Filter = factor(c("High Depth", "HWE Excess Het", p_label), levels=c("High Depth", "HWE Excess Het", p_label)),
    Mb = c(total_depth, total_hwe, total_paralog)
    )

    df_combos <- data.frame(
    Category = factor(
        c("Paralog Only", "HWE Only", "Depth Only", "HWE + Paralog", "Depth + Paralog", "Depth + HWE", "All 3 Filters"),
        levels=c("Paralog Only", "HWE Only", "Depth Only", "HWE + Paralog", "Depth + Paralog", "Depth + HWE", "All 3 Filters")
    ),
    Mb = c(combo_mb["Paralog"], combo_mb["HWE"], combo_mb["Depth"], combo_mb["HWE+Paralog"], combo_mb["Depth+Paralog"], combo_mb["Depth+HWE"], combo_mb["All_3"])
    )
    levels(df_combos\$Category) <- gsub("Paralog", p_label, levels(df_combos\$Category))

    p1 <- ggplot(df_indiv, aes(x=Filter, y=Mb, fill=Filter)) +
    geom_col(width=0.6, show.legend=FALSE) +
    geom_text(aes(label=sprintf("%.2f Mb", Mb)), vjust=-0.5, size=3.5, fontface="bold") +
    scale_fill_manual(values=c("#377EB8", "#4DAF4A", "#E41A1C")) +
    labs(title="Total Excluded Megabases (Mb) per Filter", x=NULL, y="Genomic Region Size (Mb)") +
    theme_minimal() +
    theme(panel.grid.minor=element_blank(), axis.text.x=element_text(face="bold", size=10))

    p2 <- ggplot(df_combos, aes(x=Category, y=Mb, fill=Category)) +
    geom_col(width=0.6, show.legend=FALSE) +
    geom_text(aes(label=sprintf("%.2f Mb", Mb)), vjust=-0.5, size=3.2) +
    scale_fill_brewer(palette="Set3") +
    labs(title="Filter Intersection Breakdown", x=NULL, y="Genomic Region Size (Mb)") +
    theme_minimal() +
    theme(panel.grid.minor=element_blank(), axis.text.x=element_text(angle=30, hjust=1, face="bold", size=9))

    png("paralog_filter_overlaps.png", width=12, height=5, units="in", res=300, type="cairo")
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(1, 2)))
    print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
    dev.off()

    # =========================================================================
    # PLOT 2: Multi-panel Histograms of Excluded Region Lengths
    # =========================================================================
    df_lens <- rbind(b_depth, b_hwe, b_paralog)
    if (nrow(df_lens) > 0) {
    df_lens\$filter <- factor(df_lens\$filter, levels=c("High Depth", "HWE Excess Het", p_label))

    p_hist <- ggplot(df_lens, aes(x=width, fill=filter)) +
        geom_histogram(bins=50, color="white", linewidth=0.2, alpha=0.85) +
        scale_x_log10(labels=function(x) sprintf("%g bp", x)) +
        facet_wrap(~ filter, ncol=1, scales="free_y") +
        scale_fill_manual(values=c("High Depth"="#377EB8", "HWE Excess Het"="#4DAF4A", "#E41A1C")) +
        labs(
        title="Distribution of Excluded Region Lengths per Filter",
        subtitle="Region length = end - start coordinate (log10 scale)",
        x="Region Length (bp)",
        y="Region Count"
        ) +
        theme_minimal() +
        theme(
        legend.position="none",
        strip.text=element_text(face="bold", size=10),
        panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold", size=13)
        )

    ggsave("filter_region_length_histograms.png", plot=p_hist, width=9, height=7, dpi=300, type="cairo")
    } else {
    png("filter_region_length_histograms.png", width=9, height=7, res=300, type="cairo")
    plot.new()
    text(0.5, 0.5, "No excluded regions found across filters")
    dev.off()
    }

    # =========================================================================
    # PLOT 3: Track-based Manhattan Plot of Excluded Loci & Overlaps
    # =========================================================================
    fai_file  <- list.files(pattern="[.]fai\$")[1]
    fai       <- read.table(fai_file, header=FALSE, stringsAsFactors=FALSE)
    chr_order <- fai\$V1
    chr_lens  <- setNames(as.numeric(fai\$V2), fai\$V1)
    chr_offsets <- setNames(c(0, head(cumsum(chr_lens), -1)), fai\$V1)

    add_offsets <- function(df) {
    if (nrow(df) == 0) return(df)
    df <- df[df\$chr %in% chr_order, ]
    if (nrow(df) == 0) return(df)
    df\$cum_start <- df\$start + chr_offsets[as.character(df\$chr)]
    df\$cum_end   <- df\$end   + chr_offsets[as.character(df\$chr)]
    return(df)
    }

    bd <- add_offsets(b_depth)
    bh <- add_offsets(b_hwe)
    bp <- add_offsets(b_paralog)

    overlap_segments <- list()
    for (chr_name in all_chrs) {
    sub_d <- b_depth[b_depth\$chr == chr_name, ]
    sub_h <- b_hwe[b_hwe\$chr == chr_name, ]
    sub_p <- b_paralog[b_paralog\$chr == chr_name, ]

    pts <- sort(unique(c(sub_d\$start, sub_d\$end, sub_h\$start, sub_h\$end, sub_p\$start, sub_p\$end)))
    if (length(pts) < 2) next

    for (i in 1:(length(pts) - 1)) {
        s_i <- pts[i]
        e_i <- pts[i+1]
        mid_i <- (s_i + e_i) / 2

        hits <- sum(c(in_regions(mid_i, sub_d), in_regions(mid_i, sub_h), in_regions(mid_i, sub_p)))
        if (hits >= 2) {
        overlap_segments[[length(overlap_segments) + 1]] <- data.frame(
            chr=chr_name, start=s_i, end=e_i, width=e_i-s_i, filter="2+ Filter Overlap", stringsAsFactors=FALSE
        )
        }
    }
    }

    b_over <- if (length(overlap_segments) > 0) do.call(rbind, overlap_segments) else data.frame(chr=character(), start=integer(), end=integer(), width=numeric(), filter=character())
    bo <- add_offsets(b_over)

    track_list <- list()
    if (nrow(bp) > 0) { bp\$track <- p_label; bp\$y <- 1; track_list[[1]] <- bp }
    if (nrow(bh) > 0) { bh\$track <- "HWE Excess Het"; bh\$y <- 2; track_list[[2]] <- bh }
    if (nrow(bd) > 0) { bd\$track <- "High Depth"; bd\$y <- 3; track_list[[3]] <- bd }
    if (nrow(bo) > 0) { bo\$track <- "2+ Filter Overlap"; bo\$y <- 4; track_list[[4]] <- bo }

    df_man <- if (length(track_list) > 0) do.call(rbind, track_list) else NULL

    if (!is.null(df_man) && nrow(df_man) > 0) {
    df_man\$track <- factor(df_man\$track, levels=c(p_label, "HWE Excess Het", "High Depth", "2+ Filter Overlap"))

    chr_bounds <- data.frame(
        chr = chr_order,
        start = chr_offsets,
        end = chr_offsets + chr_lens
    )

    p_man <- ggplot() +
        geom_rect(
        data=chr_bounds[seq(1, nrow(chr_bounds), 2), ],
        aes(xmin=start, xmax=end, ymin=0.4, ymax=4.6),
        fill="grey93", alpha=0.5, inherit.aes=FALSE
        ) +
        geom_rect(
        data=df_man,
        aes(xmin=cum_start, xmax=cum_end, ymin=y-0.35, ymax=y+0.35, fill=track),
        color=NA, alpha=0.9
        ) +
        scale_fill_manual(values=c("High Depth"="#377EB8", "HWE Excess Het"="#4DAF4A", "#E41A1C", "2+ Filter Overlap"="#984EA3")) +
        scale_y_continuous(
        breaks=1:4,
        labels=c(p_label, "HWE Excess Het", "High Depth", "2+ Filter Overlap"),
        limits=c(0.4, 4.6)
        ) +
        labs(
        title="Manhattan Plot of Excluded Genomic Loci Across Filters",
        subtitle="Top track highlights problematic regions flagged by two or more filters concurrently",
        x="Genome Position",
        y=NULL
        ) +
        theme_minimal() +
        theme(
        legend.position="none",
        panel.grid.minor=element_blank(),
        panel.grid.major.y=element_blank(),
        axis.text.y=element_text(face="bold", size=10),
        plot.title=element_text(face="bold", size=13)
        )

    ggsave("filter_overlaps_manhattan.png", plot=p_man, width=13, height=6, dpi=300, type="cairo")
    } else {
    png("filter_overlaps_manhattan.png", width=13, height=6, res=300, type="cairo")
    plot.new()
    text(0.5, 0.5, "No excluded loci found across filters")
    dev.off()
    }
    EOF

    ${rcmd} plot_overlaps.R
    """
}

// Plot Allele Balance vs Proportion of Heterozygotes with hexbins
process PLOT_ALLELE_BALANCE {
    tag "Plot_Allele_Balance"
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path site_stats // Tab-delimited file with header: chr, pos, allele_balance, prop_het

    output:
    path "allele_balance_vs_het.png"

    script:
    def rcmd = task.ext.rscript ?: 'Rscript'
    """
    ${rcmd} -e '
    options(bitmapType = "cairo")
    library(ggplot2)

    df <- read.table("${site_stats}", header=TRUE, stringsAsFactors=FALSE)
    colnames(df)[1:4] <- c("chr", "pos", "allele_balance", "prop_het")

    # Filter out missing data and ensure bounds [0, 1]
    df <- df[!is.na(df\$allele_balance) & !is.na(df\$prop_het), ]
    df <- df[df\$allele_balance >= 0 & df\$allele_balance <= 1, ]

    p <- ggplot(df, aes(x=allele_balance, y=prop_het)) +
      geom_bin2d(bins=100) +
      scale_fill_viridis_c(option="magma", trans="log10", name="Site Count") +
      geom_vline(xintercept=0.5, linetype="dashed", color="white", alpha=0.8) +
      scale_x_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1)) +
      scale_y_continuous(breaks=seq(0, 1, 0.1), limits=c(0, 1)) +
      labs(
        title="Allele Balance vs. Proportion Heterozygotes",
        subtitle="Dashed line = Expected diploid 0.5 allele balance",
        x="Allele Balance (Alt / Total Reads at Het Sites)",
        y="Proportion of Heterozygous Individuals"
      ) +
      theme_minimal() +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(face="bold", size=13),
        axis.title = element_text(face="bold")
      )

    png("allele_balance_vs_het.png", width=8, height=6, units="in", res=300, type="cairo")
    print(p)
    dev.off()
    '
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

    // 6. Run ANGSD HWE & Depth filtering, output allele balance and heterozygosity stats
    angsd_out = ANGSD_HWE_DEPTH(all_bams_ch, all_bais_ch, ref_bundle_ch)

    // 6.1. Parse probabilistic allele balance & heterozygosity from ANGSD outputs
    ab_stats_ch = PARSE_ALLELE_BALANCE(angsd_out.geno_gz, angsd_out.counts_gz)

    // 7. Generate candidate BED file by scaffold/contig
	contig_beds_ch = GENERATE_CONTIG_BEDS(ref_bundle_ch).contig_beds.flatten()

    // 8. Run ngsParalog across full scaffolds concurrently
    lr_chunks = NGSPARALOG_CALCLR(contig_beds_ch, all_bams_ch, all_bais_ch, ref_bundle_ch)

    // 9. Combine into final output files
    COMBINE_PARALOGS(lr_chunks.collect())

    // Optional DUPHMM run
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
        DUPHMM_INFER(ch_duphmm_inputs, DUPHMM_TRAIN.out.param_file.first())

        // Collect all scaffold output files into the final BED
        COMBINE_DUPHMM(DUPHMM_INFER.out.rf.collect().ifEmpty([]))
    }

    // --- Generate Diagnostic & Comparison Plots ---
    ch_duphmm_plot = params.run_duphmm ? COMBINE_DUPHMM.out : Channel.of(file("NO_DUPHMM"))

    PLOT_PARALOG_LR(
        COMBINE_PARALOGS.out.lr_txt,
        COMBINE_PARALOGS.out.threshold_txt,
        ch_duphmm_plot,
        ref_bundle_ch
    )

    PLOT_ANGSD_DEPTH(
        angsd_out.pos_gz,
        ref_bundle_ch
    )

    PLOT_ANGSD_HWE(
        angsd_out.hwe_gz,
        ref_bundle_ch
    )

    PLOT_ALLELE_BALANCE(ab_stats_ch.site_stats)

    ch_paralog_bed    = params.run_duphmm ? COMBINE_DUPHMM.out : COMBINE_PARALOGS.out.paralog_bed
    val_paralog_label = params.run_duphmm ? "dupHMM" : "ngsParalog LR"

    PLOT_FILTER_OVERLAPS(
        angsd_out.high_depth_bed,
        angsd_out.hwe_bed,
        ch_paralog_bed,
        val_paralog_label,
        ref_bundle_ch
    )
}