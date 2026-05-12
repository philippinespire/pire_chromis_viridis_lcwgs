nextflow.enable.dsl=2

// --- Default Parameters ---
params.samples_file = "${projectDir}/inputfiles/fastq_filenames.txt"
params.indir        = "${projectDir}/data/symlinks"
params.outdir       = "${projectDir}/results" 
params.reference    = "${projectDir}/data/reference/reference.ssl.Cvi20k_rename.fasta"
params.bed_file     = "${projectDir}/data/reference/reference.ssl.Cvi20k_rename.repma.bed"
params.split_script = "${projectDir}/scripts/split_reads.sh"
params.rmdup_script = "${projectDir}/scripts/samremovedup.py"
params.amber_script = "/home/mdehasqu/TOOLS/AMBER/AMBER" // Ignore this.
params.bwa_threads  = 4
params.bam_q        = 1 // Mapping quality. Currently set to 1 simply to remove unmapped reads. 
params.trimlength   = 121 // length of historical reads
params.run_mapdamage = false // change to true to run mapDamage (optional, can be time-consuming)

// --- Input Channel ---
// Reads the file line by line (e.g., TzoCMta031_1_22CVWFLT3L3)
// Checks if R1/R2 exist, and emits a tuple [sample_id, [r1, r2]]
Channel
    .fromPath(params.samples_file)
    .splitText()
    .map { it.trim() }
    .filter { it.length() > 0 }
    .map { sample_id ->
        def r1 = file("${params.indir}/${sample_id}_R1.fastq.gz")
        def r2 = file("${params.indir}/${sample_id}_R2.fastq.gz")
        
        if( !r1.exists() ) error "R1 file not found for ${sample_id}: ${r1}"
        if( !r2.exists() ) error "R2 file not found for ${sample_id}: ${r2}"
        
        return [ sample_id, [r1, r2] ]
    }
    .set { raw_reads_ch }

// --- Processes ---

process FASTP_MERGED {
    tag "$sample_id"
    publishDir "${params.outdir}/data/stats", mode: 'copy', pattern: "*.{html,json}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.fastq.gz"), emit: merged
    tuple val(sample_id), path("${sample_id}_R1_unmerged.fq.gz"), path("${sample_id}_R2_unmerged.fq.gz"), emit: unmerged
    path "*.json"
    path "*.html"

    script:
    """
    fastp \
      -i ${reads[0]} -I ${reads[1]} \
      -p -c --trim_poly_g --merge \
      --merged_out=${sample_id}_trimmed_merged.fastq.gz \
      -o ${sample_id}_R1_unmerged.fq.gz -O ${sample_id}_R2_unmerged.fq.gz \
      -h ${sample_id}_fastp_report.html -j ${sample_id}_fastp_report.json \
      -R "${sample_id}" -w ${task.cpus} -l 30 --overlap_diff_limit 1 --overlap_len_require 11
    """
}

process SPLIT_MERGED {
    tag "$sample_id"
    publishDir "${params.outdir}/data/fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(merged_fq)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.L${params.trimlength}.fastq.gz")

    script:
    """
    bash ${params.split_script} ${merged_fq} ${sample_id}_trimmed_merged.L${params.trimlength}.fastq.gz ${params.trimlength}
    """
}

process QC_MERGED {
    tag "$sample_id"
    publishDir "${params.outdir}/data/stats", mode: 'copy'

    input:
    tuple val(sample_id), path(merged_fq)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc -o . -t ${task.cpus} --extract ${merged_fq}
    """
}

process FASTP_UNMERGED {
    tag "$sample_id"
    publishDir "${params.outdir}/data/stats", mode: 'copy', pattern: "*.{html,json}"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fq.gz"), path("${sample_id}_R2.trimmed.fq.gz")
    path "*.json"
    path "*.html"

    script:
    """
    fastp \
      -i ${r1} -I ${r2} \
      -p -c --trim_poly_g \
      -o ${sample_id}_R1.trimmed.fq.gz -O ${sample_id}_R2.trimmed.fq.gz \
      -h ${sample_id}_unmerged_fastp_report.html \
      -j ${sample_id}_unmerged_fastp_report.json \
      -R "${sample_id}" -w ${task.cpus} -l 30
    """
}

process SEQTK_TRIM {
    tag "$sample_id"
    publishDir "${params.outdir}/data/fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_R1_trimmed.L${params.trimlength}.fastq.gz"), path("${sample_id}_R2_trimmed.L${params.trimlength}.fastq.gz")

    script:
    """
    seqtk trimfq -L ${params.trimlength} ${r1} | gzip > ${sample_id}_R1_trimmed.L${params.trimlength}.fastq.gz
    seqtk trimfq -L ${params.trimlength} ${r2} | gzip > ${sample_id}_R2_trimmed.L${params.trimlength}.fastq.gz
    """
}

process QC_UNMERGED {
    tag "$sample_id"
    publishDir "${params.outdir}/data/stats", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc -o . -t ${task.cpus} --extract ${r1}
    fastqc -o . -t ${task.cpus} --extract ${r2}
    """
}

process BWAALN_MERGED {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(merged_fq)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.L${params.trimlength}.sorted.bam")

    script:
    def fields = sample_id.split('_')
    def name = fields[0]
    def lib  = fields[1]
    def rg   = fields[2]
    
    """
    bwa aln -l 16500 -n 0.01 -o 2  -t ${task.cpus} ${params.reference} ${merged_fq} > ${sample_id}.sai
    //bwa mem -M -t ${task.cpus} -R $(cat {input.rg}) {input.ref} {input.fastq_mod_R1} {input.fastq_mod_R2} // under development

    
    bwa samse -r "@RG\\tID:${rg}\\tSM:${name}\\tPL:ILLUMINA\\tLB:${name}_${lib}\\tPU:${rg}" \
        ${params.reference} ${sample_id}.sai ${merged_fq} \
        | samtools view -q${params.bam_q} -F 4 -@ ${task.cpus} -bSh - \
        | samtools sort -m 4G -o ${sample_id}_trimmed_merged.L${params.trimlength}.sorted.bam -T ${sample_id}.sorting -@ ${task.cpus} -
    """
}

process BWAALN_UNMERGED {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.L${params.trimlength}.sorted.bam")

    script:
    def fields = sample_id.split('_')
    def name = fields[0]
    def lib  = fields[1]
    def rg   = fields[2]

    """
    bwa aln -l 16500 -n 0.01 -o 2  -t ${task.cpus} ${params.reference} ${r1} > ${sample_id}_R1.sai
    bwa aln -l 16500 -n 0.01 -o 2  -t ${task.cpus} ${params.reference} ${r2} > ${sample_id}_R2.sai

    bwa sampe \
        -r "@RG\\tID:${rg}\\tSM:${name}\\tPL:ILLUMINA\\tLB:${name}_${lib}\\tPU:${rg}" \
        ${params.reference} ${sample_id}_R1.sai ${sample_id}_R2.sai ${r1} ${r2} \
        | samtools view -q${params.bam_q} -F 4 -@ ${task.cpus} -bSh - \
        | samtools sort -m 4G -o ${sample_id}.L${params.trimlength}.sorted.bam -T ${sample_id}.sorting -@ ${task.cpus} -
    """
}

process BWAMEM_MERGED {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(merged_fq)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.L${params.trimlength}.sorted.bam")

    script:
    def fields = sample_id.split('_')
    def name = fields[0]
    def lib  = fields[1]
    def rg   = fields[2]
    
    """
    bwa mem -M -t ${task.cpus} \
        -R "@RG\\tID:${rg}\\tSM:${name}\\tPL:ILLUMINA\\tLB:${name}_${lib}\\tPU:${rg}" \
        ${params.reference} ${merged_fq} \
        | samtools view -q${params.bam_q} -F 4 -@ ${task.cpus} -bSh - \
        | samtools sort -m 4G -o ${sample_id}_trimmed_merged.L${params.trimlength}.sorted.bam -T ${sample_id}.sorting -@ ${task.cpus} -
    """
}

process BWAMEM_UNMERGED {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.L${params.trimlength}.sorted.bam")

    script:
    def fields = sample_id.split('_')
    def name = fields[0]
    def lib  = fields[1]
    def rg   = fields[2]

    """
    bwa mem -M -t ${task.cpus} \
        -R "@RG\\tID:${rg}\\tSM:${name}\\tPL:ILLUMINA\\tLB:${name}_${lib}\\tPU:${rg}" \
        ${params.reference} ${r1} ${r2} \
        | samtools view -q${params.bam_q} -F 4 -@ ${task.cpus} -bSh - \
        | samtools sort -m 4G -o ${sample_id}.L${params.trimlength}.sorted.bam -T ${sample_id}.sorting -@ ${task.cpus} -
    """
}

process MARKDUP_MERGED {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(merged_bam)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.L${params.trimlength}.sorted.rmdup.bam")

    script:
    """
    samtools view -@ ${task.cpus} -h ${merged_bam} | python3 ${params.rmdup_script} | samtools view -b -o ${sample_id}_trimmed_merged.L${params.trimlength}.sorted.rmdup.bam
    """
}

process MARKDUP_UNMERGED {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(unmerged_bam)

    output:
    tuple val(sample_id), path("${sample_id}.L${params.trimlength}.sorted.rmdup.bam")

    script:
    """
    samtools collate -o ${sample_id}.sorted.namecollate.bam ${unmerged_bam}
    samtools fixmate -m ${sample_id}.sorted.namecollate.bam ${sample_id}.sorted.fixmate.bam
    samtools sort -o ${sample_id}.sorted.fixmate.positionsort.bam ${sample_id}.sorted.fixmate.bam
    samtools markdup -r ${sample_id}.sorted.fixmate.positionsort.bam ${sample_id}.L${params.trimlength}.sorted.rmdup.bam
    
    # Cleanup
    rm ${sample_id}.sorted.namecollate.bam ${sample_id}.sorted.fixmate.bam ${sample_id}.sorted.fixmate.positionsort.bam
    """
}

process MERGE_BAMS {
    tag "$sample_name"

    input:
    tuple val(sample_name), path(bams)

    output:
    tuple val(sample_name), path("${sample_name}.merged.L${params.trimlength}.bam"), path("${sample_name}.merged.L${params.trimlength}.bam.bai")

    script:
    """
    samtools merge -@ ${task.cpus} ${sample_name}.merged.L${params.trimlength}.bam ${bams}
    samtools index ${sample_name}.merged.L${params.trimlength}.bam
    """
}

process INDEL_REALN {
    tag "$sample_name"
    
    input:
    tuple val(sample_name), path(bam), path(bai)
    path ref
    path ref_fai
    path ref_dict

    output:
    tuple val(sample_name), path("${sample_name}.merged.L${params.trimlength}.realn.bam")

    script:
    """
    # Target Creator
    java -jar /usr/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator -R ${ref} -I ${bam} -o ${sample_name}.realn_targets.list -nt ${task.cpus}

    # Indel Realigner
    java -jar /usr/GenomeAnalysisTK.jar \
        -T IndelRealigner -R ${ref} -I ${bam} -targetIntervals ${sample_name}.realn_targets.list -o ${sample_name}.merged.L${params.trimlength}.realn.bam
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
    publishDir "${params.outdir}/results/stats", mode: 'copy'

    input:
    tuple val(sample_name), path(bam), path(bai)

    output:
    path "${sample_name}.Q25.L${params.trimlength}.bam.dpstats.txt"

    script:
    """
    # Depth
    samtools depth -a -Q 30 -q 25 -b ${params.bed_file} ${bam} > ${sample_name}.Q25.bam.dp
    awk '{sum+=\$3} END { print sum/NR }' ${sample_name}.Q25.bam.dp > ${sample_name}.Q25.L${params.trimlength}.bam.dpstats.txt
    rm ${sample_name}.Q25.bam.dp
    """
}

process MAPDAMAGE {
    tag "$sample_name"
    publishDir "${params.outdir}/results/stats/mapdamage", mode: 'copy'
    module 'container_env:mapdamage2' // not sure this is needed, but it doesn't hurt to load the module

    input:
    tuple val(sample_name), path(bam), path(bai)
    path ref

    output:
    path("${sample_name}.mapdamage")

    script:
    """
    crun mapDamage -i ${bam} -r ${ref} --folder ${sample_name}.mapdamage --no-stats
    """
}

process AMBER_PREP {
    tag "$sample_name"
    
    input:
    tuple val(sample_name), path(bam), path(bai)

    output:
    tuple val(sample_name), path("amber_input.txt"), path("${sample_name}.tags.MQ25.bam")

    script:
    """
    # Adding MD tags and filtering for MQ25
    samtools calmd -b ${bam} ${params.reference} > ${sample_name}.tags.bam
    samtools view -bq 25 ${sample_name}.tags.bam > ${sample_name}.tags.MQ25.bam
    
    # Create input file for AMBER
    echo -e "${sample_name}\t${sample_name}.tags.MQ25.bam" > amber_input.txt
    """
}

process AMBER {
    tag "$sample_name"
    publishDir "${params.outdir}/results/stats", mode: 'copy'

    input:
    tuple val(sample_name), path(amber_input), path(bam_file)

    output:
    tuple val(sample_name), path("${sample_name}.amber_MQ25.pdf"), path("${sample_name}.amber_MQ25.txt")

    script:
    """
    # Run AMBER using the input file created in the previous step
    python3 ${params.amber_script} --bamfiles ${amber_input} --output ${sample_name}.amber_MQ25
    """
}

workflow {
    // 1. Trimming & Splitting
    FASTP_MERGED(raw_reads_ch)
    SPLIT_MERGED(FASTP_MERGED.out.merged)
    FASTP_UNMERGED(FASTP_MERGED.out.unmerged)
    SEQTK_TRIM(FASTP_UNMERGED.out[0])
    
    // 2. QC (New Split Logic)
    QC_MERGED(SPLIT_MERGED.out)
    QC_UNMERGED(SEQTK_TRIM.out)
    
    // 3. Mapping (Per Lane)
    def merged_mapped_ch
    def unmerged_mapped_ch

    if (params.trimlength <= 80) { // Use BWA-ALN for shorter reads, BWA-MEM for longer reads
        BWAALN_MERGED(SPLIT_MERGED.out)
        BWAALN_UNMERGED(SEQTK_TRIM.out)
        merged_mapped_ch = BWAALN_MERGED.out
        unmerged_mapped_ch = BWAALN_UNMERGED.out
    } else {
        BWAMEM_MERGED(SPLIT_MERGED.out)
        BWAMEM_UNMERGED(SEQTK_TRIM.out)
        merged_mapped_ch = BWAMEM_MERGED.out
        unmerged_mapped_ch = BWAMEM_UNMERGED.out
    }

    // 4. Mark Duplicates (Split Processes)
    MARKDUP_MERGED(merged_mapped_ch)
    MARKDUP_UNMERGED(unmerged_mapped_ch)

    // 5. Merge BAMs (Per Biological Sample)
    // Mix both streams (merged & unmerged), convert Lane ID to Sample Name, then Group
    MARKDUP_MERGED.out
        .mix(MARKDUP_UNMERGED.out)
        .map { id, bam -> 
            def sample_name = id.split('_')[0] // e.g. "TzoCMta031" from "TzoCMta031_1_22CVWFLT3L3"
            return [ sample_name, bam ]
        }
        .groupTuple() // Groups all BAMs (merged & unmerged) for "TzoCMta031"
        .set { bams_to_merge }

    MERGE_BAMS(bams_to_merge)
    
    // 6. Indel Realign & QC
    ref_ch      = Channel.fromPath(params.reference).first()
    ref_fai_ch  = Channel.fromPath("${params.reference}.fai").first()
    ref_dict_ch = Channel.fromPath(params.reference.replaceAll(/\.fasta$/, '.dict')).first()

    INDEL_REALN(MERGE_BAMS.out, ref_ch, ref_fai_ch, ref_dict_ch)
    INDEX_REALIGNED(INDEL_REALN.out)
    
    // Run Depth QC
    BAM_QC(INDEX_REALIGNED.out)

    // Optional mapDamage
    if (params.run_mapdamage) {
        MAPDAMAGE(INDEX_REALIGNED.out, ref_ch)
    }

    // Run AMBER (Prep -> Run)
    // AMBER_PREP(INDEX_REALIGNED.out)
    // AMBER(AMBER_PREP.out)
}