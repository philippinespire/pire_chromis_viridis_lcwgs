process BWA_MERGED {
    // map the merged reads
    tag "$sample_id"

    input:
    tuple val(sample_id), path(merged_fq), val(mapper)
    val bam_q
    val trimlength
    path reference_bundle // a bundle of reference files, including the fasta and index files. fasta is first.

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_merged.L${trimlength}.sorted.bam")

    script:
    def fields = sample_id.split('_')
    def name = fields[0]
    def lib  = fields[1]
    def rg   = fields[2]
    
    if (mapper == "aln") {
        """
        bwa aln -l 16500 -n 0.01 -o 2  -t ${task.cpus} ${reference_bundle[0]} ${merged_fq} > ${sample_id}.sai
        
        bwa samse -r "@RG\\tID:${rg}\\tSM:${name}\\tPL:ILLUMINA\\tLB:${name}_${lib}\\tPU:${rg}" \
            ${reference_bundle[0]} ${sample_id}.sai ${merged_fq} \
            | samtools view -q${bam_q} -F 4 -@ ${task.cpus} -bSh - \
            | samtools sort -m 4G -o ${sample_id}_trimmed_merged.L${trimlength}.sorted.bam -T ${sample_id}.sorting -@ ${task.cpus} -
        """
    } else {
        """
        bwa mem -R "@RG\\tID:${rg}\\tSM:${name}\\tPL:ILLUMINA\\tLB:${name}_${lib}\\tPU:${rg}" \
            -t ${task.cpus} ${reference_bundle[0]} ${merged_fq} \
            | samtools view -q${bam_q} -F 4 -@ ${task.cpus} -bSh - \
            | samtools sort -m 4G -o ${sample_id}_trimmed_merged.L${trimlength}.sorted.bam -T ${sample_id}.sorting -@ ${task.cpus} -
        """
    }
}

process BWA_UNMERGED {
    // map the unmerged reads
    tag "$sample_id"

    input:
    tuple val(sample_id), path(r1), path(r2), val(mapper)
    val bam_q
    val trimlength
    path reference_bundle // a bundle of reference files, including the fasta and index files. fasta is first.

    output:
    tuple val(sample_id), path("${sample_id}.L${trimlength}.sorted.bam")

    script:
    def fields = sample_id.split('_')
    def name = fields[0]
    def lib  = fields[1]
    def rg   = fields[2]

    if (mapper == "aln") {
        """
        bwa aln -l 16500 -n 0.01 -o 2  -t ${task.cpus} ${reference_bundle[0]} ${r1} > ${sample_id}_R1.sai
        bwa aln -l 16500 -n 0.01 -o 2  -t ${task.cpus} ${reference_bundle[0]} ${r2} > ${sample_id}_R2.sai

        bwa sampe \
            -r "@RG\\tID:${rg}\\tSM:${name}\\tPL:ILLUMINA\\tLB:${name}_${lib}\\tPU:${rg}" \
            ${reference_bundle[0]} ${sample_id}_R1.sai ${sample_id}_R2.sai ${r1} ${r2} \
            | samtools view -q${bam_q} -F 4 -@ ${task.cpus} -bSh - \
            | samtools sort -m 4G -o ${sample_id}.L${trimlength}.sorted.bam -T ${sample_id}.sorting -@ ${task.cpus} -
        """
    } else {
        """
        bwa mem -R "@RG\\tID:${rg}\\tSM:${name}\\tPL:ILLUMINA\\tLB:${name}_${lib}\\tPU:${rg}" \
            -t ${task.cpus} ${reference_bundle[0]} ${r1} ${r2} \
            | samtools view -q${bam_q} -F 4 -@ ${task.cpus} -bSh - \
            | samtools sort -m 4G -o ${sample_id}.L${trimlength}.sorted.bam -T ${sample_id}.sorting -@ ${task.cpus} -
        """
    }
}