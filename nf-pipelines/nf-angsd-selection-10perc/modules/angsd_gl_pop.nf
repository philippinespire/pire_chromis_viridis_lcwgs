// Identify loci passing minind_frac and maxdepth_mult per population/era subset
// Step 1: Calculate parameters once per pop/era
process CALC_POP_THRESHOLDS {
    tag "${pop}_${era}"

    input:
    tuple val(pop), val(era), path(bamlist)

    output:
    tuple val(pop), val(era), path(bamlist), env(min_ind), env(max_depth), emit: thresholds

    script:
    """
    # 1. Total target length from ANGSD-style BED file (1-based start and end). 
    target_bp=\$(awk 'BEGIN{s=0} {s += \$3 - \$2 + 1} END{print s}' ${params.bed_file})

    total_bases=0
    while read -r bam; do
        n_reads=\$(samtools idxstats "\$bam" | awk '{s+=\$3} END{print s}')
        read_len=\$(samtools view "\$bam" | head -n 1000 | awk 'NR>0 {s+=length(\$10); c++} END{if(c>0) print int(s/c); else print 100}')
        total_bases=\$(awk -v tb="\$total_bases" -v nr="\$n_reads" -v rl="\$read_len" 'BEGIN{print tb + (nr * rl)}')
    done < ${bamlist}

    n_ind=\$(wc -l < ${bamlist})
    min_ind=\$(awk -v n="\$n_ind" -v f="${params.min_ind_frac}" 'BEGIN{print int(n * f + 0.999)}')
    
    mean_combined_depth=\$(awk -v tb="\$total_bases" -v tbp="\$target_bp" 'BEGIN{print tb / tbp}')
    max_depth=\$(awk -v md="\$mean_combined_depth" -v m="${params.max_depth_mult}" 'BEGIN{print int(md * m + 0.999)}')
    """
}

// Step 2: Parallel ANGSD pass per contig
process QC_POP_SITES_CONTIG {
    tag "${pop}_${era}:${contig}"

    input:
    tuple val(pop), val(era), path(bamlist), val(min_ind), val(max_depth), val(contig)
    tuple path(ref), path(ref_fai)
    tuple path(sites_pos), path(sites_bin), path(sites_idx) // Accepts bundled files

    output:
    tuple val(pop), val(era), path("${pop}_${era}_${contig}.pos"), emit: contig_pos

    script:
    """
    angsd -bam ${bamlist} \
        -ref "$ref" \
        -r "$contig" \
        -GL 1 -doMajorMinor 4 -doMaf 1 \
        -doCounts 1 \
        -minMapQ 25 -minQ 30 \
        -minInd ${min_ind} \
        -setMaxDepth ${max_depth} \
        -sites ${sites_pos} \
        -out ${pop}_${era}_${contig}_qc \
        -P ${task.cpus}

    zcat ${pop}_${era}_${contig}_qc.mafs.gz | awk 'NR>1 {print \$1"\t"\$2}' > ${pop}_${era}_${contig}.pos
    """
}

// Step 3: Combine contig files back into single per-population position file
process GATHER_POP_SITES {
    tag "${pop}_${era}"

    input:
    tuple val(pop), val(era), path(pos_files)

    output:
    path "${pop}_${era}.passing.pos", emit: passing_pos

    script:
    """
    # 1. Force C locale for fast byte-level sorting (5-10x speedup)
    export LC_ALL=C

    # 2. Linear K-way merge across pre-sorted contig streams
    sort -m -k1,1 -k2,2n \
        -S ${task.memory ? "${task.memory.toGiga()}G" : "2G"} \
        --parallel=${task.cpus ?: 1} \
        ${pos_files} > ${pop}_${era}.passing.pos
    """
}

// Intersect passing sites across both historical and modern population/era subsets (Venn Diagram Overlap)
process INTERSECT_POP_SITES {
    tag "Loci Intersection Across Populations"
    //publishDir "${params.outdir}/sites", mode: 'copy'

    input:
    path pos_files // Collect of all passing.pos files

    output:
    path "pop_intersect.pos"    , emit: pos
    path "pop_intersect.pos.bin", emit: bin
    path "pop_intersect.pos.idx", emit: idx

    script:
    """
    num_files=\$(ls -1 ${pos_files} | wc -l)

    # Retain ONLY loci that pass QC filters in EVERY historical and modern subset
    sort -m -k1,1 -k2,2n ${pos_files} | uniq -c | awk -v n="\$num_files" '\$1 == n {print \$2"\t"\$3}' > pop_intersect.pos

    angsd sites index pop_intersect.pos
    """
}

// Final Population Execution using the shared intersected site index
// run ANGSD on each population and era combination
process ANGSD_GL_POP {
    tag { "${pop}_${era}" }
    publishDir "${params.outdir}/large_data/saf",  mode: 'copy', pattern: "*.saf*"
    publishDir "${params.outdir}/large_data/mafs", mode: 'copy', pattern: "*.mafs.gz"

    input:
    tuple val(pop), val(era), path(bamlist)
    tuple path(ref), path(ref_fai)
    tuple path(intersect_snps), path(intersect_bin), path(intersect_idx)

    output:
    tuple val(pop), val(era), path("${pop}_${era}.saf.gz"), path("${pop}_${era}.saf.idx"), path("${pop}_${era}.saf.pos.gz"), emit: saf_files
    path "${pop}_${era}.mafs.gz", emit: mafs

    script:
    """
    angsd -bam ${bamlist} \
        -doSaf 1 \
        -GL 1 \
        -doMajorMinor 4 \
        -doMaf 1 \
        -minMapQ 25 -minQ 30 \
        -doCounts 1 -doDepth 1 -dumpCounts 1 \
        -uniqueOnly 1 -remove_bads 1 \
        -P ${task.cpus} \
        -ref "$ref" \
        -anc "$ref" \
        -sites ${intersect_snps} \
        -out ${pop}_${era} \
        -noTrans 1
    """
}
