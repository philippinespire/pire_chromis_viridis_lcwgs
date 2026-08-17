process ANGSD_GL_ALL {
    tag "$contigs"

    input:
    val contigs
    path bamlist
    tuple path(sites_pos), path(sites_bin), path(sites_idx) // Accepts bundled files
    tuple path(ref), path(ref_fai)

    output:
    tuple val(contigs), path("${params.species}.${contigs}.mafs.gz"), emit: mafs
    tuple val(contigs), path("${params.species}.${contigs}.beagle.gz"), emit: beagle

    script:
    """
    # 1. Total target length from BED file
    target_bp=\$(awk 'BEGIN{s=0} {s += \$3 - \$2 +1} END{print s}' ${params.bed_file})

    # 2. Calculate total mapped bases across all BAMs in bamlist using index files
    total_bases=0
    while read -r bam; do
        n_reads=\$(samtools idxstats "\$bam" | awk '{s+=\$3} END{print s}')
        read_len=\$(samtools view "\$bam" | head -n 1000 | awk 'NR>0 {s+=length(\$10); c++} END{if(c>0) print int(s/c); else print 100}')
        total_bases=\$(awk -v tb="\$total_bases" -v nr="\$n_reads" -v rl="\$read_len" 'BEGIN{print tb + (nr * rl)}')
    done < ${bamlist}

    # 3. Calculate observed combined depth across all samples and set max_depth cutoff
    n_ind=\$(wc -l < ${bamlist})
    min_ind=\$(awk -v n="\$n_ind" -v f="${params.min_ind_frac}" 'BEGIN{print int(n * f + 0.999)}')
    
    mean_combined_depth=\$(awk -v tb="\$total_bases" -v tbp="\$target_bp" 'BEGIN{print tb / tbp}')
    max_depth=\$(awk -v md="\$mean_combined_depth" -v m="${params.max_depth_mult}" 'BEGIN{print int(md * m + 0.999)}')

    # 4. Run ANGSD
    angsd -bam ${bamlist} \
        -GL 1 \
        -doGlf 2 \
        -doMajorMinor 4 \
        -doMaf 1 \
        -minMapQ 25 -minQ 30 \
        -SNP_pval 1e-6 \
        -minInd \$min_ind \
        -setMaxDepth \$max_depth \
        -uniqueOnly 1 -remove_bads 1 \
        -skipTriallelic 0 \
        -doCounts 1 -doDepth 1 -dumpCounts 1 \
        -noTrans 1 \
        -P ${task.cpus ?: 8} \
        -ref ${ref} \
        -sites ${sites_pos} \
        -r ${contigs} \
        -out ${params.species}.${contigs}
    """
}

process ANGSD_COLLECT_OUTPUT {
    publishDir "${params.outdir}/large_data/beagle", mode: 'copy', pattern: "*.beagle.gz"
    publishDir "${params.outdir}/large_data/mafs",   mode: 'copy', pattern: "*.mafs.gz"
    
    input:
    path mafs
    path beagle

    output:
    path "${params.species}.mafs.gz", emit: all_mafs
    path "${params.species}.beagle.gz", emit: all_beagle

    script:
    """
    # Concatenate beagle: header from first, then all without headers
    zcat ${beagle[0]} | head -n 1 > ${params.species}.beagle
    zcat ${beagle} | grep -v -w marker >> ${params.species}.beagle
    gzip ${params.species}.beagle

    # Concatenate mafs: header from first, then all without headers
    zcat ${mafs[0]} | head -n 1 > ${params.species}.mafs
    zcat ${mafs} | grep -v -w chromo >> ${params.species}.mafs
    gzip ${params.species}.mafs
    """
}

process ANGSD_EXTRACT_SITES {
    publishDir "${params.outdir}/large_data/sites", mode: 'copy'
    
    input:
    path all_mafs

    output:
    path "${params.species}.sites.txt" , emit: snps
    path "${params.species}.sites.txt.idx" , emit: snps_idx
    path "${params.species}.sites.txt.bin" , emit: snps_bin
    path "${params.species}.regions" , emit: regions
    path "${params.species}.chrs" , emit: chrs

    script:
    """
    ##Create a SNP list to use in downstream analyses
    zcat ${all_mafs} | cut -f 1,2,3,4 | tail -n +2 > ${params.species}.sites.txt
    angsd sites index ${params.species}.sites.txt

    ## Also make it in regions format for downstream analyses
    cut -f1 "${params.species}.sites.txt" | awk '!seen[\$0]++'| awk '{print \$0 ":"}' > ${params.species}.regions

    ## Lastly, extract a list of chromosomes/LGs/scaffolds for downstream analysis
    cut -f1 ${params.species}.sites.txt | sort | uniq > ${params.species}.chrs
    
    """
}

process INDEX_BED_SITES {
    // Index the bed file to make a list of all callable sites for downstream analyses. 
    // This is useful for calculating diversity statistics across all callable sites, not just SNPs.
    tag "Index BED Sites"
    publishDir "${params.outdir}/large_data/sites", mode: 'copy'

    input:
    path bed_file

    output:
    path "all_sites.pos", emit: snps
    path "all_sites.pos.bin", emit: bin
    path "all_sites.pos.idx", emit: idx
    path "all_sites.regions", emit: regions

    script:
    """
    # Convert 3-column ANGSD BED (1-based start, 1-based end) to 2-column ANGSD POS (1-based pos)
    awk 'BEGIN{OFS="\t"} {for(i=\$2; i<=\$3; i++) print \$1, i}' ${bed_file} > all_sites.pos
    angsd sites index all_sites.pos

    cut -f1 all_sites.pos | awk '!seen[\$0]++' | awk '{print \$0 ":"}' > all_sites.regions
    """
}