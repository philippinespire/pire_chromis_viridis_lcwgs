process ANGSD_GL_ALL {
    tag "$contigs"

    input:
    val contigs
    path bamlist

    output:
    tuple val(contigs), path("${params.species}.${contigs}.mafs.gz"), emit: mafs
    tuple val(contigs), path("${params.species}.${contigs}.beagle.gz"), emit: beagle

    script:
    """
        angsd -bam ${bamlist} \
        -GL 1 \
        -doGlf 2 \
        -doMajorMinor 4 \
        -doMaf 1 \
        -minMapQ 25 -minQ 30 \
        -SNP_pval 1e-6 \
        -minInd ${params.minind} \
        -uniqueOnly 1 -remove_bads 1 \
        -skipTriallelic	0 \
        -doCounts 1 -doDepth 1 -dumpCounts 1 -setmaxdepth ${params.maxdepth} \
        -noTrans 1 \
        -P 1 \
        -ref ${params.reference} \
        -sites ${params.bed_file} \
        -r ${contigs} \
        -out ${params.species}.${contigs}
    """
}

process ANGSD_COLLECT_OUTPUT {
    publishDir "${params.outdir}/GL", mode: 'copy'

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
    publishDir "${params.outdir}/GL", mode: 'copy'

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