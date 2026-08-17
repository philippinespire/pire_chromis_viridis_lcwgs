process COLLECT_BAM_ALL {

    publishDir "${params.outdir}/inputfiles", mode: 'copy'

    input:
    val bams

    output:
    path "bamlist.txt"

    script:
    """
    printf '%s\\n' ${bams.join(' ')} | head -c -1> bamlist.txt
    """
}

process COLLECT_BAM_POP {

    tag { "${pop}_${era}" }
    publishDir "${params.outdir}/inputfiles", mode: 'copy'

    input:
    tuple val(pop), val(era), val(bams)

    output:
    tuple val(pop), val(era), path("${pop}.${era}.bamlist.txt")

    script:
    """
    printf '%s\\n' ${bams.join(' ')} | head -c -1 > ${pop}.${era}.bamlist.txt
    """
}