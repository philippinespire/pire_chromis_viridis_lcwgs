process ANGSD_DIVERSITY {
    tag { "${pop}_${era}" }
    publishDir "${params.outdir}/angsd_pop_theta", mode: 'copy'

    input:
    tuple val(pop), val(era),path(saf), path(saf_idx), path(pos)

    output:
    path "${pop}_${era}.pestPG", emit: pestPG
    path "${pop}_${era}.thetas.idx", emit: thetas_idx
    path "${pop}_${era}.thetas.gz", emit: thetas
    path "${pop}_${era}.sfs", emit: sfs

    script:
    """
        realSFS ${pop}_${era}.saf.idx -P 8 -fold 1 > ${pop}_${era}.sfs
        realSFS saf2theta ${pop}_${era}.saf.idx -sfs ${pop}_${era}.sfs -fold 1 -P 8 -outname ${pop}_${era}
        thetaStat do_stat ${pop}_${era}.thetas.idx -outnames ${pop}_${era}
    """
}