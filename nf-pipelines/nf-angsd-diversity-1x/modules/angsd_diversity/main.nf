process ANGSD_DIVERSITY {
    tag { "${pop}_${era}" }
    publishDir "${params.outdir}/angsd_pop_theta", mode: 'copy'

    input:
    tuple val(pop), val(era), path(saf), path(saf_idx), path(saf_pos)
    path pruned_files // Staged pruned sites [pos, bin, idx] or [] if ld_prune is false

    output:
    path "${pop}_${era}.sfs", emit: sfs
    path "${pop}_${era}.ld_pruned.sfs", optional: true, emit: ld_pruned_sfs
    path "${pop}_${era}.pestPG", emit: pestPG
    path "${pop}_${era}.thetas.idx", emit: thetas_idx
    path "${pop}_${era}.thetas.gz", emit: thetas

    script:
    """
    module load container_env ngsTools

    # 1. Generate Unpruned SFS (All callable sites, no SNP_pval filter, no LD pruning)
    realSFS ${pop}_${era}.saf.idx -P ${task.cpus ?: 32} -fold 1 > ${pop}_${era}.sfs

    # 2. Calculate Pi and Theta using the UNPRUNED SFS
    realSFS saf2theta ${pop}_${era}.saf.idx -sfs ${pop}_${era}.sfs -fold 1 -P ${task.cpus ?: 32} -outname ${pop}_${era}
    thetaStat do_stat ${pop}_${era}.thetas.idx -outnames ${pop}_${era}

    # 3. Generate LD-Pruned SFS (Only if LD pruning was enabled and pruned_sites.pos exists)
    # Useful for FST and other downstream analyses that require LD-pruned data
    if [ -f "pruned_sites.pos" ]; then
        realSFS ${pop}_${era}.saf.idx -sites pruned_sites.pos -P ${task.cpus ?: 32} -fold 1 > ${pop}_${era}.ld_pruned.sfs
    fi
    """
}