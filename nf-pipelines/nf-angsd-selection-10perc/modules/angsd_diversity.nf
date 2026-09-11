process ANGSD_DIVERSITY {
    tag "${region}_${era}"
    publishDir "${params.outdir}/diversity",          mode: 'copy', pattern: "*.{sfs,pestPG}"
    publishDir "${params.outdir}/large_data/thetas", mode: 'copy', pattern: "*.thetas*"
    
    input:
    tuple val(region), val(era), path(saf), path(saf_idx), path(saf_pos)
    path neutral_sites  // [pos, bin, idx] bundle from PREPARE_DIVERSITY_SITES

    output:
    tuple val(region), val(era), path("${region}_${era}.sfs"), emit: sfs
    tuple val(region), val(era), path("${region}_${era}.pestPG"), emit: pestPG
    path "${region}_${era}.thetas.idx", emit: thetas_idx
    path "${region}_${era}.thetas.gz", emit: thetas

    script:
    """
    # 1. Generate Unpruned SFS using neutral callable sites (includes monomorphic loci)
    realSFS ${saf_idx} -sites ${neutral_sites[0]} -P ${task.cpus ?: 32} -fold 1 > ${region}_${era}.sfs

    # 2. Calculate Pi and Theta using the Neutral SFS
    realSFS saf2theta ${saf_idx} -sites ${neutral_sites[0]} -sfs ${region}_${era}.sfs -fold 1 -P ${task.cpus ?: 32} -outname ${region}_${era}
    thetaStat do_stat ${region}_${era}.thetas.idx -outnames ${region}_${era}
    """
}

process PLOT_PI_DIVERSITY {
    tag "Plot pi: ${region}"
    publishDir "${params.outdir}/diversity", mode: 'copy'

    input:
    tuple val(region), path(hist_pestPG), path(mod_pestPG)
    path script

    output:
    path "pi_historic_vs_modern_${region}.png" , emit: plot
    path "pi_summary_${region}.tsv"            , emit: summary

    script:
    """
    Rscript ${script} \
        --historic_file="${hist_pestPG}" \
        --modern_file="${mod_pestPG}" \
        --region="${region}" \
        --out_png="pi_historic_vs_modern_${region}.png" \
        --out_tsv="pi_summary_${region}.tsv" \
        --n_boot=10000
    """
}