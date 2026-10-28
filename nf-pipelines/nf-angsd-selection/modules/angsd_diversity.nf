process ANGSD_DIVERSITY {
    tag "${region}_${era}"
    publishDir "${params.outdir}/diversity",          mode: 'copy', pattern: "*.{sfs,pestPG}"
    publishDir "${params.outdir}/large_data/thetas", mode: 'copy', pattern: "*.thetas*"
    
    input:
    tuple val(region), val(era), path(saf), path(saf_idx), path(saf_pos)
    val label
    tuple path(sites_pos), path(sites_bin), path(sites_idx)

    output:
    tuple val(region), val(era), val(label), path("${region}_${era}_${label}.sfs"), emit: sfs
    tuple val(region), val(era), val(label), path("${region}_${era}_${label}.pestPG"), emit: pestPG
    path "${region}_${era}_${label}.thetas.idx", emit: thetas_idx
    path "${region}_${era}_${label}.thetas.gz", emit: thetas

    script:
    """
    # 1. Generate Unpruned SFS using neutral callable sites (includes monomorphic loci)
    realSFS ${saf_idx} -sites ${sites_pos} -P ${task.cpus ?: 32} -fold 1 > ${region}_${era}_${label}.sfs

    # 2. Calculate Pi and Theta using the Neutral SFS
    realSFS saf2theta ${saf_idx} -sites ${sites_pos} -sfs ${region}_${era}_${label}.sfs -fold 1 -P ${task.cpus ?: 32} -outname ${region}_${era}_${label}
    thetaStat do_stat ${region}_${era}_${label}.thetas.idx -outnames ${region}_${era}_${label}
    """
}

process PLOT_PI_DIVERSITY {
    tag "Plot pi: ${region}"
    publishDir "${params.outdir}/diversity", mode: 'copy'

    input:
    tuple val(region), val(label), path(hist_pestPG), path(mod_pestPG)
    path script

    output:
    path "pi_historic_vs_modern_${region}_${label}.png" , emit: plot
    path "pi_summary_${region}_${label}.tsv"            , emit: summary

    script:
    """
    Rscript ${script} \
        --historic_file="${hist_pestPG}" \
        --modern_file="${mod_pestPG}" \
        --region="${region}" \
        --out_png="pi_historic_vs_modern_${region}_${label}.png" \
        --out_tsv="pi_summary_${region}_${label}.tsv" \
        --n_boot=10000
    """
}