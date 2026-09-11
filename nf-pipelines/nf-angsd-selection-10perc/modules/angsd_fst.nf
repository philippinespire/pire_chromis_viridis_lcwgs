process ANGSD_FST {
    tag "FST: ${region}"

    publishDir "${params.outdir}/fst",                 mode: 'copy', pattern: "*.{global_fst.txt,png}"
    publishDir "${params.outdir}/large_data/fst",      mode: 'copy', pattern: "*.{fst.gz,fst.idx,2dsfs,windowed_fst.txt}"

    input:
    tuple val(region), path(hist_saf), path(hist_idx), path(hist_pos), path(mod_saf), path(mod_idx), path(mod_pos)
    tuple path(sites_pos), path(sites_bin), path(sites_idx)
    path plot_script

    output:
    path "${region}_hist_vs_mod.2dsfs",          emit: sfs_2d
    path "${region}_hist_vs_mod.global_fst.txt", emit: global_fst
    path "${region}_hist_vs_mod.fst.gz",         emit: fst_gz
    path "${region}_hist_vs_mod.fst.idx",        emit: fst_idx
    path "${region}_windowed_fst.txt",        emit: windows
    path "${region}_fst_manhattan.png",       emit: plot

    script:
    """
    # 1. Estimate 2D-SFS using only LD-pruned neutral sites
    realSFS ${hist_idx} ${mod_idx} \
        -sites ${sites_pos} \
        -P ${task.cpus} > ${region}_hist_vs_mod.2dsfs

    # 2. Index FST per-site
    realSFS fst index ${hist_idx} ${mod_idx} \
        -sfs ${region}_hist_vs_mod.2dsfs \
        -sites ${sites_pos} \
        -fstout ${region}_hist_vs_mod \
        -P ${task.cpus}

    # 3. Calculate global weighted and unweighted FST statistics
    realSFS fst stats ${region}_hist_vs_mod.fst.idx > ${region}_hist_vs_mod.global_fst.txt

    # 4. Calculate sliding window FST and attach expected header
    echo -e "region\tchr\tmidPos\tNsites\tFst" > ${region}_windowed_fst.txt
    realSFS fst stats2 ${region}_hist_vs_mod.fst.idx \
        -win ${params.fst_window} \
        -step ${params.fst_step} | tail -n +2 >> ${region}_windowed_fst.txt

    # 5. Generate Manhattan plot
    Rscript ${plot_script} \
        ${region}_windowed_fst.txt \
        ${region}_fst_manhattan.png
    """
}