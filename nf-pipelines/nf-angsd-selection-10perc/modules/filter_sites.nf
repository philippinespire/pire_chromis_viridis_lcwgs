process PREPARE_SITE_SETS {
    tag "Generate Complete Site Set Library"
    publishDir "${params.outdir}/sites", mode: 'copy', pattern: "site_counts.tsv"
    publishDir "${params.outdir}/large_data/sites", mode: 'copy', pattern: "*.pos*"

    input:
    path pop_intersect  // pop_intersect.pos (Venn diagram output)
    path poly_snps     // sites.snps from ANGSD_GL_ALL
    path pruned_pos    // merged_sites from LD_PRUNE (or [])
    path chisq_results // ACER output (or [])

    output:
    path "all_callable.pos"             , emit: all_callable
    path "all_snps.pos"                 , emit: all_snps
    path "selected_loci.pos"            , emit: selected
    path "callable_neutral.pos"         , emit: callable_neutral_pos
    path "callable_neutral.pos.bin"     , emit: callable_neutral_bin
    path "callable_neutral.pos.idx"     , emit: callable_neutral_idx
    path "snps_neutral.pos"             , emit: snps_neutral
    path "snps_pruned.pos"              , emit: snps_pruned
    path "snps_pruned_neutral.pos"      , emit: snps_pruned_neutral_pos
    path "snps_pruned_neutral.pos.bin"  , emit: snps_pruned_neutral_bin
    path "snps_pruned_neutral.pos.idx"  , emit: snps_pruned_neutral_idx
    path "site_counts.tsv"              , emit: summary

    script:
    """
    cp ${pop_intersect} all_callable.pos
    awk 'NR==FNR {pop[\$1"\t"\$2]; next} (\$1"\t"\$2) in pop' all_callable.pos ${poly_snps} > all_snps.pos

    if [ -f "${pruned_pos}" ] && [ -s "${pruned_pos}" ]; then
        awk 'NR==FNR {pop[\$1"\t"\$2]; next} (\$1"\t"\$2) in pop' all_callable.pos ${pruned_pos} > snps_pruned.pos
    else
        cp all_snps.pos snps_pruned.pos
    fi

    if [ -f "${chisq_results}" ] && [ -s "${chisq_results}" ]; then
        awk -v fdr="${params.fdr_cutoff}" 'NR>1 && \$NF <= fdr {print \$1"\t"\$2}' ${chisq_results} > selected_loci.pos
    else
        touch selected_loci.pos
    fi

    if [ -s selected_loci.pos ]; then
        awk 'NR==FNR {sel[\$1"\t"\$2]; next} !((\$1"\t"\$2) in sel)' selected_loci.pos all_callable.pos > callable_neutral.pos
        awk 'NR==FNR {sel[\$1"\t"\$2]; next} !((\$1"\t"\$2) in sel)' selected_loci.pos all_snps.pos     > snps_neutral.pos
        awk 'NR==FNR {sel[\$1"\t"\$2]; next} !((\$1"\t"\$2) in sel)' selected_loci.pos snps_pruned.pos   > snps_pruned_neutral.pos
    else
        cp all_callable.pos callable_neutral.pos
        cp all_snps.pos     snps_neutral.pos
        cp snps_pruned.pos  snps_pruned_neutral.pos
    fi

    # Index position lists for realSFS site filtering
    angsd sites index callable_neutral.pos
    angsd sites index snps_pruned_neutral.pos

    # Generate site summary using redirected wc -l (avoids filename string parsing)
    echo -e "site_set\tcount" > site_counts.tsv
    echo -e "all_callable\t\$(wc -l < all_callable.pos)" >> site_counts.tsv
    echo -e "all_snps\t\$(wc -l < all_snps.pos)" >> site_counts.tsv
    echo -e "selected_loci\t\$(wc -l < selected_loci.pos)" >> site_counts.tsv
    echo -e "callable_neutral\t\$(wc -l < callable_neutral.pos)" >> site_counts.tsv
    echo -e "snps_neutral\t\$(wc -l < snps_neutral.pos)" >> site_counts.tsv
    echo -e "snps_pruned\t\$(wc -l < snps_pruned.pos)" >> site_counts.tsv
    echo -e "snps_pruned_neutral\t\$(wc -l < snps_pruned_neutral.pos)" >> site_counts.tsv
    """
}