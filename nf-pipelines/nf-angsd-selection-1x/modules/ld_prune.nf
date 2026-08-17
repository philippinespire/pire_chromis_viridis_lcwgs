process LD_PRUNE_CONTIG {
    tag "LD Pruning ${contig}"
    // Note: publishDir is removed here so it doesn't flood your results folder with 100s of tiny files

    input:
    tuple val(contig), path(beagle), path(maf)

    output:
    path "${contig}_pruned_sites.pos"

    script:
    """
    # 1. Create a pos file directly from the MAF file for this specific contig
    zcat ${maf} | tail -n +2 | awk 'BEGIN{OFS="\t"} {print \$1, \$2}' > ${contig}.pos

    # 2. Count fields and sites
    beagle_ncols=\$(zcat ${beagle} | head -n 1 | awk '{print NF}')
    n_ind=\$(( (beagle_ncols - 3) / 3 ))
    n_sites=\$(wc -l < ${contig}.pos)

    # 3. If contig has < 2 sites, skip ngsLD and return empty file
    if [ "\$n_sites" -lt 2 ]; then
        touch ${contig}_pruned_sites.pos
        exit 0
    fi

    # 4. Run ngsLD on just this contig
    ngsLD \\
        --geno ${beagle} \\
        --probs \\
        --pos ${contig}.pos \\
        --n_ind \$n_ind \\
        --n_sites \$n_sites \\
        --max_kb_dist ${params.max_kb_dist} \\
        --n_threads ${task.cpus ?: 1} \\
        --out ${contig}.ld

    # 5. If ngsLD found no linked edges, keep all sites
    if [ ! -s ${contig}.ld ]; then
        cp ${contig}.pos ${contig}_pruned_sites.pos
        exit 0
    fi

    # 6. Format the pairwise table (Filter out empty r2 columns)
    max_bp_dist=\$(awk -v kb="${params.max_kb_dist}" 'BEGIN{printf "%.6f", kb*1000}')

    # Extract site1 (1), site2 (2), dist (3), and r2_EM (7)
    awk 'BEGIN{OFS="\t"; print "site1","site2","dist","r2"} \$7 != "" {print \$1,\$2,\$3,\$7}' ${contig}.ld > ${contig}_prune_input.tsv

    # 7. Safety check: If the TSV only has a header (<= 1 line), skip prune_graph
    if [ \$(wc -l < ${contig}_prune_input.tsv) -le 1 ]; then
        cp ${contig}.pos ${contig}_pruned_sites.pos
        exit 0
    fi

    weight_filter="dist <= \${max_bp_dist} && r2 >= ${params.min_weight}"

    # 8. Run prune_graph on this tiny graph
    prune_graph \\
        --header \\
        --in ${contig}_prune_input.tsv \\
        --weight-field "r2" \\
        --weight-filter "\$weight_filter" \\
        --out ${contig}_pruned.raw.pos

    # 9. Reconstruct the 2-column pos file mapping
    awk '
      FNR == NR {
        chr = \$1; pos = \$2;
        lines[chr "_" pos] = chr "\t" pos;
        lines[chr ":" pos] = chr "\t" pos;
        lines[pos] = chr "\t" pos;
        next
      }
      { if (\$1 in lines) print lines[\$1] }
    ' ${contig}.pos ${contig}_pruned.raw.pos > ${contig}_pruned_sites.pos
    """
}

process MERGE_PRUNED_SITES {
    tag "Merge Pruned Sites"
    //publishDir "${params.outdir}/ld_pruning", mode: 'copy'

    input:
    path pos_files

    output:
    path "pruned_sites.pos"

    script:
    """
    cat ${pos_files} | sort -k1,1 -k2,2n > pruned_sites.pos
    """
}

process INDEX_PRUNED_SITES {
    tag "Index Pruned Sites"
    publishDir "${params.outdir}/ld_pruning", mode: 'copy'

    input:
    path pruned_sites

    output:
    path "pruned_sites.pos", emit: snps
    path "pruned_sites.pos.bin", emit: bin
    path "pruned_sites.pos.idx", emit: idx

    script:
    """
    # Index the pruned sites with ANGSD so they can be read downstream
    angsd sites index ${pruned_sites}
    """
}

process SUBSET_BEAGLE {
    tag "Subset Beagle File"
    // publishDir "${params.outdir}/large_data/beagle", mode: 'copy' // publishDir removed to save disk space. Reconstruct from sites files and full beagle if needed.

    input:
    path original_beagle
    path pruned_pos

    output:
    path "pruned.beagle.gz", emit: pruned_beagle

    script:
    """
    # Match the marker ID column (col 1 of beagle) to the pruned positions.
    # This robust parser handles chromosomes containing any number of underscores.
    zcat ${original_beagle} | awk -F'\t' -v OFS='\t' '
        FILENAME == ARGV[1] {
            split(\$0, a, /[ \t]+/)
            if (a[1] != "" && a[2] != "") keys[a[1], a[2]] = 1
            next
        }
        {
            if (FNR == 1) { print \$0; next }
            marker = \$1
            sub(/(_[A-Za-z])+\$/, "", marker)
            last_under = 0
            for (i = length(marker); i > 0; i--) {
                if (substr(marker, i, 1) == "_") { last_under = i; break }
            }
            if (last_under > 0) {
                chr = substr(marker, 1, last_under - 1)
                pos = substr(marker, last_under + 1)
                if ((chr, pos) in keys) print \$0
            }
        }
    ' ${pruned_pos} - | gzip -c > pruned.beagle.gz
    """
}