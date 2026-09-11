process RUN_ACER {
    tag "ACER Selection Scan"
    publishDir "${params.outdir}/selection", mode: 'copy'

    input:
    path mafs           // All ${region}_${era}.mafs.gz files from ANGSD_GL_POP
    path var_sites      // Polymorphic sites index (sites.snps)
    path helpers        // acer_helpers.R
    path run_script     // run_acer.R
    path samplesheet    // To map historic/modern pairs per region

    output:
    path "iteration_summary.tsv"     , emit: summary
    path "final_test_results.tsv"    , emit: chisq_results
    path "ne_bootstrap.tsv"          , emit: ne_boot
    path "*.png"                     , optional: true, emit: plots

    script:
    """
    # 1. Intersect population MAF files with polymorphic SNP list to discard monomorphic sites
    for maf in *.mafs.gz; do
        zcat \$maf | awk 'NR==FNR {a[\$1"\t"\$2]; next} FNR==1 || (\$1"\t"\$2) in a' ${var_sites} - | gzip > var_\${maf}
    done

    # 2. Parse sample CSV to construct --hist_mafs, --mod_mafs, and --region_names CLI arguments
    # Expects columns: sample, region, era, region
    HIST_MAFS=\$(Rscript -e '
        df <- read.csv("${samplesheet}")
        regs <- unique(df\$region)
        hist_files <- sapply(regs, function(r) {
            era <- df\$era[df\$region == r & tolower(df\$era) %in% c("historic", "historical", "hist", "era1")][1]
            paste0("var_", r, "_", era, ".mafs.gz")
        })
        cat(paste(hist_files, collapse=","))
    ')

    MOD_MAFS=\$(Rscript -e '
        df <- read.csv("${samplesheet}")
        regs <- unique(df\$region)
        mod_files <- sapply(regs, function(r) {
            era <- df\$era[df\$region == r & tolower(df\$era) %in% c("modern", "mod", "era2")][1]
            paste0("var_", r, "_", era, ".mafs.gz")
        })
        cat(paste(mod_files, collapse=","))
    ')

    REGIONS=\$(Rscript -e '
        df <- read.csv("${samplesheet}")
        cat(paste(unique(df\$region), collapse=","))
    ')

    # 3. Execute ACER iterative selection scan
    Rscript ${run_script} \\
        --hist_mafs="\$HIST_MAFS" \\
        --mod_mafs="\$MOD_MAFS" \\
        --region_names="\$REGIONS" \\
        --out_dir="." \\
        --helpers="${helpers}" \\
        --generations=${params.generations} \\
        --fdr_cutoff=${params.fdr_cutoff} \\
        --max_rounds=${params.max_rounds} \\
        --n_boot=${params.n_boot} \\
        --min_ind=${params.min_ind}
    """
}

process FILTER_NEUTRAL_SITES {
    tag "Extract Neutral Sites"
    publishDir "${params.outdir}/sites", mode: 'copy'

    input:
    path sites_snps          // Full polymorphic sites.snps file
    path cmh_results         // final_cmh_results_full.tsv from ACER

    output:
    path "neutral_sites.snps", emit: snps
    path "neutral_sites.bin" , emit: bin
    path "neutral_sites.idx" , emit: idx

    script:
    """
    # 1. Extract CHR and BP of candidate loci under selection (FDR <= cutoff)
    awk -v fdr="${params.fdr_cutoff}" 'NR>1 && \$NF <= fdr {print \$1"\t"\$2}' ${cmh_results} > selected_coords.txt

    # 2. Exclude selected loci from full sites file
    awk 'NR==FNR {selected[\$1"\t"\$2]; next} !((\$1"\t"\$2) in selected)' selected_coords.txt ${sites_snps} > neutral_sites.snps

    # 3. Re-index neutral sites for ANGSD
    angsd sites index neutral_sites.snps
    """
}