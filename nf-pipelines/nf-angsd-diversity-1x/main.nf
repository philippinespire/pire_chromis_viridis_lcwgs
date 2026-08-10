nextflow.enable.dsl=2

// Default Parameters
params.samplesheet = "${projectDir}/inputfiles/samplesheet.csv"
params.contigs     = "${projectDir}/inputfiles/contig_list.txt"
params.outdir      = "${projectDir}/results"
params.reference   = "/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/data/GCA_051013605.1_ASM5101360v1_genomic_20kb.fna"
params.bed_file    = "/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-generode/results-iridian/data/reference/GCA_051013605.1_ASM5101360v1_genomic_20kb.repma.angsd.txt"
params.species     = "Cvi"
params.maxdepth    = 461 // maximum depth to include in analysis. 10x the expected depth across the 26 individuals.
params.minind      = 18 // minimum number of individuals to include in analysis. 70 percent of 26 individuals.
params.ld_prune    = true  // Set to true to enable LD pruning
params.max_kb_dist = 50     // Maximum pairwise distance in kb to test for LD if pruning
params.min_weight  = 0.2    // Minimum r2 threshold for pruning filter

// import modules
include { ANGSD_GL_ALL; ANGSD_COLLECT_OUTPUT; ANGSD_EXTRACT_SITES; INDEX_BED_SITES} from './modules/angsd_gl'
include { COLLECT_BAM_ALL; COLLECT_BAM_POP } from './modules/collect_bam'
include { PCANGSD; PLOT_PCANGSD; PLOT_ADMIXTURE } from './modules/pcangsd'
include { ANGSD_GL_POP } from './modules/angsd_gl_pop'
include { ANGSD_DIVERSITY } from './modules/angsd_diversity'

// Processes
process LD_PRUNE_CONTIG {
    tag "LD Pruning ${contig}"
    // Note: publishDir is removed here so it doesn't flood your results folder with 100s of tiny files

    input:
    tuple val(contig), path(beagle), path(maf)

    output:
    path "${contig}_pruned_sites.pos"

    script:
    """
    module load container_env ngsTools

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
    crun ngsLD \\
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
    
    first_dist_field=\$(awk 'NR==1 {print \$7; exit}' ${contig}.ld)
    if [[ "\$first_dist_field" =~ ^[0-9]+([.][0-9]+)?\$ ]]; then
        awk 'BEGIN{OFS="\t"; print "site1","site2","dist","r2"} \$8 != "" {print \$1,\$4,\$7,\$8}' ${contig}.ld > ${contig}_prune_input.tsv
    else
        awk 'BEGIN{OFS="\t"; print "site1","site2","dist","r2"} NR>1 && \$8 != "" {print \$1,\$4,\$7,\$8}' ${contig}.ld > ${contig}_prune_input.tsv
    fi

    # 7. Safety check: If the TSV only has a header (<= 1 line), skip prune_graph
    if [ \$(wc -l < ${contig}_prune_input.tsv) -le 1 ]; then
        cp ${contig}.pos ${contig}_pruned_sites.pos
        exit 0
    fi

    weight_filter="dist <= \${max_bp_dist} && r2 >= ${params.min_weight}"

    # 8. Run prune_graph on this tiny graph
    crun prune_graph \\
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
    publishDir "${params.outdir}/ld_pruning", mode: 'copy'

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
    publishDir "${params.outdir}/ld_pruning", mode: 'copy'

    input:
    path original_beagle
    path pruned_pos

    output:
    path "pruned.beagle.gz", emit: pruned_beagle

    script:
    """
    # Match the marker ID column (col 1 of beagle) to the pruned positions.
    # This robust parser handles chromosomes containing any number of underscores.
    zcat ${original_beagle} | awk -F'\\t' -v OFS='\\t' '
        FNR == NR { 
            # Robustly split the pruned_pos columns by any whitespace (tabs or spaces)
            split(\$0, a, /[ \\t]+/)
            if (a[1] != "" && a[2] != "") {
                keys[a[1], a[2]] = 1 
            }
            next 
        }
        {
            if (FNR == 1) { 
                print \$0 
                next 
            }
            
            # Extract chromosome and position from marker ID (col 1)
            marker = \$1
            # Clean up trailing alleles if they exist (e.g., _A_G or _C_T)
            sub(/(_[A-Za-z])+\$/, "", marker)
            
            # Find the last underscore which separates chromosome from position
            last_under = 0
            for (i = length(marker); i > 0; i--) {
                if (substr(marker, i, 1) == "_") {
                    last_under = i
                    break
                }
            }
            
            if (last_under > 0) {
                chr = substr(marker, 1, last_under - 1)
                pos = substr(marker, last_under + 1)
                if ((chr, pos) in keys) {
                    print \$0
                }
            }
        }
    ' ${pruned_pos} - | gzip -c > pruned.beagle.gz
    """
}

// Workflow 

workflow {
    // --- Input Channels (Defined inside workflow for DSL2 scope compliance) ---
    samples = channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.pop, row.era, row.region, row.bam) }

    all_bams = samples
        .map { sample, pop, era, region, bam -> bam }
        .collect()

    contigs = channel
        .fromPath(params.contigs, checkIfExists: true)
        .splitText()
        .map { it.trim() }
        .filter { it }

    pop_bams = samples
        .groupTuple(by: [1, 2])
        .map { sample_list, pop, era, regions, bams -> tuple(pop, era, bams) }

    samplesheet_file = channel.fromPath(params.samplesheet)

    bamlist = COLLECT_BAM_ALL(all_bams)
    bamlist_pop = COLLECT_BAM_POP(pop_bams)

    // 1. Index full callable bed file for UNBIASED diversity calculations (No SNP_pval filter)
    all_sites = INDEX_BED_SITES(params.bed_file)

    // 2. Run ANGSD_GL_ALL for PCA (uses SNP_pval 1e-6)
    genotypes = ANGSD_GL_ALL(contigs, bamlist)

    mafs = genotypes.mafs.map { contig, file -> file }.collect()
    beagle = genotypes.beagle.map { contig, file -> file }.collect()

    collected = ANGSD_COLLECT_OUTPUT(mafs, beagle)
    sites = ANGSD_EXTRACT_SITES(collected.all_mafs)

    // Setup channel paths for PCA (filtered & pruned) and Diversity (unpruned)
    def final_beagle = collected.all_beagle
    def pruned_files = []

    if (params.ld_prune) {
        // Match the per-contig beagle and maf files together
        contig_ld_input = genotypes.beagle.join(genotypes.mafs)
        
        // Run ngsLD and prune_graph on a per-contig basis
        pruned_contig_sites = LD_PRUNE_CONTIG(contig_ld_input)
        
        // Merge the individual contig pos files into one complete list
        merged_sites = MERGE_PRUNED_SITES(pruned_contig_sites.collect())
                
        // Index pruned sites for downstream ANGSD runs
        indexed = INDEX_PRUNED_SITES(merged_sites)
        
        // Subset the Beagle file using the pruned positions
        subset_beagle = SUBSET_BEAGLE(collected.all_beagle, merged_sites)
        
        final_beagle = subset_beagle.pruned_beagle
        pruned_files = indexed.snps.combine(indexed.bin).combine(indexed.idx).first()
    }

    regions = sites.regions

    // 3. Run PCAngsd on filtered (+ pruned) beagle file
    pcangsd_results = PCANGSD(final_beagle)
    PLOT_PCANGSD(pcangsd_results.pcangsd_cov, file(params.samplesheet))
    PLOT_ADMIXTURE(pcangsd_results.q_matrix, file(params.samplesheet))
    
    // 4. Run ANGSD_GL_POP on ALL callable sites (No SNP_pval filter, No LD pruning)
    ANGSD_GL_POP(bamlist_pop, all_sites.snps, all_sites.bin, all_sites.idx, regions)

    // 5. Diversity Calculations: Pi & Theta calculated on full SFS; 
    // Outputs both ld-pruned and un-pruned SFS versions for downstream uses.
    ANGSD_DIVERSITY(ANGSD_GL_POP.out.saf_files, pruned_files)
}