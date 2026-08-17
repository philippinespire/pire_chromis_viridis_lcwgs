nextflow.enable.dsl=2

// Default Parameters
params.samplesheet = "${projectDir}/inputfiles/samplesheet.csv" // sample sheet csv with columns for sample name, bam path, pop, era, and region. era is historic or modern. region matches historic and modern populations for diversity calculations. pop is not used
params.contigs     = "${projectDir}/inputfiles/contig_list.txt" // text file with one contig name per line. Contig names must match those in the reference genome FASTA file.
params.outdir      = "${projectDir}/results" // output directory for all results
params.reference   = "/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/data/GCA_051013605.1_ASM5101360v1_genomic_20kb.fna" // path to reference genome FASTA file
params.bed_file    = "/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-trim-generode/results-iridian/data/reference/GCA_051013605.1_ASM5101360v1_genomic_20kb.repma.angsd.txt" // a 3-column ANGSD-format bed file. ideally repeat-masked.
params.species     = "Cvi"   // 3-letter species code for output file naming. used in ANGSD output files and ACER results.
params.min_ind_frac   = 0.70 // Minimum percentage of individuals required per site (default 0.70)
params.max_depth_mult = 10   // Max depth threshold as a multiplier of the average total expected coverage (default 10)
params.ld_prune    = true    // Set to true to enable LD pruning (default true)
params.max_kb_dist = 50      // Maximum pairwise distance in kb to test for LD if pruning
params.min_weight  = 0.2     // Minimum r2 threshold for pruning filter (default 0.2)
params.run_selection = true  // Set to false to skip ACER by default (default true)
params.generations = 114     // Number of generations to use for Ne calculations. Default is 114 for Albatross to contemporary populations with a 1 year generation time. Adjust for other species as needed.
params.fdr_cutoff  = 0.05    // FDR cutoff for iterative selection scan (default 0.05)
params.max_rounds  = 20      // Maximum number of rounds to run for iterative selection scan (default 20)
params.n_boot      = 1000    // Number of bootstrap replicates to run for Ne calculations (default 1000)
params.min_ind     = 4       // Minimum number of individuals to include in Ne calculations (default 4)
params.fst_window = 50000    // Window size in bp for sliding window FST calculations (default 50000)
params.fst_step   = 10000    // Step size in bp for sliding window FST calculations (default 10000)

// import modules
include { ANGSD_GL_ALL; ANGSD_COLLECT_OUTPUT; ANGSD_EXTRACT_SITES; INDEX_BED_SITES  } from './modules/angsd_gl.nf'
include { COLLECT_BAM_ALL; COLLECT_BAM_POP                                          } from './modules/collect_bam.nf'
include { PCANGSD; PLOT_PCANGSD; PLOT_ADMIXTURE                                     } from './modules/pcangsd.nf'
include { CALC_POP_THRESHOLDS; QC_POP_SITES_CONTIG; GATHER_POP_SITES                } from './modules/angsd_gl_pop.nf'
include { INTERSECT_POP_SITES; ANGSD_GL_POP                                         } from './modules/angsd_gl_pop.nf'
include { ANGSD_DIVERSITY; PLOT_PI_DIVERSITY                                        } from './modules/angsd_diversity.nf'
include { LD_PRUNE_CONTIG; MERGE_PRUNED_SITES; SUBSET_BEAGLE                        } from './modules/ld_prune.nf'
include { PREPARE_SITE_SETS                                                         } from './modules/filter_sites.nf'
include { RUN_ACER                                                                  } from './modules/acer.nf'
include { ANGSD_FST                                                                 } from './modules/angsd_fst.nf'

// Workflow 
workflow {
    // --- Input Channels (Defined inside workflow for DSL2 scope compliance) ---
    ch_ref = Channel.value([ file(params.reference), file("${params.reference}.fai") ])
    
    samples = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, row.pop, row.era, row.region, row.bam) }

    all_bams = samples
        .map { sample, pop, era, region, bam -> bam }
        .collect()

    contigs = Channel
        .fromPath(params.contigs, checkIfExists: true)
        .splitText()
        .map { it.trim() }
        .filter { it }

    // Group samples by region (index 3) and era (index 2)
    pop_bams = samples
        .groupTuple(by: [3, 2])
        .map { sample_list, pop, era, region, bams -> tuple(region, era, bams) }

    bamlist = COLLECT_BAM_ALL(all_bams)
    bamlist_pop = COLLECT_BAM_POP(pop_bams)

    // 1. Index full callable bed file for full-genome diversity calculations (no SNP_pval filter)
    all_sites = INDEX_BED_SITES(params.bed_file)
    ch_sites_bundle = all_sites.snps
        .combine(all_sites.bin)
        .combine(all_sites.idx)

    // 2. Calculate mindind and maxdepththresholds once per population/era
    CALC_POP_THRESHOLDS(bamlist_pop)

    // 3. Run per-population QC to identify loci passing minind_frac and maxdepth_mult in historical and modern subsets
    // Combine the calculated thresholds with the existing contigs channel
    ch_qc_inputs = CALC_POP_THRESHOLDS.out.thresholds.combine(contigs)

    // parallelized by contig and population/era 
    QC_POP_SITES_CONTIG(
        ch_qc_inputs,
        ch_ref,
        ch_sites_bundle.first()
    )
    
    // 4. Group by population and era keys [pop, era], then gather into single .pos file
    ch_grouped_pos = QC_POP_SITES_CONTIG.out.contig_pos
        .groupTuple(by: [0, 1])

    GATHER_POP_SITES(ch_grouped_pos)

    // 5. Intersect site lists across historical and modern populations
    ch_all_passing_pos = GATHER_POP_SITES.out.passing_pos.collect()
    INTERSECT_POP_SITES(ch_all_passing_pos)
    ch_intersect_bundle = INTERSECT_POP_SITES.out.pos
        .combine(INTERSECT_POP_SITES.out.bin)
        .combine(INTERSECT_POP_SITES.out.idx)

    // 6. Generate SAF and MAF files for historical and modern populations
    // using ONLY the shared intersected site set
    ANGSD_GL_POP(
        bamlist_pop,
        ch_ref,
        ch_intersect_bundle.first()
    )
    
    // 7. Run ANGSD_GL_ALL to identify SNPs (uses SNP_pval 1e-6)
    genotypes = ANGSD_GL_ALL(contigs, bamlist.first(), ch_sites_bundle.first(), ch_ref)

    mafs = genotypes.mafs.map { contig, file -> file }.collect()
    beagle = genotypes.beagle.map { contig, file -> file }.collect()

    collected = ANGSD_COLLECT_OUTPUT(mafs, beagle)
    sites = ANGSD_EXTRACT_SITES(collected.all_mafs)

    // 8. LD Pruning (optional)
    if (params.ld_prune) {
        contig_ld_input     = genotypes.beagle.join(genotypes.mafs)
        pruned_contig_sites = LD_PRUNE_CONTIG(contig_ld_input)
        ch_pruned_pos       = MERGE_PRUNED_SITES(pruned_contig_sites.collect())
    } else {
        ch_pruned_pos       = Channel.value([])
    }

    // 9. ACER Chi-Squared selection scan (optional)
    if (params.run_selection) {
        RUN_ACER(
            ANGSD_GL_POP.out.mafs.collect(),
            sites.snps,
            file("${projectDir}/scripts/acer_helpers.R"),
            file("${projectDir}/scripts/run_acer.R"),
            file(params.samplesheet)
        )
        ch_chisq_results = RUN_ACER.out.chisq_results
    } else {
        ch_chisq_results = Channel.value([])
    }

     // 10. Generate Site Set Libraries for downstream uses
    PREPARE_SITE_SETS(
        INTERSECT_POP_SITES.out.pos,
        sites.snps,
        ch_pruned_pos,
        ch_chisq_results
    )

    // 11. PCA & Admixture (Runs on LD-pruned neutral SNPs)
    subset_beagle   = SUBSET_BEAGLE(collected.all_beagle, PREPARE_SITE_SETS.out.snps_pruned_neutral_pos)
    pcangsd_results = PCANGSD(subset_beagle.pruned_beagle)
    PLOT_PCANGSD(pcangsd_results.pcangsd_cov, file(params.samplesheet))
    PLOT_ADMIXTURE(pcangsd_results.q_matrix, file(params.samplesheet))

    // 12. Diversity Calculations (Runs on all callable neutral loci)
    ch_diversity_bundle = PREPARE_SITE_SETS.out.callable_neutral_pos
        .combine(PREPARE_SITE_SETS.out.callable_neutral_bin)
        .combine(PREPARE_SITE_SETS.out.callable_neutral_idx)
    ANGSD_DIVERSITY(ANGSD_GL_POP.out.saf_files, ch_diversity_bundle.first())

    // Split diversity outputs into historic and modern channels
    ch_pest_hist = ANGSD_DIVERSITY.out.pestPG
        .filter { region, era, pestPG -> era.toLowerCase() in ['historic', 'historical', 'hist', 'era1'] }
        .map { region, era, pestPG -> tuple(region, pestPG) }

    ch_pest_mod = ANGSD_DIVERSITY.out.pestPG
        .filter { region, era, pestPG -> era.toLowerCase() in ['modern', 'mod', 'era2'] }
        .map { region, era, pestPG -> tuple(region, pestPG) }

    // Join channels on region key -> outputs [ region, hist_pestPG, mod_pestPG ]
    ch_paired_thetas = ch_pest_hist.join(ch_pest_mod)

    // Run diversity plotting per region
    PLOT_PI_DIVERSITY(
        ch_paired_thetas,
        file("${projectDir}/scripts/plot_tp.R")
    )

    // 13. Calc historic-modern FST
    // Split SAF channel into historical and modern channels
    ch_saf_hist = ANGSD_GL_POP.out.saf_files
        .filter { region, era, saf, idx, pos -> era.toLowerCase() in ['historic', 'historical', 'hist', 'era1'] }
        .map { region, era, saf, idx, pos -> tuple(region, saf, idx, pos) }

    ch_saf_mod = ANGSD_GL_POP.out.saf_files
        .filter { region, era, saf, idx, pos -> era.toLowerCase() in ['modern', 'mod', 'era2'] }
        .map { region, era, saf, idx, pos -> tuple(region, saf, idx, pos) }

    // Join on region key (region)
    ch_fst_pairs = ch_saf_hist.join(ch_saf_mod)

    // Bundle pruned neutral site files
    ch_pruned_neutral_bundle = PREPARE_SITE_SETS.out.snps_pruned_neutral_pos
        .combine(PREPARE_SITE_SETS.out.snps_pruned_neutral_bin)
        .combine(PREPARE_SITE_SETS.out.snps_pruned_neutral_idx)

    // Execute FST calculation
    ANGSD_FST(ch_fst_pairs, ch_pruned_neutral_bundle.first(), file("${projectDir}/scripts/plot_windowed_fst.R"))
}