nextflow.enable.dsl=2

// Default Parameters
params.samplesheet = "${projectDir}/inputfiles/samplesheet.csv"
params.contigs     = "${projectDir}/data/reference/reference.contigs.txt"
params.outdir      = "${projectDir}/results"
params.reference   = "${projectDir}/data/reference/reference.fasta"
params.bed_file    = "${projectDir}/data/reference/reference.repma.angsd.txt" // create awk '{print $1"\t"$2+1"\t"$3}' reference.repma.bed > reference.repma.angsd.txt
params.species     = "Sor"
params.maxdepth    = 1000
params.minind      = 20

// import modules
include { ANGSD_GL_ALL; ANGSD_COLLECT_OUTPUT; ANGSD_EXTRACT_SITES } from './modules/angsd_gl'
include { COLLECT_BAM_ALL; COLLECT_BAM_POP } from './modules/collect_bam'
include { PCANGSD; PLOT_PCANGSD } from './modules/pcangsd'
include { ANGSD_GL_POP } from './modules/angsd_gl_pop'
include { ANGSD_DIVERSITY } from './modules/angsd_diversity'


// --- Input Channels ---

// Read in samplesheet and collect metadata

samples = Channel
    .fromPath(params.samplesheet, checkIfExists: true)
    .splitCsv(header: true)
    .map { row -> tuple(row.sample, row.pop, row.era,row.region, row.bam) }

all_bams = samples
    .map { sample, pop, era, region, bam -> bam }
    .collect()

contigs = Channel
    .fromPath(params.contigs, checkIfExists: true)
    .splitText()
    .map { it.trim() }
    .filter { it }

pop_bams = samples
    .groupTuple(by: [1, 2])
    .map { samples, pop, era, regions, bams -> tuple(pop, era, bams)}

// Workflow 

workflow {

    samplesheet_file = Channel.fromPath(params.samplesheet)

    bamlist = COLLECT_BAM_ALL(all_bams)

    bamlist_pop = COLLECT_BAM_POP(pop_bams)

    genotypes = ANGSD_GL_ALL(contigs, bamlist)

    mafs = genotypes.mafs
        .map { contig, file -> file }
        .collect()

    beagle = genotypes.beagle
        .map { contig, file -> file }
        .collect()

    collected = ANGSD_COLLECT_OUTPUT(mafs, beagle)

    sites = ANGSD_EXTRACT_SITES(collected.all_mafs)
    snps = sites.snps
    bin = sites.snps_bin
    idx = sites.snps_idx
    regions = sites.regions

    // PCAngsd (optional)
     pcangsd_cov = PCANGSD(collected.all_beagle)
     PLOT_PCANGSD(pcangsd_cov, samplesheet_file)
    
    ANGSD_GL_POP(bamlist_pop,snps, bin,idx,regions)

    ANGSD_DIVERSITY(ANGSD_GL_POP.out.saf_files)


    

}
