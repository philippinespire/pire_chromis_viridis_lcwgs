process PCANGSD {
    // calculates pca covariance matrix and admixture proportions
    publishDir "${params.outdir}/PCAngsd", mode: 'copy'

    input:
    path (all_beagle)

    output:
    path "${params.species}.cov", emit: pcangsd_cov
	path "${params.species}.*.Q", emit: q_matrix

    script:
    """
    pcangsd -b ${all_beagle} --admix -t ${task.cpus ?: 8} -o ${params.species}
    """
}

process PLOT_PCANGSD {

    publishDir "${params.outdir}/PCAngsd", mode: 'copy'

    input:
    path (pcangsd_cov)
    path samplesheet_file

    output:
    path "${params.species}.pcangsd.plot.pdf"

    script:
    """
    Rscript ${projectDir}/scripts/pcangsd.R \\
        --cov ${pcangsd_cov} \\
        --samplesheet ${samplesheet_file} \\
        --species ${params.species} \\
        --out ${params.species}.pcangsd.plot.pdf
    """
}

process PLOT_ADMIXTURE {
    
    publishDir "${params.outdir}/PCAngsd", mode: 'copy'

    input:
    path q_file
    path metadata

    output:
    path "${params.species}.admixture.pdf"

    script:
    """
    # Runs the provided R script using the matrix and metadata
    Rscript ${projectDir}/scripts/plot_admixture.R \\
        ${q_file} \\
        ${metadata} \\
        ${params.species}.admixture.pdf
    """
}