process PCANGSD {

    publishDir "${params.outdir}/PCAngsd", mode: 'copy'

    input:
    path (all_beagle)

    output:
    path "${params.species}.pcangsd.cov", emit: pcangsd_cov

    script:
    """
    pcangsd -b ${all_beagle} -t 8 -o ${params.species}.pcangsd
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