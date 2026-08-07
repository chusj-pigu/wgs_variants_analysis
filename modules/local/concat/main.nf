process SEQKIT_STATS {

    input:
    tuple val(meta), path(fastqFiles)

    output:
    tuple val(meta), path("stats.txt")

    script:
    """
    cat ${fastqFiles.join(' ')} | seqkit stats > stats.txt
    """
}