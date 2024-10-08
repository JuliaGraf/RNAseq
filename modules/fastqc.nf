process FASTQC {

    input:
    tuple val(meta), val(reads)

    script:
    """
    fastqc ${reads[0]} ${reads[0]}
    ls -lahtr
    """
}