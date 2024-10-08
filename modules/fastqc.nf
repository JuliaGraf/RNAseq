process FASTQC {
    debug true
    input:
    tuple val(meta), path(reads)

    output:
    path "*_fastqc.zip"

    script:
    """
    fastqc ${reads[0]} ${reads[0]}
    ls -lahtr
    """
}