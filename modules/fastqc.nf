process FASTQC {

    input:
    tuple val(meta), val(reads)

    output:
    path "*_fastqc.zip"

    script:
    """
    fastqc ${reads[0]} ${reads[0]}
    ls -lahtr
    """
}