process FASTQC {
    debug true
    input:
    tuple val(meta), path(reads)

    output:
    path "*_fastqc.zip"

    publishDir 'fastqc_results', mode: 'copy', pattern: '*_fastqc.zip'

    script:
    """
    fastqc ${reads[0]} ${reads[1]}
    ls -lahtr
    """
}