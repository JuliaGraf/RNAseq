process TRIMGALORE {
    debug true
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*{3prime,5prime,trimmed,val}*.fq.gz"), emit: reads
    
    script:
    """
    trimgalore ${reads[0]} ${reads[1]}
    """



}
