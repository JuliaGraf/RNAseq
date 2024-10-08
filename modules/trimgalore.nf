process TRIMGALORE {
    publishDir 'trimgalore_results', mode: 'copy', pattern: "*{3prime,5prime,trimmed,val,report}*.{fq,gz,txt}"
    debug true

    input:
    tuple val(meta), path(reads)

    output:
    path "*.*"
    path  "versions.yml", emit: versions


    script:
    """
    trim_galore ${reads[0]} ${reads[1]}
    
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.version //; s/Last.\$//')
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """


}
