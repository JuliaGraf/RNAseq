process TRIMGALORE {
    publishDir 'results/trimgalore_results', mode: 'copy', pattern: "*{3prime,5prime,trimmed,val,report}*.{fq,gz,txt}"
    debug true

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path ("*{trimmed}*.{fq,gz}"), emit: trimmed
    path  "versions.yml", emit: versions

    script: 
    if (reads[1] == null){
        reads = [reads[0],'']
    }
    if (reads[0] == null){
        reads = ['',reads[1]]
    }
    """
    trim_galore ${reads[0]} ${reads[1]}
    
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.version //; s/Last.\$//')
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """


}
