process SALMON_QUANTIFICATION {
    publishDir 'results/', mode: 'copy', pattern: '[!_]*'
    debug true


    input:
    tuple val(meta), path(reads)
    path  index
    path  gtf


    output:
    path "salmon_quantification"                    
    path  "_versions.yml"                            , emit: versions

    script:
    if (reads[1] == null){
        reads = ['-r',reads[0],'','']
    }
    else if (reads[0] == null){
        reads = ['-r',reads[1],'','']
    }
    else {
        reads = ['-1',reads[0],'-2',reads[1]]
    }
    """
    salmon quant --geneMap ${gtf} --index ${index} -l A ${reads[0]} ${reads[1]} ${reads[2]} ${reads[3]} -o salmon_quantification

    cat <<-END_VERSIONS > _versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """


}