process SALMON_QUANTIFICATION {
    publishDir 'results/', mode: 'copy', pattern: '[!_]*'

    input:
    tuple val(meta), path(reads)
    path  index
    path  gtf

    output:
    path "salmon_quantification"                    
    path  "_versions.yml"            , emit: versions

    script:

    //Check if reads are single-end or paired-end
    if (reads[1] == null){
        reads = '-r '+reads[0] 
    }
    else if (reads[0] == null){
        reads = '-r '+reads[1]
    }
    else {
        reads = '-1 '+reads[0]+' -2 '+reads[1]
    }
    """
    salmon quant --geneMap ${gtf} --index ${index} -l A ${reads} -o salmon_quantification

    cat <<-END_VERSIONS > _versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}