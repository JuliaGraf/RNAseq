process  SALMON_INDEX {
    publishDir 'results/', mode: 'copy', pattern: '[!_]*'

    input:
    path(transcripts)

    output:
    path "salmon_index",          emit: index
    path "_versions.yml", emit: versions

    script:
    """
    salmon index -t ${transcripts} -i salmon_index
    
    cat <<-END_VERSIONS > _versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """


}