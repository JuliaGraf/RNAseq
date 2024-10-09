process FASTQC {
    publishDir 'results/fastqc_results', mode: 'copy', pattern: '*_fastqc.zip'

    input:
    tuple val(meta), path(reads)

    output:
    path "*_fastqc.zip"
    path  "versions.yml", emit: versions

    script:
    if (reads[1] == null){
        reads = [reads[0],'']
    }
    if (reads[0] == null){
        reads = ['',reads[1]]
    }
    """
    fastqc ${reads[0]} ${reads[1]}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}