process FASTQC {
    publishDir 'results/fastqc', mode: 'copy', pattern: '[!_]*'

    input:
    tuple val(meta), path(reads)

    output:
    path "*_fastqc.zip"   , emit: zip
    path "*_fastqc.html"  , emit: report
    path  "_versions.yml" , emit: versions

    script:

    //Check if reads are single-end or paired-end
    if (reads[1] == null){
        reads = [reads[0],'']
    }
    if (reads[0] == null){
        reads = ['',reads[1]]
    }
    """
    fastqc ${reads[0]} ${reads[1]}

    cat <<-END_VERSIONS > _versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}