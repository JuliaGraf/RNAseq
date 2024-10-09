process MARKDUPLICATES {
    publishDir 'results/markduplicates_results', mode: 'copy', pattern: "*.*"
    debug true
    input:
    tuple val(meta), path(bam)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.bam")        , emit: bam
    tuple val(meta), path("*.bai")        , optional:true, emit: bai
    tuple val(meta), path("*.metrics.txt"), emit: metrics
    path  "versions.yml"                  , emit: versions

    script:
    def prefix = "${meta.sample}"

    """
    picard MarkDuplicates \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.bam \\
        --REFERENCE_SEQUENCE $fasta \\
        --METRICS_FILE ${prefix}.MarkDuplicates.metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(echo \$(picard MarkDuplicates --version 2>&1) | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}