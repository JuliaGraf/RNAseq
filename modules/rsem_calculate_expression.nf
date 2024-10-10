process RSEM_CALCULATE_EXPRESSION{
    publishDir 'results/rsem_calculate_exp', mode: 'copy', pattern: '*'

    input:
    tuple val(meta), path(reads)
    path(in_dir)

    output:
    tuple val(meta), path("*.genes.results")   , emit: counts_gene
    tuple val(meta), path("*.isoforms.results"), emit: counts_transcript
    tuple val(meta), path("*.stat")            , emit: stat
    tuple val(meta), path("*.log")             , emit: logs
    path  "versions.yml"                       , emit: versions

    script:
    def prefix = "${meta.sample}"
    if (reads[1] == null && reads[0] != null){  
        """
        rsem-calculate-expression ${reads[0]} rsem_out ${prefix} --temporary-folder ./tmp
        """
    } else if (reads[0] == null && reads[1] != null){
        """
        rsem-calculate-expression ${reads[1]} rsem_out ${prefix} --temporary-folder ./tmp
        """
    } else if (reads[0] != null && reads[1] != null) {
        """
        rsem-calculate-expression --paired-end ${reads[0]} ${reads[1]} rsem_out ${prefix} --temporary-folder ./tmp
        """
    }
}