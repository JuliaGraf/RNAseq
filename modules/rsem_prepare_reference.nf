process RSEM_PREPARE_REFERENCE {
    publishDir 'results/rsem_prepare_reference', mode: 'copy', pattern: '*'
    debug true

    input:
    path(fasta)
    path(gtf)

    output:
    path "*rsem_out*", emit: out_dir
    path ".", emit:index
    path "versions.yml"   , emit: versions

    script:
    """
    rsem-prepare-reference --gtf $gtf $fasta rsem_out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rsem: \$(rsem-calculate-expression --version | sed -e "s/Current version: RSEM v//g")
    END_VERSIONS
    """

}