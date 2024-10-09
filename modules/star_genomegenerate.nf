process STAR_GENOMEGENERATE {
    publishDir 'results/star_genome', mode: 'copy', pattern: "*.*"
    debug true
    
    input:
    path(genomeFasta)
    path(gtfFile)
    
    output:
    path "star_genome" , emit: index
    path "versions.yml", emit: versions

    script:
    """
    mkdir star_genome
    STAR --runMode genomeGenerate --genomeDir star_genome/ --genomeFastaFiles $genomeFasta --genomeSAsparseD 2 --sjdbGTFfile $gtfFile
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
    

}