process STAR_GENOMEGENERATE {
    publishDir 'results/', mode: 'copy', pattern: "[!_]*"
    debug true
    
    input:
    path(genomeFasta)
    path(gtfFile)
    
    output:
    path "star_genome"          , emit: index
    path "star_genome/SAindex*" , emit: fai
    path "_versions.yml"        , emit: versions

    script:
    """
    mkdir star_genome
    STAR --runMode genomeGenerate --genomeDir star_genome/  --genomeFastaFiles $genomeFasta --sjdbGTFfile $gtfFile 
    
    cat <<-END_VERSIONS > _versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
    

}