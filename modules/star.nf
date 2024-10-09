process STAR {
    publishDir 'results/star_results', mode: 'copy', pattern: "*.{bam,tab,out,gz,junction,wig,bg}"
    debug true
    
    input:
    tuple val(meta), path(reads) 
    path(genomeFasta)
    path(gtfFile)
    
    output:
    tuple val(meta), path ("*.*")
    path  "versions.yml", emit: versions

    script:
    """
    mkdir star_genome
    STAR --runMode genomeGenerate --genomeDir star_genome/ --genomeFastaFiles $genomeFasta --sjdbGTFfile $gtfFile --runThreadN 4
    STAR --readFilesIn ${reads[0]} ${reads[1]} --outFileNamePrefix star_results --genomeDir star_genome/
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
    

}