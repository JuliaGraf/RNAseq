process STAR_ALIGN {
    publishDir 'results/star', mode: 'copy', pattern: "[!_]*"
    debug true
    
    input:
    tuple val(meta), path(reads)
    path gtfFile
    path index
    
    output:
    tuple val(meta), path ("*")
    tuple val(meta), path("*.bam") , emit: bam
    path  "_versions.yml"          , emit: versions

    script:
    if (reads[1] == null){
        reads = [reads[0],'']
    }
    if (reads[0] == null){
        reads = ['',reads[1]]
    }
    """
    STAR --genomeDir $index --readFilesCommand gunzip -c --readFilesIn ${reads[0]} ${reads[1]}  --outSAMtype BAM SortedByCoordinate --outTmpDir /tmp/test
    
    cat <<-END_VERSIONS > _versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
    

}