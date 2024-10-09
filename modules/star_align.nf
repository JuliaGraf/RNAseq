process STAR_ALIGN {
    publishDir 'results/star_results', mode: 'copy', pattern: "*.*"
    debug true
    
    input:
    tuple val(meta), path(reads)
    path gtfFile
    path index
    
    output:
    tuple val(meta), path ("*.*")
    path  "versions.yml", emit: versions

    script:
    if (reads[1] == null){
        reads = [reads[0],'']
    }
    if (reads[0] == null){
        reads = ['',reads[1]]
    }
    """
    STAR --genomeDir $index --readFilesCommand gunzip -c --readFilesIn ${reads[0]} ${reads[1]} --outFileNamePrefix star_results
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
    

}