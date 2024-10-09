process STAR {
    publishDir 'results/star_results', mode: 'copy', pattern: "*.{bam,tab,out,gz,junction,wig,bg}"
    debug true
    
    input:
    tuple val(meta), path(reads) 
    path(genomeFasta)
    path(gtfFile)
    
    output:
    tuple val(meta), path('*Log.final.out')   , emit: log_final
    tuple val(meta), path('*Log.out')         , emit: log_out
    tuple val(meta), path('*Log.progress.out'), emit: log_progress
    path  "versions.yml"                      , emit: versions

    tuple val(meta), path('*d.out.bam')              , optional:true, emit: bam
    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*fastq.gz')               , optional:true, emit: fastq
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab
    tuple val(meta), path('*.SJ.out.tab')            , optional:true, emit: spl_junc_tab
    tuple val(meta), path('*.ReadsPerGene.out.tab')  , optional:true, emit: read_per_gene_tab
    tuple val(meta), path('*.out.junction')          , optional:true, emit: junction
    tuple val(meta), path('*.out.sam')               , optional:true, emit: sam
    tuple val(meta), path('*.wig')                   , optional:true, emit: wig
    tuple val(meta), path('*.bg')                    , optional:true, emit: bedgraph

    
    script:
    """
    mkdir star_genome
    STAR runMode genomeGenerate --genomeDir star_genome/ --genomeFastaFiles $genomeFasta --sjdbGTFfile $gtfFile --runThreadN 4
    STAR --readFilesIn ${reads[0]} ${reads[1]} --outFileNamePrefix star_results --genomeDir star_genome/
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        gawk: \$(echo \$(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*\$//')
    END_VERSIONS
    """
    

}