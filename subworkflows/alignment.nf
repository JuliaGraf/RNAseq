include { STAR_ALIGN }          from '../modules/star_align.nf'
include { STAR_GENOMEGENERATE } from '../modules/star_genomegenerate'

workflow ALIGNMENT {
    take:
    ch_trimmed
    fasta
    gtf

    main:
    STAR_GENOMEGENERATE(fasta, gtf)
    ch_versions = STAR_GENOMEGENERATE.out.versions
    STAR_ALIGN(ch_trimmed, gtf, STAR_GENOMEGENERATE.out.index)
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    emit:
    versions = ch_versions
    bam      = STAR_ALIGN.out.bam
    fai      = STAR_GENOMEGENERATE.out.fai
}