include { RSEM_PREPARE_REFERENCE }    from '../modules/rsem_prepare_reference'
include { SALMON_INDEX }              from '../modules/salmon_index'
include { SALMON_QUANTIFICATION }     from '../modules/salmon_quantification'

workflow QUANTIFICATION {
    take:
    ch_trimmed
    fasta
    gtf

    main:
    RSEM_PREPARE_REFERENCE(fasta, gtf)
    ch_versions = RSEM_PREPARE_REFERENCE.out.versions
    SALMON_INDEX(RSEM_PREPARE_REFERENCE.out.index)
    ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
    SALMON_QUANTIFICATION(ch_trimmed,SALMON_INDEX.out.index, gtf)
    ch_versions = ch_versions.mix(SALMON_QUANTIFICATION.out.versions.first())

    emit:
    versions = ch_versions
}