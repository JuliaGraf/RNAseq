include { FASTQC } from '../modules/fastqc'

workflow RNASEQ {

    main:
    ch_versions    = Channel.empty()
    ch_samplesheet = Channel.value(file(params.input, checkIfExists: true))

    // 1. Read samplesheet
    ch_samplesheet
        .splitCsv(header: true)
        .map { it -> [[sample: it.sample, strandedness: it.strandedness], [it.fastq_1, it.fastq_2]]}
        .set { in_ch }
    in_ch.view()
    // 2. Quality control
    FASTQC(in_ch)

}