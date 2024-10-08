include { FASTQC } from '../modules/fastqc'
include {TRIMGALORE} from '../modules/trimgalore'

workflow RNASEQ {

    main:
    ch_versions    = Channel.empty()
    ch_samplesheet = Channel.value(file(params.input, checkIfExists: true))

    // 1. Read samplesheet
    ch_samplesheet
        .splitCsv ( header:true, sep:',' )
        .map { it -> [[sample: it.sample, strandedness: it.strandedness],
                     [file(it.fastq_1, checkIfExists: true), file(it.fastq_2, checkIfExists: true)]]
        }
        .set { in_ch }

    // 2. Quality Control
    FASTQC(in_ch)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // 3. Read Trimming
    TRIMGALORE(in_ch)
    ch_trimmed = TRIMGALORE.out.trimmed
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())

    // 4. Alignment


}
