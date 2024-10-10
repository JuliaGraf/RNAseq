include { FASTQC }                    from '../modules/fastqc'
include { TRIMGALORE }                from '../modules/trimgalore'
include { STAR_ALIGN }                from '../modules/star_align.nf'
include { STAR_GENOMEGENERATE }       from '../modules/star_genomegenerate'
include { MARKDUPLICATES }            from '../modules/markduplicates'
include { RSEM_PREPARE_REFERENCE }    from '../modules/rsem_prepare_reference'
include { RSEM_CALCULATE_EXPRESSION } from '../modules/rsem_calculate_expression'


import org.yaml.snakeyaml.Yaml

workflow RNASEQ {

    main:
    ch_versions    = Channel.empty()
    ch_samplesheet = Channel.value(file(params.input, checkIfExists: true))

    // 1. Read samplesheet
    ch_samplesheet
        .splitCsv ( header:true, sep:',' )
        .map { it -> 
                if (!it.fastq_2) {
                    return [[sample: it.sample, strandedness: it.strandedness], [file(it.fastq_1, checkIfExists: true) ]]
                } else {
                    return [[sample: it.sample, strandedness: it.strandedness],[file(it.fastq_1, checkIfExists: true), file(it.fastq_2, checkIfExists: true)]]
                }
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
    STAR_GENOMEGENERATE(file(params.genomeFasta, checkIfExists:true),file(params.gtfFile, checkIfExists:true))
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    STAR_ALIGN(ch_trimmed,file(params.gtfFile, checkIfExists:true), STAR_GENOMEGENERATE.out.index)
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    // 5. QUANTIFICATION
    RSEM_PREPARE_REFERENCE(file(params.genomeFasta, checkIfExists:true),file(params.gtfFile, checkIfExists:true))
    RSEM_CALCULATE_EXPRESSION(ch_trimmed, RSEM_PREPARE_REFERENCE.out.out_dir)

    // 6. Mark Duplicates
    MARKDUPLICATES(STAR_ALIGN.out.bam, file(params.genomeFasta, checkIfExists:true), STAR_GENOMEGENERATE.out.fai)
    ch_versions = ch_versions.mix(MARKDUPLICATES.out.versions.first())

    // Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "results/pipeline_info",
            name: 'pipeline_software_versions.yml',
            sort: true,
            newLine: true
        )
}


//
// Get software versions for pipeline
//
def processVersionsFromYAML(yaml_file) {
    Yaml yaml = new Yaml()
    versions = yaml.load(yaml_file).collectEntries { k, v -> [ k.tokenize(':')[-1], v ] }
    return yaml.dumpAsMap(versions).trim()
}

//
// Generate workflow version string
//
def getWorkflowVersion() {
    String version_string = ""
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }

    if (workflow.commitId) {
        def git_shortsha = workflow.commitId.substring(0, 7)
        version_string += "-g${git_shortsha}"
    }

    return version_string
}

//
// Get workflow version for pipeline
//
def workflowVersionToYAML() {
    return """
    Workflow:
        $workflow.manifest.name: ${getWorkflowVersion()}
        Nextflow: $workflow.nextflow.version
    """.stripIndent().trim()
}

//
// Get channel of software versions used in pipeline in YAML format
//
def softwareVersionsToYAML(ch_versions) {
    return ch_versions
                .unique()
                .map { processVersionsFromYAML(it) }
                .unique()
                .mix(Channel.of(workflowVersionToYAML()))
}
