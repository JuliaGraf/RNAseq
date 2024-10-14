// Import modules
include { FASTQC }                    from '../modules/fastqc'
include { TRIMGALORE }                from '../modules/trimgalore'
include { MARKDUPLICATES }            from '../modules/markduplicates'

// Import subworkflows
include { ALIGNMENT }                 from '../subworkflows/alignment.nf'
include { QUANTIFICATION }            from '../subworkflows/quantification.nf'

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
                    return [[sample: it.sample], [file(it.fastq_1, checkIfExists: true) ]]
                } else {
                    return [[sample: it.sample],[file(it.fastq_1, checkIfExists: true), file(it.fastq_2, checkIfExists: true)]]
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

    // 4. SUBWORKFLOW: Alignment
    ALIGNMENT(ch_trimmed,
              file(params.genomeFasta, checkIfExists:true),
              file(params.gtfFile, checkIfExists:true))
    ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

    // 5. SUBWORKFLOW: Quantification
    QUANTIFICATION(ch_trimmed,
                   file(params.genomeFasta, checkIfExists:true),
                   file(params.gtfFile, checkIfExists:true))
    ch_versions = ch_versions.mix(QUANTIFICATION.out.versions)

    // 6. Mark Duplicates
    MARKDUPLICATES(ALIGNMENT.out.bam,
                   file(params.genomeFasta, checkIfExists:true),
                   ALIGNMENT.out.fai)
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
