//
// Subworkflow with functionality specific to the nf-core/cnvanalysis pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    if (params.karyotype) {
        Channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map { meta, path_to_fastq, ref, karyotype_ref, karyotype_repeats, karyotype_config, karyotype_cov_cutoff ->
                def refPath      = validateRef(meta.id, ref)
                def fastqFiles   = validateFastqDir(meta.id, path_to_fastq)
                def ktRefPath    = validateKaryotypeFile(meta.id, karyotype_ref,     'karyotype_ref',     '.json')
                def ktRepPath    = validateKaryotypeFile(meta.id, karyotype_repeats,  'karyotype_repeats', '.bed')
                def ktConfigPath = validateKaryotypeFile(meta.id, karyotype_config,   'karyotype_config',  '.json')
                def isPaired     = fastqFiles.any { it =~ /(?i)(_R?2|_2)\./ }
                meta             = meta + [ single_end: !isPaired ]
                return [ meta.id, meta, fastqFiles, refPath, ktRefPath, ktRepPath, ktConfigPath, karyotype_cov_cutoff ]
            }
            .groupTuple()
            .map { validateInputSamplesheet(it, true) }
            .set { ch_samplesheet }
    } else {
        Channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map { row ->
                def meta          = row[0]
                def path_to_fastq = row[1]
                def ref           = row[2]
                def refPath    = validateRef(meta.id, ref)
                def fastqFiles = validateFastqDir(meta.id, path_to_fastq)
                def isPaired   = fastqFiles.any { it =~ /(?i)(_R?2|_2)\./ }
                meta           = meta + [ single_end: !isPaired ]
                return [ meta.id, meta, fastqFiles, refPath ]
            }
            .groupTuple()
            .map { validateInputSamplesheet(it, false) }
            .set { ch_samplesheet }
    }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                []
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate and resolve a reference FASTA path (must exist and have a .fai index)
// Returns the resolved path string.
//
def validateRef(sampleId, ref) {
    def refPath = ref.toString()
    def refFile = file(refPath)
    if (!refFile.exists()) {
        error("Reference file for sample '${sampleId}' does not exist: ${refPath}")
    }
    def r = refPath.toLowerCase()
    if (!(r.endsWith('.fa') || r.endsWith('.fasta') || r.endsWith('.fa.gz') || r.endsWith('.fasta.gz'))) {
        error("Reference file for sample '${sampleId}' must be a FASTA file (.fa/.fasta, optionally .gz): ${refPath}")
    }
    def faiFile = file("${refPath}.fai")
    if (!faiFile.exists()) {
        error("Index file (.fai) for reference of sample '${sampleId}' does not exist. " +
              "Please run: samtools faidx ${refPath}")
    }
    return refPath
}

//
// Validate a FastQ directory and return the list of FastQ file paths it contains.
//
def validateFastqDir(sampleId, path_to_fastq) {
    def fqDirPath = path_to_fastq.toString()
    def fqDir     = file(fqDirPath)
    if (!fqDir.exists() || !fqDir.isDirectory()) {
        error("FastQ directory for sample '${sampleId}' does not exist or is not a directory: ${fqDirPath}")
    }
    def fastqFiles = fqDir.list()
        .findAll { name ->
            def n = name.toLowerCase()
            n.endsWith('.fq') || n.endsWith('.fq.gz') || n.endsWith('.fastq') || n.endsWith('.fastq.gz')
        }
        .collect { name -> fqDir.resolve(name).toString() }
    if (fastqFiles.size() == 0) {
        error("No FastQ files found in directory for sample '${sampleId}': ${fqDirPath}")
    }
    return fastqFiles
}

//
// Validate a karyotype-specific file (must exist and match the expected extension).
// Returns the resolved path string.
//
def validateKaryotypeFile(sampleId, filePath, fieldName, expectedExt) {
    def pathStr  = filePath.toString()
    def fileObj  = file(pathStr)
    if (!fileObj.exists()) {
        error("Karyotype field '${fieldName}' for sample '${sampleId}' does not exist: ${pathStr}")
    }
    if (!pathStr.toLowerCase().endsWith(expectedExt)) {
        error("Karyotype field '${fieldName}' for sample '${sampleId}' must be a ${expectedExt} file " +
              "(see NASVAR documentation): ${pathStr}")
    }
    return pathStr
}

//
// Validate channels from input samplesheet
//
//
// Validate channels from input samplesheet after groupTuple().
// Checks that multiple runs of the same sample share the same endedness.
// Returns the tuple emitted into ch_samplesheet.
//
// Non-karyotype:  [ meta, fastqFiles, refPath ]
// Karyotype:      [ meta, fastqFiles, refPath, ktRefPath, ktRepPath, ktConfigPath, ktCovCutoff ]
//
def validateInputSamplesheet(input, Boolean withKaryotype) {
    def metas = input[1]
    if (metas.collect { it.single_end }.unique().size() != 1) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same " +
              "datatype (single-end or paired-end): ${metas[0].id}")
    }

    def meta        = metas[0]
    def fastqFiles  = input[2].flatten()
    def refPath     = input[3][0]

    if (!withKaryotype) {
        return [ meta, fastqFiles, refPath ]
    }

    def ktRefPath    = input[4][0]
    def ktRepPath    = input[5][0]
    def ktConfigPath = input[6][0]
    def ktCovCutoff  = input[7][0]
    return [ meta, fastqFiles, refPath, ktRefPath, ktRepPath, ktConfigPath, ktCovCutoff ]
}

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
