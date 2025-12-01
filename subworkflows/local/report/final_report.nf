/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUARTO_TEXT } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION } from '../../../modules/local/quarto/main.nf'
include { QUARTO_FIGURE } from '../../../modules/local/quarto/main.nf'
include { QUARTO_REPORT } from '../../../modules/local/quarto/main.nf'
include { softwareVersionsToYAML } from '../../../subworkflows/nf-core/utils_nfcore_pipeline'

workflow MIDNIGHT_REPORT {

    take:
    ch_id
    ch_sections
    ch_versions
    ch_mode

    main:


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOFTWARE VERSIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // Extract all versions into a single channel of values
    versions = softwareVersionsToYAML(ch_versions).view{it -> "versions = $it"}
    // Collapse the channel of versions into a single value
    versions = versions.collect().map { it.join('\n\n') }
    versions = ch_id.combine(versions)

    // Give it an ID of versions
    versions = versions
        .map {
            versions_out ->
            def section = "Versions"
            def process = "versions"

            [versions_out[0], versions_out[1]] + [section, process]
            }

    QUARTO_TEXT(
        versions
        )

    ch_section_inputs = QUARTO_TEXT.out.quarto_text

    QUARTO_SECTION(
        ch_section_inputs,
        "Software Versions"
    )

    ch_section_vers = QUARTO_SECTION.out.quarto_section
        .map{meta, section, filePaths, reports ->
            tuple (id:meta.id, section, filePaths, reports)
        }

    // // Add the versions to the channel of sections for every report

    ch_sections = ch_sections.mix(ch_section_vers)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INDIVIDUAL REPORTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    ch_report_sections = ch_sections
        .groupTuple()
        .map { id, section, filePaths, reports ->
            tuple(id, section, filePaths, reports)
        }

    ch_title = ch_report_sections
        .combine(ch_mode)
        .map { meta, section, filePaths, reports, mode ->
            def title = "Variant Analysis Pipeline ${mode} Report"
            tuple(meta, title) }

    ch_subtitle = ch_report_sections
        .combine(ch_mode)
        .map { meta, section, filePaths, reports, mode ->
            def subtitles = "Outputs for the ${mode} branch of the Variant Analysis pipeline"
            tuple(meta, subtitles)}

    ch_template = channel.fromPath(params.report_template)

    ch_report_in = ch_report_sections
        .combine(ch_template)
        .join(ch_title)
        .join(ch_subtitle)
        .map { meta, section, filePaths, reports, template, title, subtitle ->
            tuple(
                tuple(meta, section, filePaths, reports),  
                template,                                   
                title,                                      
                subtitle                                   
            )
        }

    ch_report_in.view{it -> "ch_report_in = $it"}

    QUARTO_REPORT(
        ch_report_in.map { it[0] },  
        ch_report_in.map { it[1] },  
        ch_report_in.map { it[2] },  
        ch_report_in.map { it[3] }   
    )

    ch_report = QUARTO_REPORT.out.report


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    emit:
    ch_report
}