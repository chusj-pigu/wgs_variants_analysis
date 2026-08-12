// Import modules
include { NASVAR_COVERAGE         } from '../../../modules/local/nasvar_dev/main.nf'
include { NASVAR_KARYOTYPE        } from '../../../modules/local/nasvar_dev/main.nf'
include { QUARTO_TABLE_COLNAMES as QUARTO_TABLE1 } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION          } from '../../../modules/local/quarto/main.nf'          // Quarto section

// Define the main workflow
workflow KARYOTYPE {
    take:
    bam         // channel: from mapping workflow, includes index
    ref_json    // reference json file
    repeats     // repeats file
    config      // config file for karyotype analysis

    main:
    ch_versions = Channel.empty() // For collecting version info

    //Prepare NASVAR_COVERAGE input channel
    ch_in_nasvar = bam
        .combine(ref_json)
        .combine(repeats)

    NASVAR_COVERAGE (ch_in_nasvar)

    //Prepare NASVAR_KARYOTYPE input channel
    ch_in_karyotype = NASVAR_COVERAGE.out.cov
        .combine(config)
        .combine(ref_json)
    
    NASVAR_KARYOTYPE (ch_in_karyotype)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        TABLE INPUT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Get only necessary lines from nasvar json output
    def desired_stats = [
    'karyotype_string',
    'iscn_string'
    ]

    ch_nasvar_summary = NASVAR_KARYOTYPE.out.karyo_json
    .map { meta, json -> tuple(meta.id, json) }
    .map { sample, json ->
        def project = file(params.outdir).name
        // Parse nasvar's "key\tvalue" lines into a lookup map
        def karyo_res = new groovy.json.JsonSlurper().parse(json)
        def karyo_data = karyo_res.karyotype ?: [:]
        tuple(project, sample, karyo_data)
    }
    .collectFile(sort: true) { project, sample, karyo_res ->
        def row = ([sample] + desired_stats.collect { karyo_res[it] ?: 'NA' }).join('\t')
        return [ "${project}_karyotype_results.tsv", "${row}\n" ]
    }

    ch_quarto_table = ch_nasvar_summary
        .map { table ->
            def meta_project = table.name.replace('_karyotype_results.tsv', '')
            tuple(id:meta_project,table)                    // Convert meta project to meta id
        }
        .map { meta, table ->
            def caption = "Summary karyotype stats for ${meta.id} on BAM files (no filtering)"
            def col_names = "Sample, Karyotype, ISCN"
            def section = "Karyotype"
            def process = "karyotype-${meta.id}"
            tuple(meta, table, caption, col_names, section, process)
        }
    

    QUARTO_TABLE1(
        ch_quarto_table
        )

    ch_section_inputs = QUARTO_TABLE1.out.quarto_table


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        QUARTO_SECTION
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_section_description = channel.of("Karyotype results from NASVAR")

    ch_section_inputs = ch_section_inputs
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }

    QUARTO_SECTION(
        ch_section_inputs,
        ch_section_description
    )

    ch_section = QUARTO_SECTION.out.quarto_section

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_versions = NASVAR_COVERAGE.out.versions
                    .mix(NASVAR_KARYOTYPE.out.versions)
                    .mix(QUARTO_TABLE1.out.versions)
                    .mix(QUARTO_SECTION.out.versions)

    emit:
    nasvar_cov         = NASVAR_COVERAGE.out.cov
    nasvar_karyo       = NASVAR_KARYOTYPE.out.karyo_json
    section            = QUARTO_SECTION.out.quarto_section
    versions           = ch_versions

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/