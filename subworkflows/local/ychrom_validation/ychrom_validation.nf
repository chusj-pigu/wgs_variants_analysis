// Import modules
include { SAMTOOLS_SPLIT_BY_BED  } from '../../../modules/local/samtools/main.nf'         // Split BAM by BED file
include { SAMTOOLS_INDEX         } from '../../../modules/local/samtools/main.nf'         // Index BAM file
include { MOSDEPTH_GENERAL       } from '../../../modules/local/mosdepth/main.nf'         // Compute coverage stats with mosdepth 
include { QUARTO_TABLE_COLNAMES as QUARTO_TABLE1           } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION         } from '../../../modules/local/quarto/main.nf'

// Define the main workflow
workflow YCHROM_VALIDATION {
    take:
    ch_bam
    mean_cov    // mean coverage from mapping subworkflow

    main:
    ch_versions = Channel.empty() // For collecting version info

    ch_msy = channel.fromPath(params.msy_locus, checkIfExists:true)


    // Split BAM by MSY locus
    ch_to_split = ch_bam.combine(ch_msy)

    SAMTOOLS_SPLIT_BY_BED (ch_to_split)
    SAMTOOLS_INDEX (SAMTOOLS_SPLIT_BY_BED.out.panel)


    MOSDEPTH_GENERAL(SAMTOOLS_INDEX.out.bamfile_index)

    ch_coverage_msy = MOSDEPTH_GENERAL.out.summary
        .map { meta, table ->
        // Read the file content as a list of lines
            def lines = table.readLines()
            def coverage = lines.last().tokenize('\t')[3].toFloat()    // Last line and only take mean coverage column (4th)
            tuple(meta, coverage)
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        TABLE INPUT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_table_input = ch_coverage_msy
        .join(mean_cov)
        .map { meta, coverage, mean_cov_value ->
            def mc = mean_cov_value.toString().toFloat()
            def ratio = coverage / mc
            def sex_call = ratio > 0.4 ? "Male" : (ratio >= 0.2 ? "Probable Male" : "Female")
            tuple(meta, [meta.id, coverage, mc, ratio, sex_call])
        }

    ch_ychrom_summary = ch_table_input
        .map { meta, row ->
            def project = file(params.outdir).name
            tuple(project, row)
        }
        .collectFile(sort: true) { project, row ->
            def line = row.join('\t')
            return [ "${project}_ychrom_validation.tsv", "${line}\n" ]
        }

    ch_quarto_table = ch_ychrom_summary
        .map { table ->
            def meta_project = table.name.replace('_ychrom_validation.tsv', '')
            tuple(id: meta_project, table)
        }
        .map { meta, table ->
            def caption = "Y chromosome presence validation using MSY locus coverage"
            def col_names = "Sample, MSY Coverage, Mean Coverage, Ratio, Sex Call"
            def section = "chrY_Validation"
            def process = "ychrom-validation"
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

    ch_section_description = channel.value("Y chromosome presence validation using MSY locus coverage. The mean coverage of the MSY locus is calculated and compared to the mean coverage of the run to determine if the sample is male or female. [Ratio < 0.2 : Female, 0.2 <= Ratio < 0.4 : Probable Male, Ratio >= 0.4 : Male]")

    ch_section_inputs = ch_section_inputs
        .map { meta, section, filePaths -> tuple([id: file(params.outdir).name], section, filePaths) }
        .groupTuple()
        .map { meta, section, filePaths -> [meta, section[0], filePaths] }
        .combine(ch_section_description)

    QUARTO_SECTION(
        ch_section_inputs
    )

    ch_section = QUARTO_SECTION.out.quarto_section

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT VERSIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_versions = SAMTOOLS_SPLIT_BY_BED.out.versions
                    .mix(MOSDEPTH_GENERAL.out.versions)
                    .mix(QUARTO_TABLE1.out.versions)
                    .mix(QUARTO_SECTION.out.versions)

    emit:
    chry_cov         = ch_coverage_msy
    chry_summary     = MOSDEPTH_GENERAL.out.summary
    chry_dist        = MOSDEPTH_GENERAL.out.dist
    section          = QUARTO_SECTION.out.quarto_section
    versions         = ch_versions

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/