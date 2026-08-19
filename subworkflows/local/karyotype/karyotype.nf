// Import modules
include { NASVAR_COVERAGE         } from '../../../modules/local/nasvar_dev/main.nf'
include { NASVAR_KARYOTYPE        } from '../../../modules/local/nasvar_dev/main.nf'
include { QUARTO_TABLE_COLNAMES as QUARTO_TABLE1 } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION          } from '../../../modules/local/quarto/main.nf'          // Quarto section

// Define the main workflow
workflow KARYOTYPE {
    take:
    bam         // channel: from mapping workflow, includes index
    ch_ref      // ch with necessary references to run nasvar
    mean_cov    // mean coverage from mapping subworkflow

    main:
    ch_versions = Channel.empty() // For collecting version info

    //Prepare NASVAR_COVERAGE input channel
    ch_ref_coverage = ch_ref.map { meta, ref, karyotype_ref, karyotype_repeats, karyotype_config, karyotype_cov_cutoff -> 
            tuple(meta, karyotype_ref, karyotype_repeats)
        }
    
    ch_in_nasvar = bam
        .join(ch_ref_coverage)

    NASVAR_COVERAGE (ch_in_nasvar)

    //Prepare NASVAR_KARYOTYPE input channel
    ch_ref_karyotype = ch_ref.map { meta, ref, karyotype_ref, karyotype_repeats, karyotype_config, karyotype_cov_cutoff -> 
            tuple(meta, karyotype_config, karyotype_ref)
        }

    ch_in_karyotype = NASVAR_COVERAGE.out.cov
        .join(ch_ref_karyotype)
    
    NASVAR_KARYOTYPE (ch_in_karyotype)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        TABLE INPUT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Get only necessary lines from nasvar json output
    def desired_stats = [ 'karyotype_string', 'iscn_string' ]

    // Join karyotype results with mean coverage per-sample (keyed by meta.id)
    ch_karyo_with_cov = NASVAR_KARYOTYPE.out.karyo_json
        .map { meta, json -> tuple(meta.id, meta, json) }
        .join( mean_cov.map { meta, cov -> tuple(meta.id, cov) } )
        .join( ch_ref.map { meta, ref, karyotype_ref, karyotype_repeats, karyotype_config, karyotype_cov_cutoff -> tuple(meta.id, ref, karyotype_cov_cutoff) } )
        .map { id, meta, json, cov, ref, cutoff -> tuple(meta, json, ref, cov, cutoff) }

    ch_nasvar_summary = ch_karyo_with_cov
        .map { meta, json, ref, cov, cutoff ->
            def project = file(params.outdir).name
            def refname = file(ref).name.replaceAll(/\.(fa|fasta)(\.gz)?$/, '')
            // Parse nasvar's json output into a lookup map
            def karyo_result = new groovy.json.JsonSlurper()
                .parse(json.toFile())
                .karyotype
                ?: [:]
            tuple(project, meta.id, refname, karyo_result, cov, cutoff)
        }
        .collectFile(sort: true) { List row ->
            def (project, sampleId, refname, karyo_result, cov, cutoff) = row
            def covStr   = cov?.toString()?.trim()
            def covValue = (covStr == null || covStr == 'NA' || covStr.isEmpty())
                ? null
                : covStr.replaceAll(/[^\d.]/, '') ?: null
            def isLowCov = (covValue == null) || ((covValue as Double) < (cutoff as Double))

            def isMm10 = (refname == "mm10")

            def fields = isMm10
                ? desired_stats.collect { "not_available"}
                : isLowCov
                    ? desired_stats.collect { "low_coverage(<${cutoff})" }
                    : desired_stats.collect { karyo_result[it]?.toString() ?: 'NA' }

            def tsvRow = ([sampleId] + [refname] + fields).join('\t')

            return [ "${project}_karyotype_results.tsv", "${tsvRow}\n" ]
    
        }

    ch_quarto_table = ch_nasvar_summary
        .map { table ->
            def meta_project = table.name.replace('_karyotype_results.tsv', '')
            tuple(id:meta_project,table)                    // Convert meta project to meta id
        }
        .map { meta, table ->
            def caption = "Summary karyotype stats for ${meta.id} on BAM files (no filtering)"
            def col_names = "Sample, Ref, Karyotype, ISCN"
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