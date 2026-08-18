// Import modules
include { MINIMAP2_ALIGN         } from '../../../modules/local/minimap2/main.nf'         // minimap2 alignment
include { SAMTOOLS_TOBAM         } from '../../../modules/local/samtools/main.nf'         // Convert SAM to BAM
include { SAMTOOLS_SORT          } from '../../../modules/local/samtools/main.nf'         // Sort BAM
include { SAMTOOLS_INDEX         } from '../../../modules/local/samtools/main.nf'         // Index BAM
include { CRAMINO_STATS          } from '../../../modules/local/cramino/main.nf'          // Coverage stats
include { QUARTO_TABLE_COLNAMES as QUARTO_TABLE1           } from '../../../modules/local/quarto/main.nf'
include { QUARTO_SECTION         } from '../../../modules/local/quarto/main.nf'

// Define the main workflow
workflow MAPPING {
    take:
    ch_samplesheet

    main:
    ch_versions = Channel.empty() // For collecting version info

    // Align reads to reference genome
    MINIMAP2_ALIGN (ch_samplesheet)
    
    // Convert SAM to BAM
    SAMTOOLS_TOBAM(MINIMAP2_ALIGN.out.sam)

    // Sort and index BAM
    SAMTOOLS_SORT(SAMTOOLS_TOBAM.out.bamfile)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.sortedbam)

    // Compute coverage stats
    CRAMINO_STATS(SAMTOOLS_INDEX.out.bamfile_index)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TABLE INPUT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // Get only necessary lines from cramino summary

    def desired_stats = [
    'Number of alignments',
    '% from total alignments',
    'Yield [Gb]',
    'Mean coverage',
    'N50',
    'Mean length',
    'Mean identity'
    ]

    ch_cramino_parsed = CRAMINO_STATS.out.stats
        .map { meta, table ->
            def stats = table.readLines()
                .findAll { it.contains('\t') }
                .collectEntries { line ->
                    def parts = line.split('\t', 2)
                    [(parts[0]): parts[1]]
                }
            tuple(meta, stats)
        }

    ch_mean_coverage = ch_cramino_parsed
        .map { meta, stats -> tuple(meta, stats['Mean coverage'] ?: 'NA') }

    ch_cramino_summary = ch_cramino_parsed
        .join (ch_samplesheet)
        .map { meta, stats, fastqfiles, ref ->
            def project = file(params.outdir).name
            def refname = file(ref).name.replaceAll(/\.(fa|fasta)(\.gz)?$/, '')
            tuple(project, meta.id, refname, stats)
        }
        .collectFile(sort: true) { project, sample, refname, stats ->
            def row = ([sample] + [refname] + desired_stats.collect { stats[it] ?: 'NA' }).join('\t')
            return [ "${project}_mapping_stats.tsv", "${row}\n" ]
        }

    ch_quarto_table = ch_cramino_summary
        .map { table ->
            def meta_project = table.name.replace('_mapping_stats.tsv', '')
            tuple(id: meta_project, table)
        }
        .map { meta, table ->
            def caption   = "Summary mapping stats for ${meta.id} on BAM files (no filtering)"
            def col_names = "Sample, Ref, # Alignments, % from total alignments, Yield [Gb], Mean Coverage, N50, Mean length, Mean identity"
            def section   = "Mapping_QC"
            def process   = "mapping-qc-${meta.id}"
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

    ch_section_description = channel.of("Summary Statistics for BAM files")

    ch_section_inputs = ch_section_inputs
        .groupTuple()
        .map { id, section, filePaths ->
            [id, section[0], filePaths]
        }.combine(ch_section_description)

    QUARTO_SECTION(
        ch_section_inputs
    )

    ch_section = QUARTO_SECTION.out.quarto_section

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT VERSIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    // Collect versions from all modules
    ch_versions = MINIMAP2_ALIGN.out.versions
        .mix(SAMTOOLS_TOBAM.out.versions)
        .mix(SAMTOOLS_SORT.out.versions)
        .mix(SAMTOOLS_INDEX.out.versions)
        .mix(CRAMINO_STATS.out.versions)
        .mix(QUARTO_TABLE1.out.versions)
        .mix(QUARTO_SECTION.out.versions)

    emit:
    bam      = SAMTOOLS_INDEX.out.bamfile_index       // Final sorted BAM with index
    coverage = CRAMINO_STATS.out.stats                // Coverage stats
    mean_cov = ch_mean_coverage                       // Mean coverage for each sample
    section  = QUARTO_SECTION.out.quarto_section      // Quarto section for mapping stats
    versions = ch_versions                            // All tool versions

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/