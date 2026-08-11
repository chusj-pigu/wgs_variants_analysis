// Import modules
include { CAT_FASTQ              } from '../../../modules/local/cat_fastq/main.nf'
include { SEQKIT_STATS           } from '../../../modules/local/seqkit/main.nf'
include { QUARTO_SECTION                } from '../../../modules/local/quarto/main.nf'
include { QUARTO_TABLE_COLNAMES as QUARTO_TABLE1 } from '../../../modules/local/quarto/main.nf'

// Define the main workflow
workflow QC {
    take:
    ch_samplesheet

    main:
    ch_versions = Channel.empty() // For collecting version info

    CAT_FASTQ (
        ch_samplesheet
    )

    //
    // Run seqkit on the input fastq files to get basic stats
    //
    SEQKIT_STATS(
        CAT_FASTQ.out.reads
    )


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        TABLE INPUT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Get only necessary columns from seqkit table
    ch_seqkit_summary = SEQKIT_STATS.out.stats
        .map { meta, table -> tuple(meta.id, table) }
        .map { sample, table ->
            def project = file(params.outdir).name
            //Read the file content as a list of lines
            def lines = table.readLines()
            // Second line and 4th column, in million of reads
            def n_reads = (lines[1].tokenize('\t')[3].toDouble() / 1000000) .toString()
            // In billion of bases
            def n_bases = (lines[1].tokenize('\t')[4].toDouble() / 1000000000) .toString()
            def n50 =  lines[1].tokenize('\t')[12]
            def av_qual =  lines[1].tokenize('\t')[16]
            def av_len =  lines[1].tokenize('\t')[6]
            def max_len =  lines[1].tokenize('\t')[7]
            def q1 =  lines[1].tokenize('\t')[8]
            def q2 =  lines[1].tokenize('\t')[9]
            def q3 =  lines[1].tokenize('\t')[10]
            def q20 = lines[1].tokenize('\t')[14]
            def q30 = lines[1].tokenize('\t')[15]
            def gc = lines[1].tokenize('\t')[17]
           tuple(project,sample,n_reads,n_bases,n50,av_qual,av_len,max_len,q1,q2,q3,q20,q30,gc)
            }
        .collectFile { table ->
            def content = [table[1] + '\t' + table[2] + '\t' + table[3] +
             '\t' + table[4] + '\t' + table[5] + '\t' + table[6] + '\t' +
             table[7] + '\t' + table[8] + '\t' + table[9] + '\t' + table[10] +
            '\t' + table[11] + '\t' + table[12] + '\t' + table[13]
            ].join('\n')
            return [ "${table[0]}_reads_stats.tsv", content + '\n' ]
        }

    ch_quarto_table = ch_seqkit_summary
        .map { table ->
            def meta_project = table.name.replace('_reads_stats.tsv', '')
            tuple(id:meta_project,table)                    // Convert meta project to meta id
        }
        .map { meta, table ->
            def caption = "Summary sequence stats for ${meta.id} on FASTQ files (no filtering)"
            def col_names = "Sample, # Reads (M), # Bases (GB), N50, Mean QS, Mean len, Max len, Q1, Q2, Q3, Q20 (%), Q30 (%), GC (%)"
            def section = "Reads_QC"
            def process = "reads-qc-${meta.id}"
            tuple(meta, table, caption, col_names, section, process)
        }


    QUARTO_TABLE1(
        ch_quarto_table
        )

    ch_section_inputs = QUARTO_TABLE1.out.quarto_table

// /*
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//     QUARTO_SECTION
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// */

    ch_section_description = channel.of("Summary Statistics for FASTQ files")

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

    ch_versions = QUARTO_TABLE1.out.versions
    ch_versions = ch_versions.mix(QUARTO_SECTION.out.versions, SEQKIT_STATS.out.versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    emit:
    reads      = CAT_FASTQ.out.reads                   // concatenated reads per sample
    stats      = SEQKIT_STATS.out.stats                // seqkit stats
    section    = ch_section                            // quarto section for reads QC
    versions   = ch_versions                           // All tool versions

}