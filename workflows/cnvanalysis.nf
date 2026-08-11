/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_FASTQ              } from '../modules/local/cat_fastq/main.nf'
include { SEQKIT_STATS           } from '../modules/local/seqkit/main.nf'
include { MAPPING                } from '../subworkflows/local/mapping/mapping.nf'
include { CNV_CHECK              } from '../subworkflows/local/cnv_check/cnv_check.nf'
include { SNP_CHECK              } from '../subworkflows/local/snp_check/snp_check.nf'
include { QUARTO_TEXT            } from '../modules/local/quarto/main.nf'
include { QUARTO_SECTION         } from '../modules/local/quarto/main.nf'
include { QUARTO_REPORT          } from '../modules/local/quarto/main.nf'
include { QUARTO_TABLE           } from '../modules/local/quarto/main.nf'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_cnvanalysis_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CNVANALYSIS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()

    CAT_FASTQ (
        ch_samplesheet.map { meta, fastqFiles, ref ->
            tuple(meta, fastqFiles instanceof List ? fastqFiles.flatten() : [fastqFiles])
        }
    )

    //
    // Run seqkit on the input fastq files to get basic stats
    //
    SEQKIT_STATS(
        CAT_FASTQ.out.reads
    )

    //
    // Run Mapping modules
    //
    MAPPING (
       ch_samplesheet
    )

    // Prep reference channel with refID only
    ch_ref = ch_samplesheet.map { meta, fastqFiles, ref -> 
    def refName = file(ref).name.replaceAll(/.fa/, '').replaceAll(/_.*$/, '')
    tuple(meta, refName)
    }

    // Run CNV_check modules
    CNV_CHECK (
       MAPPING.out.bam,
       ch_ref
    )

    // Prep reference channel and bam channel for SNP_check
    ch_refinfo = ch_samplesheet.map { meta, fastqFiles, ref -> 
    def refName = file(ref).name.replaceAll(/.fa/, '').replaceAll(/_.*$/, '')
    tuple(meta, refName, ref, ref + '.fai')
    }
    ch_baminfo = ch_samplesheet.map{item -> item[0]}
        .join(MAPPING.out.bam) //contains both paths to .bam and .bai files
    // Run SNP_check modules
    SNP_CHECK (
       ch_baminfo,
       ch_refinfo
    )

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOFTWARE VERSIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    // Collect versions
    ch_versions = ch_versions
        .mix(SEQKIT_STATS.out.versions)
        .mix(MAPPING.out.versions)
        .mix(CNV_CHECK.out.versions)
        .mix(SNP_CHECK.out.versions)

    ch_samples = CAT_FASTQ.out.reads
    .map { meta, _stats ->
        def meta_reduced = file(params.outdir).name
        tuple(id: meta_reduced)
    }
    .unique()

    // Extract all versions into a single channel of values
    versions = softwareVersionsToYAML(ch_versions)
    // Collapse the channel of versions into a single value
    versions = versions.collect().map { it.join('\n\n') }
    versions = ch_samples
        .combine(versions)

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

    ch_section_description = channel.of("Software Versions")

    ch_section_inputs = QUARTO_TEXT.out.quarto_text
        .combine(ch_section_description)

    QUARTO_SECTION(
        ch_section_inputs
    )
    // // Add the versions to the channel of sections for every report

    ch_sections = ch_sections.mix(QUARTO_SECTION.out.quarto_section)


    
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

    // Emit outputs so the parent workflow can reference them
    emit:
    multiqc_report = Channel.empty()   // placeholder - actual MultiQC step may set this
    versions = versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
