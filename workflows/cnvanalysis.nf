/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MAPPING                } from '../subworkflows/local/mapping/mapping.nf'
include { CNV_CHECK              } from '../subworkflows/local/cnv_check/cnv_check.nf'
include { SNP_CHECK              } from '../subworkflows/local/snp_check/snp_check.nf'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_cnvanalysis_pipeline'
include { VARIANTS_REPORT        } from '../subworkflows/local/report/variants_report.nf'
include { MIDNIGHT_REPORT        } from '../subworkflows/local/report/final_report.nf'

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

    //
    // Collate and save software versions
    //
    ch_versions = ch_versions
        .mix(MAPPING.out.versions)
        //.mix(CNV_CHECK.out.versions)
        //.mix(SNP_CHECK.out.versions)
    
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE REPORT
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    VARIANTS_REPORT(
            MAPPING.out.coverage
        )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        COLLECT SECTIONS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Collect sections from all analysis steps
    ch_sections = VARIANTS_REPORT.out.sections
    
    ch_mode = channel.of("WGS")

    ch_id = ch_samplesheet.map{meta, fastqFiles, ref -> meta}

    
    MIDNIGHT_REPORT(
        ch_id,
        ch_sections,
        ch_versions,
        ch_mode
    )


    // Emit outputs so the parent workflow can reference them
    emit:
    multiqc_report = Channel.empty()   // placeholder - actual MultiQC step may set this
    versions = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
