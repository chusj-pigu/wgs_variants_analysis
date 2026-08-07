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

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'cnvanalysis_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    
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
