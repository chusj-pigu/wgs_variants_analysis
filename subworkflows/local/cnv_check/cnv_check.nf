// Import modules
include { QDNASEQ_CALL         } from '../../../modules/local/qdnaseq/main.nf'         // QDNAseq CNV calling

// Define the main workflow
workflow CNV_CHECK {
    take:
    bam         // channel: from mapping workflow, includes index
    ref         // reference ID (hg19 or hg38)

    main:
    ch_versions = Channel.empty() // For collecting version info

    //Prepare QDNAseq input channel
    ch_in_qdnaseq = bam
        .join(ref)

    // Identify CNV with QDNAseq
    QDNASEQ_CALL (ch_in_qdnaseq)


    ch_versions = QDNASEQ_CALL.out.versions

    emit:
    qdnaseq_vcf         = QDNASEQ_CALL.out.call_vcf
    qdnaseq_plot        = QDNASEQ_CALL.out.cov_png
    versions            = ch_versions

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/