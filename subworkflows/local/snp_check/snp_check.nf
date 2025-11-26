// Import modules
include { SNIFFLES_CALL         } from '../../../modules/local/sniffles/main.nf'         // Sniffles SNP calling

// Define the main workflow
workflow SNP_CHECK {
    take:
    baminfo         // channel: tuples meta, bam_path, bai_path
    refinfo         // channel: tuples meta, ref_ID, ref_fastq_path, ref_fai_path

    main:
    ch_versions = Channel.empty() // For collecting version info

    //Prepare Sniffles input channel
    ch_in_sniffles = baminfo
        .join(refinfo)

    // Identify CNV with QDNAseq
    SNIFFLES_CALL(ch_in_sniffles)


    ch_versions = SNIFFLES_CALL.out.versions

    emit:
    sniffles_vcf         = SNIFFLES_CALL.out.vcf
    versions            = ch_versions

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/