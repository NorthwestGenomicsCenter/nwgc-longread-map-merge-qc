include { MAP_ONT_FASTQ } from '../modules/map_ont_fastq.nf'
include { MERGE_MAPPED_BAMS } from '../modules/merge_mapped_bams.nf'
include { CHECKSUM_BAM } from '../modules/checksum_bam.nf'

workflow ONT_MAP_MERGE_FASTQS {

    main:
        ontFastqs = Channel.empty()
        for (ontFastqFolder in params.ontFastqFolders) {
            fastqFolderChannel = Channel.fromPath(ontFastqFolder + "/*.fastq.gz")
            ontFastqs = ontFastqs.concat(fastqFolderChannel)
        }
        MAP_ONT_FASTQ(ontFastqs)

        def outFolder = "${params.sampleDirectory}"
        def outPrefix = "${params.sampleId}.${params.sequencingTarget}"
        // Merge
        MERGE_MAPPED_BAMS(MAP_ONT_FASTQ.out.mapped_bam.collect(), outFolder, outPrefix)

        // checksum
        CHECKSUM_BAM(MERGE_MAPPED_BAMS.out.merged_sorted_bam, outFolder)

        // Versions
        ch_versions = Channel.empty()
        ch_versions = ch_versions.mix(MAP_ONT_FASTQ.out.versions)
        ch_versions = ch_versions.mix(MERGE_MAPPED_BAMS.out.versions)
        ch_versions = ch_versions.mix(CHECKSUM_BAM.out.versions)
        ch_versions.unique().collectFile(name: 'map_merge_software_versions.yaml', storeDir: "${params.sampleDirectory}")

    emit:
        bam = MERGE_MAPPED_BAMS.out.merged_sorted_bam
        bai = MERGE_MAPPED_BAMS.out.bai
        bam_md5sum = CHECKSUM_BAM.out.md5sum
}