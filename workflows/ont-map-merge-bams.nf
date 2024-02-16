include { MAP_ONT_BAM } from '../modules/map_ont_bam.nf'
include { MERGE_MAPPED_BAMS } from '../modules/merge_mapped_bams.nf'
include { CHECKSUM_BAM } from '../modules/checksum_bam.nf'

workflow ONT_MAP_MERGE_BAMS {
    take:
        arguments
    main:
        ontBams = Channel.empty()
        for (ontBamFolder in params.ontBamFolders) {
            bamFolderChannel = Channel.fromPath(ontBamFolder + "/*.bam")
            ontBams = ontBams.concat(bamFolderChannel)
        }
        MAP_ONT_BAM(ontBams)

        ////def outFolder = "${params.sampleDirectory}"
        ////def outPrefix = "${params.sampleId}.${params.sequencingTarget}"
        def outFolder = "${arguments['outFolder']}"
        def outPrefix = "${arguments['outPrefix']}"

        // Merge
        MERGE_MAPPED_BAMS(MAP_ONT_BAM.out.mapped_bam.collect(), outFolder, outPrefix)

        // checksum
        CHECKSUM_BAM(MERGE_MAPPED_BAMS.out.merged_sorted_bam, outFolder)

        // Versions
        ch_versions = Channel.empty()
        ch_versions = ch_versions.mix(MAP_ONT_BAM.out.versions)
        ch_versions = ch_versions.mix(MERGE_MAPPED_BAMS.out.versions)
        ch_versions = ch_versions.mix(CHECKSUM_BAM.out.versions)
        ch_versions.unique().collectFile(name: 'map_merge_software_versions.yaml', storeDir: "${params.sampleDirectory}")

    emit:
        bam = MERGE_MAPPED_BAMS.out.merged_sorted_bam
        bai = MERGE_MAPPED_BAMS.out.bai
        bam_md5sum = CHECKSUM_BAM.out.md5sum
}
