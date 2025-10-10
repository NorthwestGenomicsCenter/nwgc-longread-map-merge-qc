include { STRIP_KINETICS } from '../modules/strip_kinetics.nf'
include { MAP_HIFI_BAM } from '../modules/map_hifi_bam.nf'
include { MAP_PACBIO_FASTQ } from '../modules/map_pacbio_fastq.nf'
include { MERGE_MAPPED_BAMS } from '../modules/merge_mapped_bams.nf'
include { ADD_NM_TAGS } from '../modules/add_nm_tags.nf'
include { ADD_NM_TAGS as ADD_NM_TAGS_FASTQS } from '../modules/add_nm_tags.nf'
include { CHECKSUM_BAM } from '../modules/checksum_bam.nf'

workflow PACBIO_MAP_MERGE {

    main:
        ch_hiFiBams = Channel.empty()
        if (params.hiFiBams) {
            ch_hiFiBams = Channel.fromPath(hiFiBams)
        }

        ch_fastqs = Channel.empty()
        if (params.pacbioFastqsFolder) {
           fastqFolderChannel = Channel.fromPath( params.pacbioFastqsFolder + "/*.fastq")
           ch_fastqs = ch_fastqs.concat(fastqFolderChannel)
           zippedFastqFolderChannel = Channel.fromPath( params.pacbioFastqsFolder + "/*.fastq.gz")
           ch_fastqs = ch_fastqs.concat(zippedFastqFolderChannel)
        }

        // Map HiFi Bams
        if (params.stripKinetics) {
            STRIP_KINETICS(ch_hiFiBams)
            MAP_HIFI_BAM(STRIP_KINETICS.out.bam)
        }
        else {
            MAP_HIFI_BAM(ch_hiFiBams)
        }

        // Map PacBio Fastqs
        MAP_PACBIO_FASTQ(ch_fastqs)

        // NM TAGS
        ADD_NM_TAGS(MAP_HIFI_BAM.out.mapped_bam)
        ADD_NM_TAGS_FASTQS(MAP_PACBIO_FASTQ.out.mapped_bam)

        def outFolder = "${params.sampleDirectory}"
        def outPrefix = "${params.sampleId}.${params.sequencingTarget}"
        // Merge
        ch_mapOut = Channel.empty()
        ch_mapOut = ch_mapOut.mix(ADD_NM_TAGS.out.nm_bam.collect())
        ch_mapOut = ch_mapOut.mix(ADD_NM_TAGS_FASTQS.out.nm_bam.collect())
        MERGE_MAPPED_BAMS(ch_mapOut, outFolder, outPrefix)

        // checksum
        CHECKSUM_BAM(MERGE_MAPPED_BAMS.out.merged_sorted_bam, MERGE_MAPPED_BAMS.out.bai, outFolder)

        // Versions
        ch_versions = Channel.empty()
        ch_versions = ch_versions.mix(MAP_HIFI_BAM.out.versions)
        ch_versions = ch_versions.mix(ADD_NM_TAGS.out.versions)
        ch_versions = ch_versions.mix(MERGE_MAPPED_BAMS.out.versions)
        ch_versions = ch_versions.mix(CHECKSUM_BAM.out.versions)
        ch_versions.unique().collectFile(name: 'map_merge_software_versions.yaml', storeDir: "${params.sampleDirectory}")

    emit:
        bam = MERGE_MAPPED_BAMS.out.merged_sorted_bam
        bai = MERGE_MAPPED_BAMS.out.bai
        bam_md5sum = CHECKSUM_BAM.out.bammd5sum
        bai_md5sum = CHECKSUM_BAM.out.baimd5sum
}
