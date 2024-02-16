include { SAMTOOLS_STATS } from '../modules/samtools_stats.nf'
include { PICARD_QUALITY_METRICS } from '../modules/picard_quality_metrics.nf'
include { PICARD_COVERAGE_METRICS } from '../modules/picard_coverage_metrics.nf'
include { NANO_PLOT } from '../modules/nano_plot.nf'
include { CONTAMINATION_CHECK } from '../modules/contamination_check.nf'
include { CREATE_FINGERPRINT_VCF } from '../modules/create_fingerprint_vcf.nf'

ch_versions = Channel.empty()

workflow LONGREAD_QC {

    // Input is a merged bam (and corresponding bai) that has NM tags
    take:
        bam
        bai
        sampleFolder
        qcFolder
        qcToRun

    main:

        def runAll = qcToRun.contains("All")
        def ch_qcouts = Channel.empty()
        def ch_nanoplotqcouts = Channel.empty()

        log.info("LONGREAD_QC: qcFolder = ${qcFolder}")

        if (runAll || qcToRun.contains("samtools_stats")) {
            SAMTOOLS_STATS(bam, bai, qcFolder)
            ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
            ch_qcouts = ch_qcouts.mix(SAMTOOLS_STATS.out.stats)
        }

        if (runAll || qcToRun.contains("quality")) {
            PICARD_QUALITY_METRICS(bam, bai, qcFolder)
            ch_versions = ch_versions.mix(PICARD_QUALITY_METRICS.out.versions)
            ch_qcouts = ch_qcouts.mix(PICARD_QUALITY_METRICS.out.stats)
        }

        if (runAll || qcToRun.contains("coverage")) {
            PICARD_COVERAGE_METRICS(bam, bai, qcFolder)
            ch_versions = ch_versions.mix(PICARD_COVERAGE_METRICS.out.versions)
            ch_qcouts = ch_qcouts.mix(PICARD_COVERAGE_METRICS.out.stats)
        }

        if (runAll || qcToRun.contains("nanoplot")) {
            NANO_PLOT(bam, bai, qcFolder)
            ch_versions = ch_versions.mix(NANO_PLOT.out.versions)
            ch_nanoplotqcouts = ch_nanoplotqcouts.mix(NANO_PLOT.out.stats, NANO_PLOT.out.html, NANO_PLOT.out.png)
        }

        if (runAll || qcToRun.contains("contamination")) {
            CONTAMINATION_CHECK(bam, bai, qcFolder)
            ch_versions = ch_versions.mix(CONTAMINATION_CHECK.out.versions)
            ch_qcouts = ch_qcouts.mix(CONTAMINATION_CHECK.out.verifications)
        }

        if (runAll || qcToRun.contains("fingerprint")) {
            CREATE_FINGERPRINT_VCF(bam, bai, qcFolder)
            ch_versions = ch_versions.mix(CREATE_FINGERPRINT_VCF.out.versions)
            ch_qcouts = ch_qcouts.mix(CREATE_FINGERPRINT_VCF.out.vcfs, CREATE_FINGERPRINT_VCF.out.vcfindexes)
        }

        ch_versions.unique().collectFile(name: 'qc_software_versions.yaml', storeDir: "${sampleFolder}")

    emit:
        qcouts = ch_qcouts
        nanoplotqcouts = ch_nanoplotqcouts
}