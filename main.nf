include { PACBIO_MAP_MERGE } from './workflows/pacbio-map-merge.nf'
include {
    ONT_BASECALL ;
    ONT_SETUP_BASECALL_ENVIRONMENT ;
    ONT_MERGE_QC_SUP_BAMS ;
    ONT_BACKUP_LIVEMODEL ;
    PUBLISH_RELEASE ;
    PUBLISH_RELEASE_QC as PUBLISH_RELEASE_QC_ROOT;
    PUBLISH_RELEASE_QC as PUBLISH_RELEASE_QC_NANAPLOT;
} from './workflows/ont-basecall-signals.nf'
include { ONT_MAP_MERGE_FASTQS } from './workflows/ont-map-merge-fastqs.nf'
include { ONT_MAP_MERGE_BAMS } from './workflows/ont-map-merge-bams.nf'
include { LONGREAD_QC } from './workflows/qc.nf'

workflow {
    NwgcCore.init(params)

    // Map-Merge
    if (params.mergedBam == null) {
        if (params.sequencingPlatform.equalsIgnoreCase("PacBio")) {
            PACBIO_MAP_MERGE()
            LONGREAD_QC(PACBIO_MAP_MERGE.out.bam, PACBIO_MAP_MERGE.out.bai, 
                params.sampleDirectory, params.sampleQCDirectory, params.qcToRun)
        }
        else if (params.sequencingPlatform.equalsIgnoreCase("ONT")) {
            if (params.containsKey('ontReleaseData') && params.ontReleaseData) {
                // --ontReleaseData true
                log.info("ONT workspace framework: release rebasecall data")
                ONT_BACKUP_LIVEMODEL()
                ONT_MERGE_QC_SUP_BAMS()
                PUBLISH_RELEASE(
                    ONT_MERGE_QC_SUP_BAMS.out.bam, ONT_MERGE_QC_SUP_BAMS.out.bai, 
                    ONT_MERGE_QC_SUP_BAMS.out.bam_md5sum, 
                    params.sampleDirectory)

                // where to write qc
                ontDataFolder = NwgcONTCore.getONTDataFolder(params)
                releaseLiveModelQCFolder = NwgcONTCore.getReleaseSupQCFolder(ontDataFolder)
                log.info("ontReleaseData: releaseLiveModelQCFolder = ${releaseLiveModelQCFolder}")

                log.info("ontReleaseData: params.qcToRun = ${params.qcToRun}")
                def fundamentalQCs = params.qcToRun
                //// FIXME: commented out for speed testing
                //// def fundamentalQCs = ['samtools_stats','quality','nanoplot','fingerprint']
                LONGREAD_QC(ONT_MERGE_QC_SUP_BAMS.out.bam, ONT_MERGE_QC_SUP_BAMS.out.bai, 
                    params.sampleDirectory, releaseLiveModelQCFolder, fundamentalQCs)
                PUBLISH_RELEASE_QC_ROOT(LONGREAD_QC.out.qcouts.flatten(), params.sampleQCDirectory)
                PUBLISH_RELEASE_QC_NANAPLOT(LONGREAD_QC.out.nanoplotqcouts.flatten(), "${params.sampleQCDirectory}/nanoPlot")
                
            } else if (params.containsKey('ontBaseCall') && params.ontBaseCall) {
                // --ontBaseCall true
                // perform the basecalling
                log.info("ONT workspace framework: perform re-basecall")
                ONT_BASECALL()
                // all QCs: might be unnecessary
                // def fundamentalQCs = ['All']
                // all QCs: skip contamination until we handle insufficient coverage
                //          at least perform the samstat [error rate] and coverage [mean cov>=25]
                //          quality can be used for better estimates
                //          nanoplot can be used for trouble-shooting
                //          fingerprint - for additional identity trouble-shooting
                def fundamentalQCs = ['samtools_stats','coverage','quality','nanoplot','fingerprint']
                LONGREAD_QC(ONT_BASECALL.out.bam, ONT_BASECALL.out.bai, 
                    params.sampleDirectory, params.sampleQCDirectory, fundamentalQCs)

            } else if (params.containsKey('ontSetupBasecall') && params.ontSetupBasecall) {
                // --ontSetupBasecall true
                // setup the basecalling environment only
                // and submit basecalling job handled via "--ontBaseCall true"
                log.info("ONT workspace framework: rebasecall workspace setup")
                ONT_SETUP_BASECALL_ENVIRONMENT()
            } else if (params.containsKey('ontRebasecall') && params.ontRebasecall) {
                // use the framework of ONT workspace
                // set up framework environment
                //   + run original HAC operations (but write to ONT workspace too)
                ontDataFolder = NwgcONTCore.getONTDataFolder(params)
                releaseLiveModelFolder = NwgcONTCore.getReleaseLiveModelFolder(ontDataFolder)
                // releaseLiveModelFolder = "${params.sampleDirectory}"
                log.info("ONT workspace framework: rebasecall and releaseLiveModelFolder = ${releaseLiveModelFolder}")

                ONT_SETUP_BASECALL_ENVIRONMENT()

                ONT_MAP_MERGE_BAMS([
                    "outFolder": "${releaseLiveModelFolder}"
                    , "outPrefix": "${params.sampleId}.${params.sequencingTarget}"
                    ])
                PUBLISH_RELEASE(
                    ONT_MAP_MERGE_BAMS.out.bam, ONT_MAP_MERGE_BAMS.out.bai, 
                    ONT_MAP_MERGE_BAMS.out.bam_md5sum, 
                    params.sampleDirectory)

                releaseLiveModelQCFolder = NwgcONTCore.getReleaseLiveModelQCFolder(ontDataFolder)
                log.info("releaseLiveModelQCFolder = ${releaseLiveModelQCFolder}")

                def fundamentalQCs = params.qcToRun
                ////FIXME: commented out for speed testing
                //// def fundamentalQCs = ['samtools_stats','quality','nanoplot','fingerprint']
                LONGREAD_QC(ONT_MAP_MERGE_BAMS.out.bam, ONT_MAP_MERGE_BAMS.out.bai, 
                    releaseLiveModelFolder, releaseLiveModelQCFolder, fundamentalQCs)
                PUBLISH_RELEASE_QC_ROOT(LONGREAD_QC.out.qcouts.flatten(), params.sampleQCDirectory)
                PUBLISH_RELEASE_QC_NANAPLOT(LONGREAD_QC.out.nanoplotqcouts.flatten(), "${params.sampleQCDirectory}/nanoPlot")
            } else {
                // equivalent: --no_ontReleaseData --no_ontBaseCall --no_ontSetupBasecall --no_ontRebasecall
                // this is the environment from <= v1.7.0
                // run as standalone; i.e. does not use the ONT workspace framework
                log.info("Standalone mode: writing results directly to ${params.sampleDirectory}")
                ONT_MAP_MERGE_BAMS([
                    "outFolder": "${params.sampleDirectory}"
                    , "outPrefix": "${params.sampleId}.${params.sequencingTarget}"
                    ])

                def fundamentalQCs = params.qcToRun
                LONGREAD_QC(ONT_MAP_MERGE_BAMS.out.bam, ONT_MAP_MERGE_BAMS.out.bai, 
                    params.sampleQCDirectory, "${params.sampleQCDirectory}/nanoPlot", fundamentalQCs)
            }
        }
        else {
            error "Error:  Unknown sequencingPlatform: ${params.sequencingPlatform}."
        }
    }
    else {
        LONGREAD_QC(params.mergedBam, "${params.mergedBam}.bai", params.sampleDirectory, params.sampleQCDirectory, params.qcToRun)
    }
}

workflow.onError {
    if (params.containsKey('rabbitHost') && params.rabbitHost!="") {
        NwgcCore.error(workflow, "$params.sampleId")
    }
}

workflow.onComplete {
    if (params.containsKey('rabbitHost') && params.rabbitHost!="") {
        NwgcCore.processComplete(workflow, "$params.sampleId")
    }
}
