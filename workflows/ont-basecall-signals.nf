import java.util.regex.Pattern
import java.util.regex.Matcher
include { MAP_ONT_BAM } from '../modules/map_ont_bam.nf'
include { MERGE_MAPPED_BAMS } from '../modules/merge_mapped_bams.nf'
include { CHECKSUM_BAM } from '../modules/checksum_bam.nf'

// DONE
process CONVERT_FAST5_TO_POD5 {
    label "CONVERT_FAST5_TO_POD5_${params.sampleId}_${params.userId}"

    // keep these pod5 for future use
    publishDir "${signalParentFolder}", mode: 'link', pattern: 'fast5_to_pod5/*.pod5'

    input:
        tuple val(chunk_idx), path('*')
        val(signalParentFolder)
    output:
        tuple val(chunk_idx), path("fast5_to_pod5/*.pod5")
    script:
        def eff_signal_path = "fast5_to_pod5/"
        """
            pod5 convert fast5 ./*.fast5 --output ${eff_signal_path} --threads 1 --one-to-one ./
        """
}

// DONE
process BASECALL_ONLY_DORADO {
    label "BASECALL_ONLY_DORADO_${params.sampleId}_${params.userId}"

    input:
        tuple val(chunk_idx), path('*')
        val(basecaller_model)
        val(basecaller_mods)
        val(signalExtensions)
    output:
        path("${params.sampleId}_${chunk_idx}.bam"), emit: uabams
        path("*.pod5"), emit: converted_pod5s, optional: true
    script:
        def cudaDevice = (null == params.ontCudaDevice) ? 'cuda:all' : params.ontCudaDevice
        // TODO: decide if R9 or R10?
        def model_arg = (null == basecaller_model) ? "\${MODELS_ROOT}/dna_r10.4.1_e8.2_400bps_sup@v3.5.2" : "${basecaller_model}"
        def mods_arg = (null == basecaller_mods) ? '' : "--modified-bases ${basecaller_mods}"

        def eff_signal_path = "."
        """
            dorado basecaller \
                --device ${cudaDevice} \
                ${model_arg} \
                ${eff_signal_path} \
                ${mods_arg} \
            > ${params.sampleId}_${chunk_idx}.bam
        """
}

// DONE
process SUMMARIZE_DORADO {
    label "SUMMARIZE_DORADO_${params.sampleId}_${params.userId}"

    input:
        path bam
    output:
        path("${bam.baseName}.summary.tsv.gz"), emit: summary
    script:
    """
        dorado summary ${bam} \
        | gzip -c \
        > ${bam.baseName}.summary.tsv.gz
    """
}

// FIXME: parameterize publishDir!
process COLLATE_DORADO_SUMMARIES {
    label "COLLATE_DORADO_SUMMARIES_${params.sampleId}_${params.userId}"

    publishDir "${params.ontBaseCallOutFolder}", mode: 'link', pattern: '*.summary.tsv.gz'

    input:
        path tsvs
    output:
        path("${params.sampleId}.${params.sequencingTarget}.summary.tsv.gz"), emit: summary
    shell:
    '''
        n=0
        for fname in !{tsvs}; do
            if [ $n -eq 0 ]; then
                n=1
                zcat $fname
            else
                zcat $fname | tail -n +2
            fi
        done \
        | gzip -c \
        > !{params.sampleId}.!{params.sequencingTarget}.summary.tsv.gz
    '''
}

// DONE
process QSFILTER {
    label "QSFILTER_${params.sampleId}_${params.userId}"

    publishDir "${pass_folder}", mode: 'link', pattern: '*.pass.bam'
    publishDir "${fail_folder}", mode: 'link', pattern: '*.fail.bam'
    
    input:
        path reads
        val(qscore_filter)
        val(pass_folder)
        val(fail_folder)
    output:
        path("${reads.baseName}.pass.bam"), emit: pass
        path("${reads.baseName}.fail.bam"), emit: fail
    script:
    """
        samtools view \
            -@ ${task.cpus} \
            -e '[qs] >= ${qscore_filter}' \
            ${reads} \
            -o ${reads.baseName}.pass.bam \
            -U ${reads.baseName}.fail.bam
    """
}

// DONE
workflow ONT_BASECALL {
    main:
        def finalMergedPath = "${params.ontBaseCallOutFolder}/mapped_pass_sup"
        String runAcqID = new File("${params.ontBaseCallOutFolder}").name
        def finalMergedFilePrefix = "${params.sampleId}.${runAcqID}"
        File finalMergedFileBam = new File("${finalMergedPath}/${finalMergedFilePrefix}.bam")
        File finalMergedFileBai = new File("${finalMergedPath}/${finalMergedFilePrefix}.bam.bai")
        File finalMergedFileMD5 = new File("${finalMergedPath}/${finalMergedFilePrefix}.bam.md5sum")
        ////FIXME: by pass checking
        if (finalMergedFileBam.exists() && finalMergedFileBai.exists() && finalMergedFileMD5.exists()) {
        ////if (finalMergedFileBam.exists() && finalMergedFileBai.exists() && finalMergedFileMD5.exists() && false) {
            log.info("Skipping re-basecalling for ${params.ontBaseCallOutFolder}")
            log.info("${finalMergedFileMD5.toString()} already exists")
            log.info("Assuming that ${finalMergedFileBam.toString()} is current and intact")
            log.info("Assuming that ${finalMergedFileBai.toString()} is current and intact")
            exit 0
        } else {

            // our default naming is: PAM05070_pass_c2f39984_a1f291a3_0
            def pattern_FC_PF_runid_acqid_chunkid = /[A-Z]{3}[0-9]{5}_(pass|fail)_[0-9a-zA-Z]{8}_[A-Za-z0-9]{8}_[0-9]*/

            // TODO: null
            // ${MODELS_ROOT}
            basecaller_model = file(params.ontBaseCallModel, type: "dir")
            log.info "Basecaller model = ${basecaller_model}"

            basecaller_mods = params.ontBaseCallBaseModifications
            log.info "Basecaller modifications = ${basecaller_mods}"

            Integer chunk_idx = 0
            ontSignals_chunks = Channel.empty()
            for (ontSignalFolder in params.ontSignalFolders) {
                signalFolderChannel = Channel.fromPath(ontSignalFolder + "/*.${params.signalExtensions}", checkIfExists: true)
                ontSignals_chunks = ontSignals_chunks.concat(signalFolderChannel)
            }
            ontSignals_chunks
                .toSortedList()
                .flatten()
                .branch{
                    pattern_FC_PF_runid_acqid_chunkid: it.baseName ==~ pattern_FC_PF_runid_acqid_chunkid
                    no_pattern: true
                }.set{signal_name_types}
            signal_name_types.no_pattern.first().subscribe { 
                log.warn "Unexpected signal file naming convention detected. ${it.baseName}"
            }

            // Split the name keeping the flow cell, run, pass/fail and pod5 index.
            signal_name_types.pattern_FC_PF_runid_acqid_chunkid
                .map{filename -> 
                    fields = filename.baseName.split("_")
                    [fields[-1] as int, fields[0], fields[3], fields[1], filename]
                }
                .map{ pod5_index, cell_id, run_id, pass, pod5 ->
                    [Math.floor(pod5_index/params.basecaller_num_chunk_per_batch), cell_id, run_id, pass, pod5]
                }
                .groupTuple(
                    by: [0, 1, 2, 3]
                )
                .map{pod5_index, cell_id, run_id, pass, pod5s -> pod5s}
                .mix(
                    signal_name_types.no_pattern
                        .buffer(size:params.basecaller_num_chunk_per_batch, remainder: true)
                )
                .map { pod5s -> [chunk_idx++, pod5s] }
                .set{ontSignalsBatches}

            if ((null == basecaller_model)) {
                // TODO: decide if R9 or R10?
                basecaller_model = "\${MODELS_ROOT}/dna_r10.4.1_e8.2_400bps_sup@v3.5.2"
                log.warn("No basecaller model specified. Defaul to ${basecaller_model}")
            } else {
                if (!basecaller_model.exists()) {
                    log.warn("Specified basecaller model ${basecaller_model} not found.")
                    basecaller_model = "\${MODELS_ROOT}/${basecaller_model.name}"
                    log.warn("Override to ${basecaller_model}")
                    // FIXME: informative but delay as it requires module loading
                    // log.warn("\${MODELS_ROOT} = ${MODELS_ROOT}")
                } else {
                    log.info("Using specified basecaller model ${basecaller_model}")
                }
            }

            if (params.signalExtensions == "fast5") {
                pod5conversion = CONVERT_FAST5_TO_POD5(ontSignalsBatches, "${params.ontBaseCallOutFolder}")
                called_bams = BASECALL_ONLY_DORADO(
                    pod5conversion, 
                    basecaller_model,
                    basecaller_mods,
                    params.signalExtensions)
            } else {
                called_bams = BASECALL_ONLY_DORADO(
                    ontSignalsBatches,
                    basecaller_model,
                    basecaller_mods,
                    params.signalExtensions
                )
            }

            // Compute summary
            SUMMARIZE_DORADO(called_bams.uabams) | collect | COLLATE_DORADO_SUMMARIES    
            summary = COLLATE_DORADO_SUMMARIES.out.summary

            // collate pass and fail
            collatedPassFailBAMs = QSFILTER(
                called_bams.uabams, 
                params.basecall_qscore_filter, 
                "${params.ontBaseCallOutFolder}/bam_pass_sup",
                "${params.ontBaseCallOutFolder}/bam_fail_sup")


            // FIXME: collatedPassFailBAMs.fail can be useful for SVs detection!
            //        we skip collatedPassFailBAMs.fail to save time
            // perform the mapping on the pass
            MAP_ONT_BAM(collatedPassFailBAMs.pass)

            // Merge
            MERGE_MAPPED_BAMS(
                MAP_ONT_BAM.out.mapped_bam.collect(), 
                finalMergedPath, 
                finalMergedFilePrefix)

            // checksum
            CHECKSUM_BAM(MERGE_MAPPED_BAMS.out.merged_sorted_bam, finalMergedPath)

            // Versions
            ch_versions = Channel.empty()
            ch_versions = ch_versions.mix(MAP_ONT_BAM.out.versions)
            ch_versions = ch_versions.mix(MERGE_MAPPED_BAMS.out.versions)
            ch_versions = ch_versions.mix(CHECKSUM_BAM.out.versions)
            ch_versions.unique().collectFile(name: 'map_merge_software_versions.yaml', storeDir: "${params.ontBaseCallOutFolder}")
        }

    emit:
        // chunked_pass_bams = collatedPassFailBAMs.pass
        // chunked_fail_bams = collatedPassFailBAMs.fail
        // summary = summary
        // converted_pod5s = called_bams.converted_pod5s
        // from second phase
        bam = MERGE_MAPPED_BAMS.out.merged_sorted_bam
        bai = MERGE_MAPPED_BAMS.out.bai
        bam_md5sum = CHECKSUM_BAM.out.md5sum
}

// DONE
workflow ONT_SETUP_BASECALL_ENVIRONMENT {
    main:
        NwgcONTCore.setLog(log)
        NwgcONTCore.setWorkflowMetadata(workflow)
        ontDataFolder = NwgcONTCore.getONTDataFolder(params)
        outPrefix = NwgcONTCore.getReleaseSupPrefix(params)
        ontSubmitBaseCallJob = true
        if (params.containsKey('ontSubmitBaseCallJob')) {
            ontSubmitBaseCallJob = params.ontSubmitBaseCallJob
        }
        def ontRunAcqs = Channel
            .fromList(params.ontBamFolders)
            .map{runAcqBam -> 
                File path = new File(runAcqBam)
                [path.parent, runAcqBam]
            }
            .groupTuple(
                by: [0]
            )
            .map{runAcqFolder, runAcqBams -> 
                NwgcONTCore.setupRunAcquisition(runAcqFolder, runAcqBams, params.sampleId, ontDataFolder, ontSubmitBaseCallJob)
            }

        /*
        // inform user how-to submit
        ontRunAcqs.subscribe {
            println ""
            println "##########"
            println "### To basecall for ${params.sampleId} - ${it.runAcqID}"
            println "##########"
            println "cd ${it.runacq.workdir}"
            def scriptFPN = new File("${it.basecall.script}")
            println "qsub -P dna -q il.q ${scriptFPN.name}"
            println "# submitted jobid: ${it.basecall.jobid}"
            println ""
        }
        */

        // TODO: consider => write the release version of the file?
        ontRunAcqs
            .toList()
            .subscribe {
                NwgcONTCore.setupRelease(it, params.sampleId, ontDataFolder, outPrefix)
            }
}

// DONE
process PUBLISH_RELEASE {
    label "PUBLISH_RELEASE_${params.sampleId}_${params.userId}"

    publishDir "${pubdir}", mode: 'copy'
    input:
        file bam
        file bai
        file md5sum
        val pubdir
    output:
        file bam
        file bai
        file md5sum
    script:
    """
    echo ${bam}
    echo ${bai}
    echo ${md5sum}
    echo ${pubdir}
    """
}

// DONE
process PUBLISH_RELEASE_QC {
    label "PUBLISH_RELEASE_QC_${params.sampleId}_${params.userId}"

    publishDir "${outdir}", mode: 'copy'
    input:
        path qcouts
        val outdir
    output:
        path qcouts
    script:
    """
    echo ${qcouts}
    echo ${outdir}
    """
}


workflow ONT_BACKUP_LIVEMODEL {
    main:
        NwgcONTCore.setLog(log)
        NwgcONTCore.setWorkflowMetadata(workflow)
        NwgcONTCore.handleLiveModelBackup(params)
        // TODO: latch the next step?!
}

// DONE
workflow ONT_MERGE_QC_SUP_BAMS {
    main:
        NwgcONTCore.setLog(log)
        NwgcONTCore.setWorkflowMetadata(workflow)
        // rehash the release bam files from SUP specific directories
        ontDataFolder = NwgcONTCore.getONTDataFolder(params)
        outFolder = NwgcONTCore.getReleaseSupFolder(ontDataFolder)
        outPrefix = NwgcONTCore.getReleaseSupPrefix(params)

        def ontBams = Channel.empty()
        def bamFile = null
        for (ontBamFolder in params.ontBamFolders) {
            // NOTE: we default to PASS reads ONLY!
            bamFile = NwgcONTCore.getUnchunkedMappedBamFileFromBamFolder(ontBamFolder, params.sampleId, ontDataFolder, true)
            ontBams = ontBams.concat(Channel.of(bamFile))
        }

        ontBams
            .flatten()
            .unique()
            .set{ontUniqueBams}

        ontUniqueBams
            .branch{
                missing: !(new File("${it}.md5sum").exists())
                present: true
            }.set{ontBamsBinned}

        // inform user of all missing bam file(s)
        int numMisses = 0
        int numOutdateds = 0
        ontBamsBinned
            .missing
            .subscribe onNext: { 
                    if (new File("${it}").exists()) {
                        log.warn("Outdated bam: ${it}")
                        numOutdateds++;
                    } else {
                        log.warn("Missing bam: ${it}")
                        numMisses++;
                    }
                }, onComplete: { 
                    if (numMisses>0 || numOutdateds>0) {
                        log.error("${numOutdateds} outdated bam file(s).")
                        log.error("${numMisses} missing bam file(s).")
                        throw new RuntimeException("${numMisses} missing and ${numOutdateds} outdated bam file(s)! Please perform basecalling with SUP model.")
                    }
                }

        // Merge
        MERGE_MAPPED_BAMS(ontBamsBinned.present.flatten().unique().collect(), outFolder, outPrefix)

        // checksum
        CHECKSUM_BAM(MERGE_MAPPED_BAMS.out.merged_sorted_bam, outFolder)

        // TODO: copy? hardlink? 

        // Versions
        ch_versions = Channel.empty()
        // ch_versions = ch_versions.mix(MAP_ONT_BAM.out.versions)
        ch_versions = ch_versions.mix(MERGE_MAPPED_BAMS.out.versions)
        ch_versions = ch_versions.mix(CHECKSUM_BAM.out.versions)
        ch_versions.unique().collectFile(name: 'merge_sup_software_versions.yaml', storeDir: "${params.sampleDirectory}")

    emit:
        bam = MERGE_MAPPED_BAMS.out.merged_sorted_bam
        bai = MERGE_MAPPED_BAMS.out.bai
        bam_md5sum = CHECKSUM_BAM.out.md5sum
        outFolder = outFolder
}
