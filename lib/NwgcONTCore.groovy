//
// This file holds several functions specific to ONT SUP data processing framework
//

import nextflow.Nextflow
import java.nio.file.Files
import org.yaml.snakeyaml.Yaml
import java.nio.file.NoSuchFileException
import java.util.regex.Matcher;
import java.util.regex.Pattern;

class NwgcONTCore {

    private static log = null
    private static wfMeta = null
    final private static String MAPPED_PASS_SUP = "mapped_pass_sup"
    final private static String MAPPED_FAIL_SUP = "mapped_fail_sup"
    final private static String PARENT_BASH_SCRIPT = "nwgc-longread-map-merge-qc.sh"
    final private static String BASECALL_BASH_SCRIPT = "nwgc-ont-basecall-map-merge.sh"
    final private static String RELEASE_BASH_SCRIPT = "nwgc-ont-merge-qc.sh"
    final private static String SIGNAL_FAST5_EXT = "fast5"
    final private static String SIGNAL_POD5_EXT = "pod5"
    final private static String DORADO_BASE_MODIFICATIONS = '5mCG'
    final private static String BASECALL_R9_SUP_MODEL = 'dna_r9.4.1_e8_sup@v3.3'
    final private static String BASECALL_R10_SUP_MODEL = 'dna_r10.4.1_e8.2_400bps_sup@v3.5.2'
    final private static String RELEASE_SUP_FOLDER = 'merge-qc_sup'
    final private static String RELEASE_LIVEMODEL_FOLDER = 'map-merge-qc_livemodel'
    final private static String ONT_BACKEND_USERID = 'cheehong'
    final private static String ONT_BACKEND_USEREMAIL = 'cheehong@uw.edu'
    final private static String MODULE_NEXTFLOW_VERSION = 'nextflow/23.10.1'
    final private static Map FlowCellModels = [
        'FLO-PRO002@4000': ['supmodel': 'dna_r9.4.1_e8_sup@v3.3', 'modifications': '5mCG'], 
        'FLO-PRO114M@4000': ['supmodel': 'dna_r10.4.1_e8.2_400bps_sup@v3.5.2', 'modifications': '5mCG'],
        'FLO-PRO114M@5000': ['supmodel': 'dna_r10.4.1_e8.2_400bps_sup@v4.2.0', 'modifications': '5mCG_5hmCG'],
        '': ['supmodel': 'dna_r10.4.1_e8.2_400bps_sup@v3.5.2', 'modifications': '5mCG']
        ]
    // use internal compute subnet reference
    final private static String BASECALL_MQ_HOST = 'burnet'
    final private static int BASECALL_MQ_PORT = 5671
    final private static String BASECALL_MQ_VHOST = 'nextflow'
    final private static String BASECALL_MQ_EXCHANGE = ''
    final private static String BASECALL_MQ_ROUTE_KEY = 'OntSupRebaseCalling'

    NwgcONTCore() {}

    public static setLog(extlog) { log = extlog }
    public static setWorkflowMetadata(existingWFMeta) { wfMeta = existingWFMeta }

    public static String getONTDataFolder(params) {
        if (params.ontDataFolder) {
            return params.ontDataFolder
        } else {
            String parent = new File(params.sampleDirectory).parent
            return new File("${parent}/ont").toString()
        }
    }

    private static getONTFolder(runAcqFolder, subFolder) {
        if (subFolder.equals("")) {
            return new File(runAcqFolder).toString()
        } else {
            return new File("${runAcqFolder}/${subFolder}").toString()
        }
    }

    private static getONTInstrumentFolder(runAcqFolder, subfolder) {
        return getONTFolder(runAcqFolder, subfolder)
    }

    private static getONTWorkspaceFolder(parentDir, runAcqId, subfolder) {
        return getONTFolder("${parentDir}/${runAcqId}", subfolder)
    }

    // extract the param-file specified in the nextflow run
    private static getParamFile(launchDir, commandLine) {
        // FIXME: handle spaces in parameters value!
        String paramFile = ""
        String[] commandTokens = commandLine.split("\\s+");
        for(int i=0; i<commandTokens.length ; i++) {
            if (commandTokens[i].equals("-params-file")) {
                if ((i+1)<commandTokens.length) {
                    File templateFPN = new File(commandTokens[i+1])
                    paramFile = templateFPN.getName()
                }
                break;
            }
        }
        return "${launchDir}/${paramFile}"
    }

    // hardlink instrument data folder to the framework working folder
    private static Map setupDataInFolder(source, destination) {
        // foreach files in source, create hardlink in destination if non-existent
        // FIXME: does not take care of subfolder! (I.e. demux barcoded samples are not handled!!)
        if (log!=null) { log.info("Setup data subfolder ${source} --> ${destination}") }

        // println "${item} --> ${destFile}"
        File destFolder = new File(destination)
        if (!destFolder.exists()) {
            destFolder.mkdirs()
            if (log!=null) { log.info("Setup data subfolder ${destination} created.") }
        }

        int totalFiles = 0
        int totalSetup = 0
        int totalSkipped = 0
        File sourceFile = new File(source)
        sourceFile.eachFile {
            item ->
            if (item.isFile()) {
                File destFile = new File("${destination}/${item.getName()}")
                if (!destFile.exists()) {
                    try {
                        // hardlink cannot work across device!
                        Files.createLink(destFile.toPath(), item.toPath())
                        totalSetup++
                    } catch (IOException e) {
                        if (e.message.contains('Invalid cross-device link')) {
                            if (log!=null) { log.warn("Attempt to softlink after hardlink error for ${item.name}") }
                            Files.createSymbolicLink(destFile.toPath(), item.toPath())
                            totalSetup++
                        } else {
                            log.error("${item.toString()} ioexception ${e.message}")
                        }
                    }
                }
                totalFiles++
            } else if (item.isDirectory()) {
                if (log!=null) { log.info("Setup data subfolder ${item.getName()} - folder skipped") }
                totalSkipped++;
            }
        }
        if (log!=null) { log.info("Setup data subfolder ${totalFiles} data file(s), ${totalSetup} created link(s), ${totalSkipped} folder(s) skipped.") }

        return ['total': totalFiles, 'setup': totalSetup, 'skipped': totalSkipped, 'source': source, 'destination': destination]
    }

    private static String getRunAcquisitionMappedBamFile(ontDataFolder, runAcqID, sampleId, passed=true) {
        String mapp_sup_fpn = getONTWorkspaceFolder(ontDataFolder, runAcqID, passed ? MAPPED_PASS_SUP : MAPPED_FAIL_SUP)
        return new File("${mapp_sup_fpn}/${sampleId}.${runAcqID}.bam").toString()
    }

    // get the framework ("cached") single mapped bam for a single run-acquisition
    // (run-acquisition data folder contains chunked unmapped bam files)
    public static String getUnchunkedMappedBamFileFromBamFolder(ontBamFolder, sampleId, ontDataFolder, passed = true) {
        String parent = new File(ontBamFolder).parent
        String runAcqID = new File(parent).name
        return getRunAcquisitionMappedBamFile(ontDataFolder, runAcqID, sampleId, passed)
    }

    // write a new bash job wrapper for the basecalling task
    // (base off the trigger from merge action generated bash job wrapper)
    private static setupBasecallBashScript(template, derivative, sampleId, ontDataFolder, runAcqID, totalSetup, yamlFile, runMeta, backendEmailReroute=false) {
        File templateFile = new File(template)
        assert templateFile.exists() : "Source script missing! ${template}"
        File derivativeFile = new File(derivative)
        def writer = derivativeFile.newWriter()

        boolean notesWritten = false
        boolean nxfuuidWritten = false
        boolean workingDirWritten = false
        templateFile.eachLine {
            line -> 
            if (line.startsWith("#!")) {
                writer.writeLine "${line}"
            } else if (line.startsWith("#\$ ")) {
                if (line.startsWith("#\$ -N ")) { // update batch job name
                    // writer.writeLine "${line}_${runAcqID}"
                    writer.writeLine "#\$ -N basecall_${sampleId}_${runAcqID}"
                } else if (line.startsWith("#\$ -o ")) { // update stdout folder
                    writer.writeLine "#\$ -o ${ontDataFolder}/${runAcqID}/logs/"
                } else if (line.startsWith("#\$ -e ")) { // update stderr folder
                    writer.writeLine "#\$ -e ${ontDataFolder}/${runAcqID}/logs/"
                } else if (line.startsWith("#\$ -M ") && backendEmailReroute) { // update notification email
                    writer.writeLine "#\$ -M ${ONT_BACKEND_USEREMAIL}"
                } else { // pass-thru
                    writer.writeLine "${line}"
                }
            } else {
                if (!notesWritten) {
                    // inform user : autogen; will be overwritten
                    writer.writeLine ""
                    writer.writeLine "##########################################"
                    writer.writeLine "# IMPORTANT: auto-generated; any changes WILL BE OVERWRITTEN by pipeline with Samplify's merge action"
                    writer.writeLine "#            source="
                    writer.writeLine "#            ${template}"
                    writer.writeLine "#"
                    writer.writeLine "# Bash script for SUP basecalling on a single run-acquisition signal data"
                    writer.writeLine "#"
                    writer.writeLine "##########################################"
                    writer.writeLine ""

                    notesWritten = true
                }
                if (line.startsWith("PUBLISH_DIR=")) { // move to run-acquisition folder
                    writer.writeLine "PUBLISH_DIR=${ontDataFolder}/${runAcqID}"
                } else if (line.startsWith("WORKING_DIR=")) { // for nextflow framework
                    writer.writeLine "${line}_${runAcqID}"
                } else if (line.startsWith("PIPELINE_YAML=")) { // for nextflow framework
                    writer.writeLine "PIPELINE_YAML=${yamlFile}"
                } else if (line.startsWith("JOB_NAME=")) { // for nextflow framework
                    //writer.writeLine "${line}_${runAcqID}_${totalSetup}"
                    writer.writeLine "JOB_NAME=basecall_${sampleId}_${runAcqID}"
                } else if (line.startsWith("USER_EMAIL=") && backendEmailReroute) { // for nextflow framework
                    writer.writeLine "USER_EMAIL=${ONT_BACKEND_USEREMAIL}"
                } else if (line.equals("module load nextflow/22.10.7")) { // module load nextflow/22.10.7 - unavailable!
                    writer.writeLine "module load ${MODULE_NEXTFLOW_VERSION}"
                } else if (line.startsWith("export NXF_UUID=")) { // should have a new session id for samplify database
                    def uuid = UUID.randomUUID()
                    writer.writeLine "export NXF_UUID=${uuid}"
                    nxfuuidWritten = true
                } else if (line.startsWith('cd ${PUBLISH_DIR}')) { // cd ${PUBLISH_DIR}
                    // use the working directory per launchify
                    writer.writeLine 'cd ${WORKING_DIR}'
                    workingDirWritten = true
                } else if (line.startsWith("nextflow ")) { // nextflow run ...
                    if (!nxfuuidWritten) {
                        def uuid = UUID.randomUUID()
                        writer.writeLine "export NXF_UUID=${uuid}"
                        nxfuuidWritten = true
                    }
                    if (!workingDirWritten) {
                        writer.writeLine 'mkdir -p ${WORKING_DIR}'
                        writer.writeLine 'cd ${WORKING_DIR}'
                        workingDirWritten = true
                    }
                    writer.writeLine "${line}"
                } else { // pass-thru
                    writer.writeLine "${line}"
                }
            }
        }
        // publish to rabbitMQ
        writer.writeLine "exitcode=\$?"
        writer.writeLine "# publish to message queue"
        writer.writeLine "amqps_publish ${BASECALL_MQ_HOST} ${BASECALL_MQ_PORT} ${BASECALL_MQ_VHOST} '${BASECALL_MQ_EXCHANGE}' ${BASECALL_MQ_ROUTE_KEY} \"sampleId: \\\"${sampleId}\\\""
        writer.writeLine "protocolGroupId: \\\"${runMeta.protocolGroupId}\\\""
        writer.writeLine "flowCellPosition: \\\"${runMeta.flowCellPosition}\\\""
        writer.writeLine "flowCellId: \\\"${runMeta.flowCellId}\\\""
        writer.writeLine "protocolRunIdShort: \\\"${runMeta.protocolRunIdShort}\\\""
        writer.writeLine "supBam: \\\"${runMeta.supBam}\\\""
        writer.writeLine "exitCode: \${exitcode}\""
        writer.flush()
        writer.close()

        if (log!=null) { log.info("Setup basecall bash script ${template} --> ${derivativeFile}") }
    }

    // write a new YAML parameters for the basecalling task
    // (base off the trigger from merge action generated YAML parameters)
    private static setupSUPBasecallParamsYAML(template, derivative, ontDataFolder, runAcqID, settings, backendEmailReroute=false) {

        File templateFile = new File(template)
        assert templateFile.exists() : "Parameters file missing! ${template}"
        File derivativeFile = new File(derivative)
        def writer = derivativeFile.newWriter()

        try {
            String content = templateFile.text
            final yaml = (Map)new Yaml().load(content)

            // sanity: restricted to basecalling only, no release of data!
            if (yaml.containsKey('ontReleaseData')) {
                if (yaml['ontReleaseData']) {
                    yaml.remove('ontReleaseData')
                }
            }
            if (yaml.containsKey('ontSetupBasecall')) {
                if (yaml['ontSetupBasecall']) {
                    yaml.remove('ontSetupBasecall')
                }
            }
            if (yaml.containsKey('ontRebasecall')) {
                if (yaml['ontRebasecall']) {
                    yaml.remove('ontRebasecall')
                }
            }
            if (yaml.containsKey('ontSubmitBaseCallJob')) {
                if (yaml['ontSubmitBaseCallJob']) {
                    yaml.remove('ontSubmitBaseCallJob')
                }
            }

            // ontDataFolder may need a default (for older files)
            if (null!=log) { log.info("setupSUPBasecallParamsYAML: settings['workdir_root'] = ${settings['workdir_root']}") }
            if (!yaml.containsKey('ontDataFolder')) {
                yaml['ontDataFolder'] = settings['workdir_root'].toString()
            }

            yaml['ontBaseCall'] = true // to perform basecalling only
            yaml['signalExtensions'] = settings['signalExtension'] // signal file extension
            yaml['ontBaseCallOutFolder'] = derivativeFile.parent // folder to write the basecall results
            // yaml['ontSignalFolders'] = [settings['signalPass']['source'],settings['signalFail']['source']] // folder(s) of signal files
            yaml['ontSignalFolders'] = []
            if (null!=log) { log.info("setupSUPBasecallParamsYAML: yaml['ontSignalFolders'] = ${yaml['ontSignalFolders']}") }
            if (settings.containsKey('signalPass')) {
                yaml['ontSignalFolders'].add(settings['signalPass']['source'])
                if (null!=log) { log.info("setupSUPBasecallParamsYAML: add signalPass, yaml['ontSignalFolders'] = ${yaml['ontSignalFolders']}") }
            }
            if (settings.containsKey('signalFail')) {
                yaml['ontSignalFolders'].add(settings['signalFail']['source'])
                if (null!=log) { log.info("setupSUPBasecallParamsYAML: add signalFail, yaml['ontSignalFolders'] = ${yaml['ontSignalFolders']}") }
            }

            // FIXME: the default model and mod can change with time / project
            // FIXME: the default model and mod are tied to caller
            if (yaml.containsKey('ontBaseCallModel')) {
                // FIXME: what we do if different?
                //if (FlowCellModels.containsKey(settings.runacq.FlowCellProductCode)) {
                    String value = FlowCellModels[settings.runacq.FlowCellProductCode]['supmodel'];
                    if (!yaml['ontBaseCallModel'].equals(value)) {
                        log.warn("Specified basecall model '${yaml.ontBaseCallModel}' differs from expected '${value}'")
                        log.warn("NOT overriding; effective basecall model = '${yaml.ontBaseCallModel}'")
                    }
                //} else {
                //    log.warn("Specified basecall model '${yaml.ontBaseCallModel}' for unknown flow cell product '${settings.runacq.FlowCellProductCode}'")
                //}
            } else {
                // FIXME: use external text file look up configuration!
                //if (FlowCellModels.containsKey("${settings.runacq.FlowCellProductCode}")) {
                    String value = FlowCellModels[settings.runacq.FlowCellProductCode]['supmodel'];
                    yaml['ontBaseCallModel'] = value
                    log.warn("No basecall model specified. Using ${value} per '${settings.runacq.FlowCellProductCode}'")
                //} else {
                //    String value = FlowCellModels['']['supmodel'];
                //    yaml['ontBaseCallModel'] = value
                //    log.warn("No basecall model specified. Using ${value} per default '${settings.runacq.FlowCellProductCode}'. ${FlowCellModels}")
                //    value = FlowCellModels[settings.runacq.FlowCellProductCode]['supmodel'];
                //    log.warn("DEBUG ${value}")
                //}
            }

            if (yaml.containsKey('ontBaseCallBaseModifications')) {
                //if (FlowCellModels.containsKey(settings.runacq.FlowCellProductCode)) {
                    String value = FlowCellModels[settings.runacq.FlowCellProductCode]['modifications'];
                    if (!yaml['ontBaseCallBaseModifications'].equals(value)) {
                        log.warn("Specified base modifications '${yaml.ontBaseCallBaseModifications}' differs from expected '${value}'")
                        log.warn("NOT overriding; effective basecall model = '${yaml.ontBaseCallBaseModifications}'")
                    }
                //} else {
                //    log.warn("Specified base modifications '${yaml.ontBaseCallBaseModifications}' for unknown flow cell product '${settings.runacq.FlowCellProductCode}'")
                //}
            } else {
                String value = FlowCellModels[settings.runacq.FlowCellProductCode]['modifications'];
                yaml['ontBaseCallBaseModifications'] = value
                log.warn("No basecall model specified. Using ${value} per '${settings.runacq.FlowCellProductCode}'")
                //String value = FlowCellModels['']['modifications'];
                //yaml['ontBaseCallBaseModifications'] = value
                //log.warn("No basecall model specified. Using ${value} per default.")
            }


            // record the signal folder inference source
            yaml['ontBamFolders'] = [settings['bamPass']['source'],settings['bamFail']['source']]

            // message to pass to Samplify?
            // we remove the update message as individual runs do not keep metric
            if (yaml.containsKey('rabbitHost')) {
                if (!(yaml['rabbitHost'].equals(''))) {
                    yaml['rabbitHost'] = ''
                }
            }
            if (yaml.containsKey('registration_url')) {
                if (!(yaml['registration_url'].equals(''))) {
                    yaml['registration_url'] = ''
                }
            }

            // localize results
            yaml['sampleDirectory'] = "${ontDataFolder}/${runAcqID}".toString()
            yaml['sampleQCDirectory'] = "${ontDataFolder}/${runAcqID}/qc".toString()

            // FIXME: set from params
            //        prevent spamming user
            if (backendEmailReroute) {
                yaml['userId'] = ONT_BACKEND_USERID
                yaml['userEmail'] = ONT_BACKEND_USEREMAIL
            }

            // inform user : autogen; will be overwritten
            writer.writeLine "##########################################"
            writer.writeLine "# IMPORTANT: auto-generated; any changes WILL BE OVERWRITTEN by pipeline with Samplify's merge action"
            writer.writeLine "#            source="
            writer.writeLine "#            ${template}"
            writer.writeLine "#"
            writer.writeLine "# Configuration SUP basecalling on a single run-acquisition signal data"
            writer.writeLine "#"
            writer.writeLine "##########################################"
            writer.writeLine ""
            writer.writeLine new Yaml().dump(yaml).toString()
            writer.flush()
            writer.close()
        }
        catch( NoSuchFileException e ) {
            throw new IllegalArgumentException("${templateFile} does not exist", e)
        }
        catch( Exception e ) {
            throw new IllegalArgumentException("Error parsing  YAML file: ${templateFile} -- Check the log file for details", e)
        }

        if (log!=null) { log.info("Setup basecall bash param-file ${template} --> ${derivativeFile}") }
    }

    private static getFlowcellProductCode(runAcquisitionPath) {
        String productCode = ''
        String sampleRate = 4000
        File folder = new File(runAcquisitionPath)
        if (folder.exists()) {

            boolean found = false
            boolean foundSR = false
            for (item in folder.listFiles()) {
                if (item.isFile()) {
                    if (item.name.startsWith('report_') && item.name.endsWith('.json')) {
                        // "product_code":"FLO-PRO114M"
                        // "sample_rate":5000
                        Pattern patternPC = Pattern.compile(/"product_code":"([^"]*)"/);
                        Pattern patternSR = Pattern.compile(/"sample_rate":(\d+)/);
                        //def regexPC = /"product_code":"([^"]*)"/;
                        //def regexSR = /"sample_rate":(\d+)/;
                        File file = new File(item.toString())
                        def line
                        file.withReader { reader ->
                            while ((line = reader.readLine()) != null) {
                                Matcher matcher = patternPC.matcher(line);
                                if (matcher.find()) {
                                    productCode = matcher.group(1)
                                    found = true
                                }
                                matcher = patternSR.matcher(line);
                                if (matcher.find()) {
                                    sampleRate = matcher.group(1)
                                    foundSR = true
                                }
                                if (found && foundSR) {
                                    break;
                                }
                            }
                        }
                    } else if (item.name.startsWith('report_') && item.name.endsWith('.md')) {
                        //    "flow_cell_product_code": "FLO-PRO114M", 
                        if (!found) {
                            Pattern pattern = Pattern.compile("\s*\"([^\"]+)\"\s*:\s*\"([^\"]+)\"\s*");
                            File file = new File(item.toString())
                            def line
                            file.withReader { reader ->
                                while ((line = reader.readLine()) != null) {
                                    Matcher matcher = pattern.matcher(line);
                                    if (matcher.find()) {
                                        if (matcher.group(1).equalsIgnoreCase('flow_cell_product_code')) {
                                            productCode = matcher.group(2)
                                            found = true
                                            break
                                        }
                                    }
                                }
                            }
                        }
                    } else if (item.name.startsWith('final_summary_') && item.name.endsWith('.txt')) {
                        //protocol=sequencing/sequencing_PRO114_DNA_e8_2_400T:FLO-PRO114M:SQK-LSK114:400
                        if (!found) {
                            def protocolPrefix = 'protocol='
                            File file = new File(item.toString())
                            def line
                            file.withReader { reader ->
                                while ((line = reader.readLine()) != null) {
                                    if (line.startsWith(protocolPrefix)) {
                                        def keyValue = line.split('=', 2)
                                        def bits = keyValue[1].split(/:/)
                                        productCode = bits[1]
                                        found = true
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if (null!=log) { log.info("getFlowcellProductCode: Product code ${productCode} in ${item.name}? ${found}") }
                    if (found && foundSR) {
                        break
                    }
                }
            }
            if (productCode.equals('')) {
                if (null!=log) { 
                    log.warn("getFlowcellProductCode : Final summary and Report files not found @ ${runAcquisitionPath.toString()}")
                    log.warn("getFlowcellProductCode : To ensure that flow cell product code is correct, please ensure the instrument's run aquisition path is specified.")
                }
            }
        } else {
            if (null!=log) { 
                log.warn("getFlowcellProductCode : The specified instrument's run aquisition path does not exist.")
                log.warn("getFlowcellProductCode : Please check the path ${runAcquisitionPath.toString()}")
            }
        }
        if (productCode.equals('')) {
            return ""
        } else {
            return "${productCode}@${sampleRate}"
        }
    }

    private static getRunAcquisitionMeta(runAcquisitionPath) {
        File folder = new File(runAcquisitionPath)
        String[] folders = runAcquisitionPath.split("/");
        String runAcqID = folder.name
        String[] bits = runAcqID.split("_");
        def meta = [
            'protocolGroupId': folders[-3],
            'flowCellPosition': bits[2],
            'flowCellId': bits[3],
            'protocolRunIdShort': bits[4]
            ]
        return meta
    }

    public static Map setupRunAcquisition(runAcquisitionPath, runAcqBamFolders, sampleId, sampleDirectory, ontDataFolder, ontSubmitBaseCallJob=true, backendEmailReroute=false) {
        /*
            set up the folder if not existing
            set up the subfolder for pod5_{pass,fail} hardlink if exists, else fast5_{pass,fail}
            set up the bam_{pass,fail} hardlink if present and error out if only fastq_{pass,fail} are present!
            set up nwgc-longread-map-merge-qc.yaml if not present
            set up nwgc-longread-map-merge-qc.sh if not present
            Note that overwriting of some marker needed if directory, dorado ver, model or mod are different
        */

        assert wfMeta!=null : "Please call NwgcONTCore.setWorkflowMetadata(workflow) in your workflow definition."

        if (log!=null) { log.info("ontDataFolder = '${ontDataFolder}'") }
        def results = [
            'workdir_root': ontDataFolder
        ]
        File folder = new File(runAcquisitionPath)
        String runAcqID = folder.name
        results['runAcqID'] = runAcqID
        results['basecall'] = ['newsetup' : 0, 'totalsetup' : 0]  // initialize for early termination

        results['runacq'] = ['source': runAcquisitionPath, 'workdir': getONTWorkspaceFolder(ontDataFolder, runAcqID, "")]
        results['runacq']['FlowCellProductCode'] = getFlowcellProductCode(runAcquisitionPath)
        if (log!=null) { log.info("Flowcell Product Code '${results.runacq.FlowCellProductCode}'") }
        if ("" == results['runacq']['FlowCellProductCode']) {
            if (log!=null) { log.info("Skipping run workspace setup for ${runAcqID}; Flowcell Product Code unavailable.") }
            return results
        }

        String signalExtension = SIGNAL_FAST5_EXT
        def unifiedPod5Present = false
        def unifiedPod5PassPresent = false
        def unifiedPod5SkipPresent = false
        def pod5Present = false
        def pod5PassPresent = false
        def pod5FailPresent = false
        def fast5Present = false
        def fast5PassPresent = false
        def fast5FailPresent = false
        folder.eachDir{item ->
            if (item.name == "pod5") {
                unifiedPod5Present = true
                unifiedPod5PassPresent = true
                signalExtension = SIGNAL_POD5_EXT
            } else if (item.name == "pod5_skip") {
                unifiedPod5Present = true
                unifiedPod5SkipPresent = true
                signalExtension = SIGNAL_POD5_EXT
            } else if (item.name == "pod5_pass") {
                pod5Present = true
                pod5PassPresent = true
                signalExtension = SIGNAL_POD5_EXT
            } else if (item.name == "pod5_fail") {
                pod5Present = true
                pod5FailPresent = true
                signalExtension = SIGNAL_POD5_EXT
            } else if (item.name == "fast5_pass") {
                fast5Present = true
                fast5PassPresent = true
                signalExtension = SIGNAL_FAST5_EXT
            } else if ( item.name == "fast5_fail") {
                fast5Present = true
                fast5FailPresent = true
                signalExtension = SIGNAL_FAST5_EXT
            }
        }
        results['signalExtension'] = signalExtension

        if (unifiedPod5Present) {
            assert (unifiedPod5PassPresent || unifiedPod5SkipPresent) : "Missing pod5 or pod5_skip signal folder."
            if (unifiedPod5PassPresent) {
                results['signalPass'] = ['source': getONTInstrumentFolder(runAcquisitionPath, "pod5"), 'unified': true]
                results['signalPass']['workdir'] = getONTWorkspaceFolder(ontDataFolder, runAcqID, "pod5")
            }
            if (unifiedPod5SkipPresent) {
                results['signalFail'] = ['source': getONTInstrumentFolder(runAcquisitionPath, "pod5_skip"), 'unified': true]
                results['signalFail']['workdir'] = getONTWorkspaceFolder(ontDataFolder, runAcqID, "pod5_skip")
            }
        } else if (pod5Present) {
            assert (pod5PassPresent || pod5FailPresent) : "Missing pod5_pass or pod5_fail signal folder."
            if (pod5PassPresent) {
                results['signalPass'] = ['source': getONTInstrumentFolder(runAcquisitionPath, "pod5_pass")]
                results['signalPass']['workdir'] = getONTWorkspaceFolder(ontDataFolder, runAcqID, "pod5_pass")
            }
            if (pod5FailPresent) {
                results['signalFail'] = ['source': getONTInstrumentFolder(runAcquisitionPath, "pod5_fail")]
                results['signalFail']['workdir'] = getONTWorkspaceFolder(ontDataFolder, runAcqID, "pod5_fail")
            }
        } else if (fast5Present) {
            assert (fast5PassPresent || fast5FailPresent) : "Missing fast5_pass or fast5_fail signal folder."
            if (fast5PassPresent) {
                results['signalPass'] = ['source': getONTInstrumentFolder(runAcquisitionPath, "fast5_pass")]
                results['signalPass']['workdir'] = getONTWorkspaceFolder(ontDataFolder, runAcqID, "fast5_pass")
            }
            if (fast5FailPresent) {
                results['signalFail'] = ['source': getONTInstrumentFolder(runAcquisitionPath, "fast5_fail")]
                results['signalFail']['workdir'] = getONTWorkspaceFolder(ontDataFolder, runAcqID, "fast5_fail")
            }
        }

        results['finalPassBam'] = getRunAcquisitionMappedBamFile(ontDataFolder, runAcqID, sampleId, true)

        results['bamFolders'] = runAcqBamFolders
        results['bamPass'] = ['source': getONTInstrumentFolder(runAcquisitionPath, "bam_pass")]
        results['bamPass']['workdir'] = getONTWorkspaceFolder(ontDataFolder, runAcqID, "bam_pass")
        results['bamFail'] = ['source': getONTInstrumentFolder(runAcquisitionPath, "bam_fail")]
        results['bamFail']['workdir'] = getONTWorkspaceFolder(ontDataFolder, runAcqID, "bam_fail")

        // create the <ontDataFolder> if not exist
        File ontDataRunAcqFolder = new File(getONTWorkspaceFolder(ontDataFolder, runAcqID, ""))
        if (!ontDataRunAcqFolder.exists()) {
            ontDataRunAcqFolder.mkdirs()
        }

        def totalSetup = 0
        def newSetup = 0
        def sessionSetups = [:]
        // keep the (hardlinked) copy of signalFolders
        if (unifiedPod5Present) {
            if (unifiedPod5PassPresent) {
                sessionSetups['signalPass'] = setupDataInFolder(results.signalPass.source, results.signalPass.workdir)
                newSetup += sessionSetups['signalPass']['setup']
                totalSetup += sessionSetups['signalPass']['total']
            }
            if (unifiedPod5SkipPresent) {
                sessionSetups['signalFail'] = setupDataInFolder(results.signalFail.source, results.signalFail.workdir)
                newSetup += sessionSetups['signalFail']['setup']
                totalSetup += sessionSetups['signalFail']['total']
            }
        } else if (pod5Present) {
            if (pod5PassPresent) {
                sessionSetups['signalPass'] = setupDataInFolder(results.signalPass.source, results.signalPass.workdir)
                newSetup += sessionSetups['signalPass']['setup']
                totalSetup += sessionSetups['signalPass']['total']
            }
            if (pod5FailPresent) {
                sessionSetups['signalFail'] = setupDataInFolder(results.signalFail.source, results.signalFail.workdir)
                newSetup += sessionSetups['signalFail']['setup']
                totalSetup += sessionSetups['signalFail']['total']
            }
        } else if (fast5Present) {
            if (fast5PassPresent) {
                sessionSetups['signalPass'] = setupDataInFolder(results.signalPass.source, results.signalPass.workdir)
                newSetup += sessionSetups['signalPass']['setup']
                totalSetup += sessionSetups['signalPass']['total']
            }
            if (fast5FailPresent) {
                sessionSetups['signalFail'] = setupDataInFolder(results.signalFail.source, results.signalFail.workdir)
                newSetup += sessionSetups['signalFail']['setup']
                totalSetup += sessionSetups['signalFail']['total']
            }
        }
        
        // keep the (hardlinked) copy of bamFolders
        File workDir = new File("${results.bamPass.source}")
        if (workDir.exists()) {
            sessionSetups['bamPass'] = setupDataInFolder(results.bamPass.source, results.bamPass.workdir)
        }
        workDir = new File("${results.bamFail.source}")
        if (workDir.exists()) {
            sessionSetups['bamFail'] = setupDataInFolder(results.bamFail.source, results.bamFail.workdir)
        }

        // set up bash script and config
        def basecallScript = getONTWorkspaceFolder(ontDataFolder, runAcqID, BASECALL_BASH_SCRIPT)
        results['basecall']['script'] = basecallScript;

        def paramFileFPN = getParamFile(sampleDirectory, wfMeta.commandLine)
        File paramFileUsed = new File(paramFileFPN)
        def basecallParamFile = getONTWorkspaceFolder(ontDataFolder, runAcqID, paramFileUsed.getName())
        results['basecall']['params'] = basecallParamFile

        results['basecall']['newsetup'] = newSetup
        results['basecall']['totalsetup'] = totalSetup

        results['runMeta'] = getRunAcquisitionMeta(runAcquisitionPath)
        results['runMeta']['sampleId'] = sampleId
        results['runMeta']['supBam'] = results['finalPassBam']

        // FIXME: [low risk] need to write directory, dorado-version, model , mod

        // newly setup file(s), force re-processing
        if (newSetup>0) {
            // remove the checksum and bam file
            File workFile = new File("${results.finalPassBam}.md5sum")
            if (workFile.exists()) {
                if (log!=null) { log.info("Re-basecall needed. Deleted ${workFile.toString()}") }
                workFile.delete()
            }
            workFile = new File("${results.finalPassBam}")
            if (workFile.exists()) {
                if (log!=null) { log.info("Re-basecall needed. Deleted ${workFile.toString()}") }
                workFile.delete()
            }
            workFile = new File("${results.finalPassBam}.bai")
            if (workFile.exists()) {
                if (log!=null) { log.info("Re-basecall needed. Deleted ${workFile.toString()}") }
                workFile.delete()
            }
        }

        boolean toSetupScriptAndParameters = false
        if (totalSetup>0) {
            File file = new File(basecallScript)
            if (!file.exists()) {
                if (log!=null) { log.info("${totalSetup} data files detected but the script is missing. To create script.") }
                toSetupScriptAndParameters = true
            }
            file = new File(basecallParamFile)
            if (!file.exists()) {
                if (log!=null) { log.info("${totalSetup} data files detected but the parameter file is missing. To create parameter file.") }
                toSetupScriptAndParameters = true
            }
        }

        if (newSetup>0 || toSetupScriptAndParameters) {
            // results['bamPass']['source'] results['bamFail']['source']

            setupBasecallBashScript(
                "${sampleDirectory}/${PARENT_BASH_SCRIPT}", 
                basecallScript, sampleId, ontDataFolder, runAcqID, totalSetup, basecallParamFile, results['runMeta'], backendEmailReroute)

            setupSUPBasecallParamsYAML(
                paramFileFPN, 
                basecallParamFile, ontDataFolder, runAcqID, results, backendEmailReroute)
        }

        if (ontSubmitBaseCallJob) {
            // FIXME: submit basecalling job to cluster
            //        more appropriate to be done externally
            //        Question: how is local execution supported?
            //defaultClusterOptions = "-S /bin/bash -P $clusterProject -m as -r yes -R yes";
            def proc = "qsub -S /bin/bash -P dna -m as -r yes -R yes -terse ${basecallScript}".execute()
            proc.waitFor()
            def jobid = proc.in.text
            results['basecall']['jobid'] = jobid
            if (log!=null) { log.info("HPC assigned job id ${jobid} to ${basecallScript}") }
        } else {
            results['basecall']['jobid'] = 0
            if (log!=null) {
                log.info("${basecallScript} has not been submitted.")
                log.info("You may submit the job script manually later.")
            }
        }

        return results
    }

    public static String getReleaseLiveModelFolder(ontDataFolder) {
        return "${ontDataFolder}/${RELEASE_LIVEMODEL_FOLDER}"
    }

    public static String getReleaseLiveModelQCFolder(ontDataFolder) {
        // FIXME: synchronize with "param.sampleQCDirectory" ?
        return "${ontDataFolder}/${RELEASE_LIVEMODEL_FOLDER}/qc"
    }

    public static String getReleaseSupFolder(ontDataFolder) {
        return "${ontDataFolder}/${RELEASE_SUP_FOLDER}"
    }

    public static String getReleaseSupQCFolder(ontDataFolder) {
        // FIXME: synchronize with "param.sampleQCDirectory" ?
        return "${ontDataFolder}/${RELEASE_SUP_FOLDER}/qc"
    }

    public static String getReleaseSupPrefix(params) {
        return "${params.sampleId}.${params.sequencingTarget}"
    }

    // write a new bash job wrapper for the re-release task
    // (base off the trigger from "release" action generated bash job wrapper)
    private static setupReleaseDataBashScript(template, derivative, sampleId, ontDataFolder, totalSetup, yamlFile, backendEmailReroute=false) {
        File templateFile = new File(template)
        assert templateFile.exists() : "Source script missing! ${template}"
        File derivativeFile = new File(derivative)
        def writer = derivativeFile.newWriter()

        boolean notesWritten = false
        boolean nxfuuidWritten = false
        boolean workingDirWritten = false
        templateFile.eachLine {
            line -> 
            if (line.startsWith("#!")) {
                writer.writeLine "${line}"
            } else if (line.startsWith("#\$ ")) {
                if (line.startsWith("#\$ -N ")) { // update batch job name
                    // writer.writeLine "${line}_${RELEASE_SUP_FOLDER}"
                    writer.writeLine "#\$ -N s${sampleId}_${RELEASE_SUP_FOLDER}"
                } else if (line.startsWith("#\$ -o ")) { // update stdout folder
                    writer.writeLine "#\$ -o ${ontDataFolder}/${RELEASE_SUP_FOLDER}/logs/"
                } else if (line.startsWith("#\$ -e ")) { // update stderr folder
                    writer.writeLine "#\$ -e ${ontDataFolder}/${RELEASE_SUP_FOLDER}/logs/"
                } else if (line.startsWith("#\$ -M ") && backendEmailReroute) { // update notification email
                    writer.writeLine "#\$ -M ${ONT_BACKEND_USEREMAIL}"
                } else { // pass-thru
                    writer.writeLine "${line}"
                }
            } else {
                if (!notesWritten) {
                    // inform user : autogen; will be overwritten
                    writer.writeLine ""
                    writer.writeLine "##########################################"
                    writer.writeLine "# IMPORTANT: auto-generated; any changes WILL BE OVERWRITTEN by pipeline with Samplify's release action"
                    writer.writeLine "#            source="
                    writer.writeLine "#            ${template}"
                    writer.writeLine "#"
                    writer.writeLine "# Bash script for SUP releasing on the merged run(s)"
                    writer.writeLine "#"
                    writer.writeLine "##########################################"
                    writer.writeLine ""

                    notesWritten = true
                }
                if (line.startsWith("PUBLISH_DIR=")) { // move to run-acquisition folder
                    writer.writeLine "PUBLISH_DIR=${ontDataFolder}/${RELEASE_SUP_FOLDER}"
                } else if (line.startsWith("WORKING_DIR=")) { // for nextflow framework
                    writer.writeLine "${line}_${RELEASE_SUP_FOLDER}"
                } else if (line.startsWith("PIPELINE_YAML=")) { // for nextflow framework
                    writer.writeLine "PIPELINE_YAML=${yamlFile}"
                } else if (line.startsWith("JOB_NAME=")) { // for nextflow framework
                    //writer.writeLine "${line}_${RELEASE_SUP_FOLDER}"
                    String[] bits = line.split('_')
                    String suffix = bits[-1]
                    writer.writeLine "JOB_NAME=s${sampleId}_${RELEASE_SUP_FOLDER}_${suffix}"
                } else if (line.startsWith("USER_EMAIL=") && backendEmailReroute) { // for nextflow framework
                    writer.writeLine "USER_EMAIL=${ONT_BACKEND_USEREMAIL}"
                } else if (line.startsWith("module load ") && -1!=line.indexOf("modules-init")) {
                    writer.writeLine "${line}"
                    writer.writeLine "module load samtools/1.17"
                } else if (line.equals("module load nextflow/22.10.7")) { // module load nextflow/22.10.7 - unavailable!
                    writer.writeLine "module load ${MODULE_NEXTFLOW_VERSION}"
                } else if (line.startsWith("export NXF_UUID=")) { // should have a new session id for samplify database
                    def uuid = UUID.randomUUID()
                    writer.writeLine "export NXF_UUID=${uuid}"
                    nxfuuidWritten = true
                } else if (line.startsWith('cd ${PUBLISH_DIR}')) { // cd ${PUBLISH_DIR}
                    // use the working directory per launchify
                    writer.writeLine 'cd ${WORKING_DIR}'
                    workingDirWritten = true
                } else if (line.startsWith("nextflow ")) { // nextflow run ...
                    if (!nxfuuidWritten) {
                        def uuid = UUID.randomUUID()
                        writer.writeLine "export NXF_UUID=${uuid}"
                        nxfuuidWritten = true
                    }
                    if (!workingDirWritten) {
                        writer.writeLine 'mkdir -p ${WORKING_DIR}'
                        writer.writeLine 'cd ${WORKING_DIR}'
                        workingDirWritten = true
                    }
                    writer.writeLine "${line}"
                } else { // pass-thru
                    writer.writeLine "${line}"
                }
            }
        }
        writer.flush()
        writer.close()

        if (log!=null) { log.info("Setup release bash script ${template} --> ${derivativeFile}") }
    }

    // write a new YAML parameters for the re-release task
    // (base off the trigger from "release" action generated YAML parameters)
    private static setupSUPReleaseDataParamsYAML(template, derivative, settings, backendEmailReroute=false) {

        File templateFile = new File(template)
        assert templateFile.exists() : "Parameters file missing! ${template}"
        File derivativeFile = new File(derivative)
        def writer = derivativeFile.newWriter()

        try {
            String content = templateFile.text
            final yaml = (Map)new Yaml().load(content)

            // sanity: retricted to release only
            yaml['ontReleaseData'] = true
            
            // sanity: no basecalling
            if (yaml.containsKey('ontBaseCall')) {
                if (yaml['ontBaseCall']) {
                    yaml.remove('ontBaseCall')
                }
            }
            if (yaml.containsKey('ontSetupBasecall')) {
                if (yaml['ontSetupBasecall']) {
                    yaml.remove('ontSetupBasecall')
                }
            }
            if (yaml.containsKey('ontRebasecall')) {
                if (yaml['ontRebasecall']) {
                    yaml.remove('ontRebasecall')
                }
            }
            if (yaml.containsKey('ontSubmitBaseCallJob')) {
                if (yaml['ontSubmitBaseCallJob']) {
                    yaml.remove('ontSubmitBaseCallJob')
                }
            }

            // ontDataFolder may need a default (for older files)
            // if one isn't available, workflow will have defaulted

            // let 'ontBamFolders' pass-thru ;
            // No better with overwriting as it must tally with the 'merge' anyway!

            // message to pass to Samplify?
            // pass-thru as release should let Samplify pick up new metrics

            // FIXME: set from params
            //        prevent spamming user
            if (backendEmailReroute) {
                yaml['userId'] = ONT_BACKEND_USERID
                yaml['userEmail'] = ONT_BACKEND_USEREMAIL
            }


            // inform user : autogen; will be overwritten
            writer.writeLine "##########################################"
            writer.writeLine "# IMPORTANT: auto-generated; any changes WILL BE OVERWRITTEN by pipeline with Samplify's release action"
            writer.writeLine "#            source="
            writer.writeLine "#            ${template}"
            writer.writeLine "#"
            writer.writeLine "# Configuration SUP releasing on the merged run(s)"
            writer.writeLine "#"
            writer.writeLine "##########################################"
            writer.writeLine ""
            writer.writeLine new Yaml().dump(yaml).toString()
            writer.flush()
            writer.close()
        }
        catch( NoSuchFileException e ) {
            throw new IllegalArgumentException("${templateFile} does not exist", e)
        }
        catch( Exception e ) {
            throw new IllegalArgumentException("Error parsing  YAML file: ${templateFile} -- Check the log file for details", e)
        }

        if (log!=null) { log.info("Setup release bash param-file ${template} --> ${derivativeFile}") }
    }

    public static Map setupRelease(runAcquisitionsList, sampleId, sampleDirectory, ontDataFolder, outPrefix, backendEmailReroute=false) {
        /*
            set up the folder if not existing
            set up the subfolder for <ontDataFolder>/release_sup
            set up nwgc-longread-map-merge-qc.yaml if not present [ontReleaseData:true]
            set up nwgc-longread-map-merge-qc.sh if not present [publishDir+workingDir+Name]
            TODO: test hardlink file content change effect?
            Note that overwriting of some marker needed if directory, dorado ver, model or mod are different
            Remove the md5sum,bam,bai if newSetup>0 or the bamfolder changes
        */
        // log.info("setupRelease: ${runAcquisitionsList}")
        // TODO:

        File ontDataReleaseFolder = new File(getReleaseSupFolder(ontDataFolder))
        if (!ontDataReleaseFolder.exists()) {
            ontDataReleaseFolder.mkdirs()
        }
        if (log!=null) { log.info("ontDataReleaseFolder = '${ontDataReleaseFolder.toString()}'") }

        def releaseScript = "${ontDataReleaseFolder.toString()}/${RELEASE_BASH_SCRIPT}"
        def paramFileFPN = getParamFile(sampleDirectory, wfMeta.commandLine)
        File paramFileUsed = new File(paramFileFPN)
        def releaseParamFile = "${ontDataReleaseFolder.toString()}/${paramFileUsed.getName()}"

        // newly setup file(s), force re-processing
        def triggerAcqRuns = []
        def newSetup = 0
        def totalSetup = 0
        for (runAcq in runAcquisitionsList) {
            if (runAcq.basecall.newsetup>0) {
                triggerAcqRuns.add(runAcq.runAcqID)
            }
            newSetup += runAcq.basecall.newsetup
            totalSetup += runAcq.basecall.totalsetup
        }
        if (triggerAcqRuns.size()>0) {
            log.info("The following ${triggerAcqRuns.size()} run(s) has changed and will trigger re-release: ${triggerAcqRuns}")
        }

        // TODO: name of output file!
        // this has no effect as we link to the sampleDir
        if (newSetup>0) {
            // remove the checksum and bam file
            def releaseBam = "${ontDataReleaseFolder}/${outPrefix}.bam"
            File workFile = new File("${releaseBam}.md5sum")
            if (workFile.exists()) {
                if (log!=null) { log.info("Re-release needed. Deleted ${workFile.toString()}") }
                workFile.delete()
            }
            workFile = new File("${releaseBam}")
            if (workFile.exists()) {
                if (log!=null) { log.info("Re-release needed. Deleted ${workFile.toString()}") }
                workFile.delete()
            }
            workFile = new File("${releaseBam}.bai")
            if (workFile.exists()) {
                if (log!=null) { log.info("Re-release needed. Deleted ${workFile.toString()}") }
                workFile.delete()
            }
        }

        // release should NOT be done if there is no signal file!
        // To ensure that every "Merge" request on Samplify will trigger a new copy of .yaml & .sh for "Final Merge"
        boolean toSetup = true
        /*
        boolean toSetup = false
        if (totalSetup>0) {
            File file = new File(releaseScript)
            if (!file.exists()) {
                if (log!=null) { log.info("${totalSetup} data files detected but the script is missing. To create script.") }
                toSetup = true
            }
            file = new File(releaseParamFile)
            if (!file.exists()) {
                if (log!=null) { log.info("${totalSetup} data files detected but the parameter file is missing. To create parameter file.") }
                toSetup = true
            }
        }
        */

        if (newSetup>0 || toSetup) {
            setupReleaseDataBashScript(
                "${sampleDirectory}/${PARENT_BASH_SCRIPT}", 
                releaseScript, sampleId, ontDataFolder, totalSetup, releaseParamFile, backendEmailReroute)

            setupSUPReleaseDataParamsYAML(
                paramFileFPN, 
                releaseParamFile, runAcquisitionsList, backendEmailReroute)

            // NOTE: do not auto-submit; to be triggered by user
        }

        // TODO: return settings?
    }

    private static boolean isAlignFileSUPModel(alignFile) {
        File sambamFile = new File(alignFile)
        assert sambamFile.exists() : "Alignment file missing! ${alignFile}"

        def RGID = ''
        // FIXME: detect samtools presence
        def p = "samtools view -H ${alignFile}".execute()
        def outThread = p.consumeProcessOutputStream(new LineOutputStream({line -> 
            if (RGID.equals('') && line.startsWith('@RG')) {
                def cols = line.split(/\t/)
                for(col in cols) { 
                    if (col.startsWith('ID:')) {
                        RGID = col.substring(3)
                        break
                    }
                }
            }}))
        def errThread = p.consumeProcessErrorStream(new LineOutputStream({line -> println "isAlignFileSUPModel: error -> $line"}))
        try { outThread.join(); } catch (InterruptedException ignore) {}
        try { errThread.join(); } catch (InterruptedException ignore) {}
        try { p.waitFor(); } catch (InterruptedException ignore) {}

        if (null!=log) { log.info("isAlignFileSUPModel: RGID = ${RGID}") }
        def results = (-1!=RGID.toLowerCase().indexOf("_sup"))
        return results
    }

    private static boolean isToBackup(folder, alignFilename, qcfolder, bkupFolder, bkupQCFolder) {
        File srcFile = new File("${folder}/${alignFilename}")
        // nothing to backup
        if (!srcFile.exists()) {
            if (null!=log) { log.info("NO backup: source ${alignFilename} non-existent") }
            return false
        }

        File destFile = new File("${bkupFolder}/${alignFilename}")
        // no backup
        if (!destFile.exists()) {
            if (null!=log) { log.info("TO BACKUP: destination ${alignFilename} non-existent") }
            return true
        }

        // different size
        if (srcFile.length() != destFile.length()) {
            if (null!=log) { log.info("TO BACKUP: ${alignFilename} source size (${srcFile.length()}) destination size (${destFile.length()})") }
            return true
        }

        srcFile = new File("${folder}/${alignFilename}.md5sum")
        destFile = new File("${bkupFolder}/${alignFilename}.md5sum")
        // same
        if (srcFile.exists() && destFile.exists()) {
            def srcMD5 = srcFile.text.split(/\n/)
            def destMD5 = destFile.text.split(/\n/)
            if (srcMD5[0].equals(destMD5[0])) {
                if (null!=log) { log.info("NO BACKUP: ${alignFilename} identical md5 (${srcMD5[0]})") }
                return false
            } else {
                if (null!=log) { log.info("TO BACKUP: ${alignFilename} source md5 (${srcMD5[0]}) destination md5 (${destMD5[0]})") }
                return true
            }
        } else {
            if (null!=log) { log.info("TO BACKUP: ${alignFilename} incomplete/missing") }
        }


        // TODO: future check the qc folder!
        return true
    }

    public static String handleLiveModelBackup(params) {
        def outPrefix = getReleaseSupPrefix(params)
        log.info("handleLiveModelBackup: ${params.sampleDirectory}/${outPrefix}.bam")

        File alignFile = new File("${params.sampleDirectory}/${outPrefix}.bam")
        if (!alignFile.exists()) {
            if (null!=log) {
                log.info("No ${outPrefix}.bam @ ${params.sampleDirectory}")
                log.info("1. No processing on the livemodel data done?")
                log.info("2. Partial processing?")
                log.info("3. Processing in-transit?")
            }
        } else {
            if (isAlignFileSUPModel(alignFile.toString())) {
                if (null!=log) {
                    log.info("SUP model identified for ${outPrefix}.bam @ ${params.sampleDirectory}")
                    log.info("No livemodel data backup necessary.")
                }
            } else {
                // to perform backup
                if (null!=log) {
                    log.info("Non-sup model identified for ${outPrefix}.bam @ ${params.sampleDirectory}")
                    log.info("Assume livemodel data. Backing up..")
                }

                def ontDataFolder = getONTDataFolder(params)
                def releaseLiveModelFolder = getReleaseLiveModelFolder(ontDataFolder)
                def releaseLiveModelQCFolder = getReleaseLiveModelQCFolder(ontDataFolder)
                if (isToBackup(params.sampleDirectory, "${outPrefix}.bam", params.sampleQCDirectory,
                    releaseLiveModelFolder, releaseLiveModelQCFolder)) {
                    // longread-map-merged-qc(hardlink/copy) --> ont/map-merge-qc_livemodel
                    log.info("Change(s) detected. Initiating backup.")

                    // handle QC folder first
                    def sourceFolderLength = params.sampleQCDirectory.length()+1
                    def sourceFolder = new File(params.sampleQCDirectory)
                    // let's create all folder(s) first!
                    sourceFolder.eachFileRecurse { 
                        item ->
                        if (item.isDirectory()) {
                            // make destination directory
                            def suffix = item.toString().substring(sourceFolderLength)
                            def destFolder = new File("${releaseLiveModelQCFolder}/${suffix}")
                            if (!destFolder.exists()) {
                                log.info("DIR ${suffix} created.")
                                destFolder.mkdirs()
                            }
                        }
                    }
                    // then, work on the file(s)
                    sourceFolder.eachFileRecurse { 
                        item ->
                        if (item.isFile()) {
                            // move file
                            def suffix = item.toString().substring(sourceFolderLength)
                            def status = item.renameTo(new File("${releaseLiveModelQCFolder}/${suffix}"))
                            if (!status && item.exists()) {
                                log.warn("Renaming failed, file still at ${item.toString()}")
                            } else {
                                log.info("FILE ${suffix} moved.")
                            }
                        }
                    }

                    // handle top level folder
                    sourceFolderLength = params.sampleDirectory.length()+1
                    sourceFolder = new File(params.sampleDirectory)
                    // let's create all folder(s) first!
                    sourceFolder.eachFileRecurse { 
                        item ->
                        if (item.isDirectory()) {
                            // make destination directory
                            // omit folder: .nextflow, monitor
                            def suffix = item.toString().substring(sourceFolderLength)
                            def destFolder = new File("${releaseLiveModelFolder}/${suffix}")
                            if (!destFolder.exists()) {
                                log.info("DIR ${suffix} created.")
                                destFolder.mkdirs()
                            }
                        }
                    }
                    // then, work on the file(s)
                    sourceFolder.eachFileRecurse { 
                        item ->
                        if (item.isFile()) {
                            // move file:
                            // copy: *.sh and *.yaml
                            def suffix = item.toString().substring(sourceFolderLength)
                            if (suffix.endsWith(".sh") || suffix.endsWith(".yaml")) {
                                def status = item.copyTo(new File("${releaseLiveModelFolder}/${suffix}"))
                                if (!status) {
                                    log.warn("Copy failed, file still at ${item.toString()}")
                                } else {
                                    log.info("FILE ${suffix} copied.")
                                }
                            } else {
                                def status = item.renameTo(new File("${releaseLiveModelFolder}/${suffix}"))
                                if (!status && item.exists()) {
                                    log.warn("Renaming failed, file still at ${item.toString()}")
                                } else {
                                    log.info("FILE ${suffix} moved.")
                                }
                            }
                        }
                    }
                } else {
                    log.info("No change detected. Not performing backup.")
                }
            }
        }
    }

}

public class LineOutputStream extends OutputStream {
    private Closure callback;
    private StringBuilder s = new StringBuilder();

    public LineOutputStream(Closure callback) {
        this.callback = callback
    }

    public synchronized void close() {
        if(s.size() > 0) {
            String line = s.toString();
            if(callback != null) {
                callback(line);
            }
        }
        s = null;
    }

    public synchronized void write(int b) {
        if((char)b == '\r') {}
        else if((char)b == '\n') {
            String line = s.toString();
            s = new StringBuilder();
            if(callback != null) {
                callback(line);
            }
        } else {
            s.append((char)b);
        }
    }
}

