process CONTAMINATION_CHECK {

    label "CONTAMINATION_CHECK_${params.sampleId}_${params.userId}"

    publishDir "${qcFolder}", mode: 'link', pattern: '*.VerifyBamID.selfSM'

    input:
        path bam
        path bai
        val(qcFolder)

    output:
        path "*.VerifyBamID.selfSM", emit: verifications
        path "versions.yaml", emit: versions

    script:
        def disableSanityCheck = params.mode == 'test' ? '--DisableSanityCheck' : ''

        """
        # mkdir -p $qcFolder

        VERIFYBAMID_RESOURCE=\$MOD_GSVERIFYBAMID_DIR/resource

        UDPATH=\$VERIFYBAMID_RESOURCE/1000g.phase3.100k.b38.vcf.gz.dat.UD
        BEDPATH=\$VERIFYBAMID_RESOURCE/1000g.phase3.100k.b38.vcf.gz.dat.bed
        MEANPATH=\$VERIFYBAMID_RESOURCE/1000g.phase3.100k.b38.vcf.gz.dat.mu

        samtools \
            mpileup \
            -f $params.referenceGenome \
            -B \
            -l \
            \$BEDPATH \
            $bam \
            -o ${params.sampleId}.pileup

        VerifyBamID \
            --PileupFile ${params.sampleId}.pileup \
            --UDPath \$UDPATH \
            --BedPath \$BEDPATH \
            --MeanPath \$MEANPATH \
            --Reference $params.referenceGenome \
            --Verbose \
            $disableSanityCheck \
            --Output ${params.sampleId}.VerifyBamID

        cat <<-END_VERSIONS > versions.yaml
        '${task.process}':
            samtools: \$(samtools --version | grep ^samtools | awk '{print \$2}')
            VerifyBamId: \$(module list 2>&1 | awk '{for (i=1;i<=NF;i++){if (\$i ~/^VerifyBamID/) {print \$i}}}' | awk -F/ '{print \$NF}')
        END_VERSIONS

        """

}
