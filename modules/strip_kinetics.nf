process STRIP_KINETICS {

    label "STRIP_KINETICS_${params.sampleId}_${params.userId}"

    input:
        path bam

    output:
        path "*.methylation_prediction.bam",  emit: bam
        path "*.methylation_prediction.bam.bai",  emit: bai
        path "versions.yaml", emit: versions

    script:

        def methylationPredictionProgram = params.sequencerType == 'Revio' ? 'jasmine' : 'primrose'

        """
        METHYLATION_PREDICTION_WAS_RUN=true
        if ! samtools view -H $bam | grep ^@PG | grep -q ID:$methylationPredictionProgram ; then METHYLATION_PREDICTION_WAS_RUN=false; fi

        METHYLATION_PREDICTION_KEPT_KINETECS=false
        if samtools view -H $bam | grep ^@PG | grep ID:$methylationPredictionProgram | grep -q keep-kinetics; then METHYLATION_PREDICTION_KEPT_KINETECS=true; fi

        if [[ "\$METHYLATION_PREDICTION_WAS_RUN" = false || "\$METHYLATION_PREDICTION_KEPT_KINETECS" = true ]] ; then
            $methylationPredictionProgram \
                -j $task.cpus \
                $bam \
                ${bam}.methylation_prediction.bam;
        else
                ln $bam ${bam}.methylation_prediction.bam;
        fi

        samtools \
            index \
            -@ $task.cpus \
            ${bam}.methylation_prediction.bam

        cat <<-END_VERSIONS > versions.yaml
        '${task.process}_${task.index}':
            $methylationPredictionProgram: \$($methylationPredictionProgram --version | grep ^methylation_prediction | awk '{print \$2}')
            samtools: \$(samtools --version | grep ^samtools | awk '{print \$2}')
        END_VERSIONS
        """

}
