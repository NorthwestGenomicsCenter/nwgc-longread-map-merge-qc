process PICARD_COVERAGE_METRICS {

    label "PICARD_COVERAGE_METRICS_${params.sampleId}_${params.userId}"

    Boolean isOnt = params.sequencingPlatform.equalsIgnoreCase("ont")

    publishDir "${qcFolder}", mode: 'link', pattern: '*.picard.coverage.txt'
 
    input:
        path bam
        path bai
        val(qcFolder)

    output:
        path "*.picard.coverage.txt", emit: stats
        path "versions.yaml", emit: versions

    script:

        def minimumBaseQuality = params.sequencingPlatform.equalsIgnoreCase("ont") ? '10' : '20'

        """
        # mkdir -p $qcFolder

        java \
            -XX:InitialRAMPercentage=80 \
            -XX:MaxRAMPercentage=85 \
            -jar \$PICARD_DIR/picard.jar CollectWgsMetrics \
            --INPUT $bam \
            --COUNT_UNPAIRED true \
            --READ_LENGTH 17000 \
            --INCLUDE_BQ_HISTOGRAM \
            --REFERENCE_SEQUENCE $params.referenceGenome \
            --VALIDATION_STRINGENCY LENIENT \
            --MINIMUM_BASE_QUALITY $minimumBaseQuality \
            --OUTPUT ${params.sampleId}.picard.coverage.txt

        cat <<-END_VERSIONS > versions.yaml
        '${task.process}':
            java: \$(java -version 2>&1 | grep version | awk '{print \$3}' | tr -d '"''')
            picard: \$(java -jar \$PICARD_DIR/picard.jar CollectRawWgsMetrics --version 2>&1 | awk '{split(\$0,a,":"); print a[2]}')
        END_VERSIONS

        """

}
