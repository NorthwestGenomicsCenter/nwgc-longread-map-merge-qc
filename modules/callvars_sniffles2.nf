process CALLVARS_SNIFFLES2 {

    label "CALLVARS_ONT_SNIFFLES2_${params.sampleId}_${params.userId}"

    publishDir "$params.sampleDirectory", mode: 'link', pattern: "*"

    input:
        path bam
        path bai
        path referenceGenome

    output:
        path "*.sniffles.vcf"
        path "*.sniffles.snf"
    script:

    """

        /usr/local/bin/sniffles \
            --input ${params.bam} \
            --ref_fn=${params.referenceGenome} \
            --threads=${task.cpus} \
            --vcf ${params.sampleId}.sniffles.vcf \
            --snf ${params.sampleId}.sniffles.snf \
            --tandem-repeats ${params.tandemRepeatBed} \
            --reference ${params.referenceGenome}

    """
}