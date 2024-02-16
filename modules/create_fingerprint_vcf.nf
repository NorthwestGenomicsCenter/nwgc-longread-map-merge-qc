process CREATE_FINGERPRINT_VCF {

    label "CREATE_FINGERPRINT_VCF_${params.sampleId}_${params.userId}"

    publishDir "${qcFolder}", mode: 'link', pattern: '*.fingerprint.vcf.gz'
    publishDir "${qcFolder}", mode: 'link', pattern: '*.fingerprint.vcf.gz.tbi'
 
    input:
        path bam
        path bai
        val(qcFolder)

    output:
        path "*.fingerprint.vcf.gz", emit: vcfs
        path "*.fingerprint.vcf.gz.tbi", emit: vcfindexes
        path "versions.yaml", emit: versions

    script:
    
        def bcftools_platform_specific_configuration_profile = params.sequencingPlatform.equalsIgnoreCase("pacbio") ? "pacbio-ccs" : "ont"

        """
        # mkdir -p $qcFolder

        bcftools \
            mpileup \
            -X $bcftools_platform_specific_configuration_profile \
            -f $params.referenceGenome \
            -R $params.fingerprintBed \
            $bam \
        | \
        bcftools \
            call \
            -m \
            -Ov \
            --annotate=GQ \
            -o ${params.sampleId}.fingerprint.vcf

        bgzip -f ${params.sampleId}.fingerprint.vcf

        tabix -f ${params.sampleId}.fingerprint.vcf.gz

        cat <<-END_VERSIONS > versions.yaml
        '${task.process}':
            bcftools: \$(bcftools --version | grep ^bcftools | awk '{print \$2}')
            htslib: \$(bcftools --version | grep htslib | awk '{print \$3}')
            tabix: \$(tabix 2>&1  | grep Version | awk '{print \$2}')
        END_VERSIONS

        """

}
