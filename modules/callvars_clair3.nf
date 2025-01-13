process CALLVARS_CLAIR3 {

    label "CALLVARS_CLAIR3_${params.sampleId}_${params.userId}"

    publishDir "$params.sampleDirectory", mode: 'link', pattern: "*"

    input:
        tuple (
        path bam,
        path bai,
        path referenceGenome,
        )

    output:
        path "*.clair3.vcf"

    script:

    """

        /opt/bin/run_clair3.sh \
            --bam_fn=${params.bam} \
            --sample_name=${params.sampleID} \
            --ref_fn=${params.referenceGenome} \
            --threads=${task.cpus} \
            --platform=${params.platform} \
            --model_path=/opt/models/${params.modelName} \
            --enable_phasing \
            --longphase_for_phasing \
            --use_whatshap_for_final_output_haplotagging \
            --gvcf \
            --output=${params.sampleId}.${params.sequencingTarget}.clair3

    """
}