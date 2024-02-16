
process MERGE_MAPPED_BAMS {

    label "MERGE_MAPPED_BAMS_${params.sampleId}_${params.userId}"

    publishDir "${mergedPath}", mode:  'link', pattern: "${saveAsPrefix}.bam"
    publishDir "${mergedPath}", mode:  'link', pattern: "${saveAsPrefix}.bam.bai"

    input:
        path bamList
        val(mergedPath)
        val(saveAsPrefix)

    output:
        path "${saveAsPrefix}.bam",  emit: merged_sorted_bam
        path "${saveAsPrefix}.bam.bai",  emit: bai
        path "versions.yaml", emit: versions

    script:

        """
        samtools \
            merge \
            --threads $task.cpus \
            -c \
            $bamList \
            -o ${saveAsPrefix}.bam

        samtools \
            index \
            -@ $task.cpus \
            ${saveAsPrefix}.bam

        cat <<-END_VERSIONS > versions.yaml
        '${task.process}':
            samtools: \$(samtools --version | grep ^samtools | awk '{print \$2}')
        END_VERSIONS

        """

}
