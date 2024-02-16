process ADD_NM_TAGS {

    label "ADD_NM_TAGS_${params.sampleId}_${params.userId}"

    input:
        path mappedBam

    output:
        path "*.nm.bam",  emit: nm_bam
        path "versions.yaml", emit: versions

    script:

        """
        samtools \
            calmd \
            -b \
            --threads $task.cpus \
            $mappedBam \
            $params.referenceGenome \
            > ${mappedBam}.nm.bam

        cat <<-END_VERSIONS > versions.yaml
        '${task.process}':
            samtools: \$(samtools --version | grep ^samtools | awk '{print \$2}')
        END_VERSIONS

        """

}
