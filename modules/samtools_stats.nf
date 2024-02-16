process SAMTOOLS_STATS {

    label "SAMTOOLS_STATS_${params.sampleId}_${params.userId}"

    publishDir "${qcFolder}", mode: 'link', pattern: '*.stats.txt'

    input:
        path bam
        path bai
        val(qcFolder)

    output:
        path "*.stats.txt", emit: stats
        path "versions.yaml", emit: versions

    script:
        """
        # mkdir -p $qcFolder

        samtools \
            stats \
            $bam \
            --threads $task.cpus \
            > ${params.sampleId}.samtools.stats.txt \

        cat <<-END_VERSIONS > versions.yaml
        '${task.process}':
            samtools: \$(samtools --version | grep ^samtools | awk '{print \$2}')
        END_VERSIONS

        """

}
