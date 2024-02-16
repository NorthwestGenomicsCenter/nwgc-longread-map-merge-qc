process MAP_ONT_FASTQ {

    label "MAP_ONT_FASTQ_${params.sampleId}_${params.userId}"

    input:
        path fastq

    output:
        path "*.mapped.bam",  emit: mapped_bam
        path "versions.yaml", emit: versions

    script:
        """
        minimap2 \\
            -a \\
            -x map-ont \\
            -t ${task.cpus} \\
            --MD \\
            $params.referenceGenome \\
            $fastq \\
        | \\
        samtools \\
            sort \\
            - \\
            -@ 3 -m 1G \\
            -o ${fastq}.mapped.bam


        cat <<-END_VERSIONS > versions.yaml
        '${task.process}_${task.index}':
            minimap2: \$(minimap2 --version)
            samtools: \$(samtools --version | grep ^samtools | awk '{print \$2}')
        END_VERSIONS
        """

}
