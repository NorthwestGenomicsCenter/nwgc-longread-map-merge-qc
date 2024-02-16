process MAP_HIFI_BAM {

    label "MAP_HIFI_BAM_${params.sampleId}_${params.userId}"

    input:
        path hiFiBam

    output:
        path "*.mapped.bam",  emit: mapped_bam
        path "versions.yaml", emit: versions

    script:
        def numCPUs = Integer.valueOf("$task.cpus")
        def align_threads = numCPUs
        def sort_threads = Math.min(8, Math.ceil(numCPUs/4).intValue())

        """
        pbmm2 \\
            align \\
            --num-threads ${task.cpus} \\
            --unmapped \\
            --sort -J $sort_threads -m 2G \\
            $params.referenceGenome \\
            $hiFiBam \\
            ${hiFiBam}.mapped.bam

        cat <<-END_VERSIONS > versions.yaml
        '${task.process}_${task.index}':
            pbmm2: \$(pbmm2 --version | awk '{print \$2}')
        END_VERSIONS
        """

}
