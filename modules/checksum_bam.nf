process CHECKSUM_BAM {

    label "CHECKSUM_BAM_${params.sampleId}_${params.userId}"

    publishDir "${checksumPath}", mode:  'link', pattern: "${inputBam}.md5sum"

    input:
        path inputBam
        val(checksumPath)

    output:
        path "${inputBam}", emit: bam
        path "${inputBam}.md5sum", emit: md5sum
        path "versions.yaml", emit: versions

    script:
        """
        md5sum $inputBam | awk '{print \$1}' > ${inputBam}.md5sum

        cat <<-END_VERSIONS > versions.yaml
        '${task.process}':
            md5sum: \$(md5sum --version | head -n 1 | awk '{print \$4}')
        END_VERSIONS

        """

}
