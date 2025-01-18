process CHECKSUM_BAM {

    label "CHECKSUM_BAM_${params.sampleId}_${params.userId}"

    // For simplicity, ONT's merge_HAC and merge_SUP publish to the same sample directory
    // Thus, overwrite is allowed
    publishDir "${checksumPath}", mode: 'link', pattern: "${bam}.md5sum", overwrite: true
    publishDir "${checksumPath}", mode: 'link', pattern: "${bai}.md5sum", overwrite: true

    input:
        path bam
        path bai
        val(checksumPath)

    output:
        path "${bam}", emit: bam
        path "${bam}.md5sum", emit: bammd5sum
        path "${bai}", emit: bai
        path "${bai}.md5sum", emit: baimd5sum
        path "versions.yaml", emit: versions

    script:
        """
        md5sum $bam | awk '{print \$1}' > ${bam}.md5sum
        md5sum $bai | awk '{print \$1}' > ${bai}.md5sum

        cat <<-END_VERSIONS > versions.yaml
        '${task.process}':
            md5sum: \$(md5sum --version | head -n 1 | awk '{print \$4}')
        END_VERSIONS

        """

}
