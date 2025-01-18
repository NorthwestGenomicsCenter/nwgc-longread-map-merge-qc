process MOSDEPTH {

    label "MOSDEPTH_${params.sampleId}_${params.userId}"

    // For simplicity, ONT's merge_HAC and merge_SUP publish to the same sample directory
    // Thus, overwrite is allowed
    publishDir "${qcFolder}", mode: 'link', pattern: '*.mosdepth.summary.txt', overwrite: true
    publishDir "${qcFolder}", mode: 'link', pattern: '*.thresholds.bed.gz', overwrite: true
    publishDir "${qcFolder}", mode: 'link', pattern: '*.thresholds.bed.gz.csi', overwrite: true
 
    input:
        path bam
        path bai
        val(qcFolder)

    output:
        path "*.mosdepth.summary.txt", emit: summary
        path "*.thresholds.bed.gz", emit: thresholds
        path "*.thresholds.bed.gz.csi", emit: thresholds_index
        path "versions.yaml", emit: versions

    script:

        """
        # mkdir -p $qcFolder

        mosdepth \
            --mapq 20 \
            --thresholds 1,5,10,15,20,25,30,40,50,60,70,80,90,100 \
            --no-per-base \
            --fast-mode \
            --by ${params.sequencingTargetBed} \
            ${params.sampleId}.mosdepth \
            $bam

        cat <<-END_VERSIONS > versions.yaml
        '${task.process}':
            mosdepth: \$(mosdepth --version | awk '{print \$2}')
        END_VERSIONS

        """

}
