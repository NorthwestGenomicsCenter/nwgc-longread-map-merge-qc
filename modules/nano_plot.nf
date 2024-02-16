process NANO_PLOT {

    label "NANO_PLOT_${params.sampleId}_${params.userId}"

    publishDir "${qcFolder}", mode: 'link', pattern: 'nanoPlot/*.html'
    publishDir "${qcFolder}", mode: 'link', pattern: 'nanoPlot/*.png'
    publishDir "${qcFolder}", mode: 'link', pattern: 'nanoPlot/NanoStats.txt'
 
    input:
        path bam
        path bai
        val(qcFolder)

    output:
        path "nanoPlot/*.html", emit: html
        path "nanoPlot/*.png", emit: png
        path "nanoPlot/NanoStats.txt", emit: stats
        path "versions.yaml", emit: versions

    script:
        """
        # echo qcFolder ${qcFolder}
        # mkdir -p ${qcFolder}/nanoPlot

        NanoPlot \
            -t $task.cpus \
            -o nanoPlot \
            --bam $bam

        cat <<-END_VERSIONS > versions.yaml
        '${task.process}':
            python: \$(python --version | awk '{print \$2}')
            NanoPlot: \$(NanoPlot --version | awk '{print \$2}')
        END_VERSIONS

        """

}
