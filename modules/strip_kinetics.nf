process STRIP_KINETICS {

    label "STRIP_KINETICS_${params.sampleId}_${params.userId}"

    input:
        path bam

    output:
        path "*.primrose.bam",  emit: bam
        path "*.primrose.bam.bai",  emit: bai
        path "versions.yaml", emit: versions

    script:

        """
        JASMINE_WAS_RUN=false
        if samtools view -H $bam | grep ^@PG | grep -q ID:jasmine ; then JASMINE_WAS_RUN=true; fi
        JASMINE_KEPT_KINETECS=false
        if samtools view -H $bam | grep ^@PG | grep ID:jasmine | grep -q keep-kinetics; then JASMINE_KEPT_KINETECS=true; fi

        PRIMROSE_WAS_RUN=false
        if samtools view -H $bam | grep ^@PG | grep -q ID:primrose ; then PRIMROSE_WAS_RUN=true; fi
        PRIMROSE_KEPT_KINETECS=false
        if samtools view -H $bam | grep ^@PG | grep ID:primrose | grep -q keep-kinetics; then PRIMROSE_KEPT_KINETECS=true; fi

        RUN_PRIMROSE=false
        if [[ "\$JASMINE_WAS_RUN" = true ]] ; then
            if [[ "\$JASMINE_KEPT_KINETECS" = true ]] ; then
                RUN_PRIMROSE=true
            else
                echo "Jasmine already ran and removed kinetics"
            fi
        elif  [[ "\$PRIMROSE_WAS_RUN" = true ]] ; then
            if [[ "\$PRIMROSE_KEPT_KINETECS" = true ]] ; then
                RUN_PRIMROSE=true
            else
                echo "Primrose already ran and removed kinetics"
            fi
        else 
            echo "Neither Jasmine or Primrose has run"
            RUN_PRIMROSE=true
        fi
            
        if [[ "\$RUN_PRIMROSE" = true ]] ; then
            primrose \
                -j $task.cpus \
                $bam \
                ${bam}.primrose.bam;
        else
            ln $bam ${bam}.primrose.bam;
        fi

        samtools \
            index \
            -@ $task.cpus \
            ${bam}.primrose.bam

        cat <<-END_VERSIONS > versions.yaml
        '${task.process}_${task.index}':
            primrose: \$(primrose --version | grep ^primrose | awk '{print \$2}')
            samtools: \$(samtools --version | grep ^samtools | awk '{print \$2}')
        END_VERSIONS
        """

}
