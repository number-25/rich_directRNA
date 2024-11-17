process BAM_TO_BED12 {
    tag "$meta.id"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://brookslab/flair:2.0.0' :
        'biocontainers/flair:2.0.0--pyhdfd78af_1' }"
        ///'quay.io/biocontainers/flair:2.0.0--pyhdfd78af_1' }"
    input:
    tuple val(meta), path(mapped_bam)
    tuple val(meta), path(mapped_bam_index)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}"
    """
    bam2Bed12 \\
        -i ${mapped_bam} \\
        --keep_supplementary \\
        | sort -k 1,1 -k2,2n \\
        > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam_to_bed12: \$(flair --version |& sed 's/flair //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bam_to_bed12: \$(flair --version |& sed 's/flair //')
    END_VERSIONS
    """
}
