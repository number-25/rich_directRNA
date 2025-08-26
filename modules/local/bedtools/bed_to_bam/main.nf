process BED_TO_BAM {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2':
        'biocontainers/bedtools:2.31.1--hf5e1c6e_2' }"
    //publishDir "

    input:
    tuple val(meta), path(collapsed_bed)
    path(genome_sizes)

    output:
    tuple val(meta), path("*.bam"), emit: collapsed_bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_collapsed"
    """
    sort -k 1,1 -k2,2n ${collapsed_bed} \\
    | bedToBam \\
        -i - \\
        -bed12 \\
        -g ${genome_sizes} \\
        > ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version |& sed 's/bedtools v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_collapsed"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version |& sed 's/bedtools v//')
    END_VERSIONS
    """
}
