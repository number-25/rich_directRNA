process BAM_TO_BEDGRAPH {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_2':
        'biocontainers/bedtools:2.31.1--hf5e1c6e_2' }"
    //publishDir "

    input:
    tuple val(meta), path(bam)
    path genome_sizes
    val strand

    output:
    tuple val(meta), path("*.bedgraph"), emit: bedgraph
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_${strand}"
    """
    bedtools \
        genomecov \
        -ibam $bam \
        -bg \
        -strand $strand \
        | bedtools sort \
        > $prefix.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version |& sed 's/bedtools v//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_${strand}"
    """
    touch ${prefix}.bedgraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version |& sed 's/bedtools v//')
    END_VERSIONS
    """
}
