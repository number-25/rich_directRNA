process GFFCOMPARE {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffcompare:0.12.6--h9f5acd7_1':
        'biocontainers/gffcompare:0.12.6--h9948957_4' }"

    input:
    tuple path(genome_fasta), path(genome_fasta_index)
    tuple val(meta), path(reconstructed_gtf)
    path annotation_gtf
    val origin

    output:
    tuple val(meta), path("*.stats"),   emit: gffcompare_stats
    path "*.annotated.gtf"
    path "*.tracking"
    path "*.refmap"
    path "*.tmap"
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}_${meta.replicate}_${origin}"

    """
    gffcompare \\
        ${args} \\
        -o ${prefix} \\
        -r ${annotation_gtf} \\
        -s ${genome_fasta} \\
        ${reconstructed_gtf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffcompare: \$(gffcompare --version |& sed '1!d ; s/ //')
    END_VERSIONS
    """

    stub:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}_${meta.replicate}_${origin}"
    """
    touch ${prefix}.stats
    touch ${prefix}.annotated.gtf
    touch ${prefix}.tracking
    touch ${prefix}.refmap
    touch ${prefix}.tmap

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffcompare: \$(gffcompare --version |& sed '1!d ; s/ //')
    END_VERSIONS
    """
}
