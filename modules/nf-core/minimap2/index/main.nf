process MINIMAP2_INDEX {
    tag "$fasta"
    label 'process_high'

    // Note: the versions here need to match the versions used in minimap2/align
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.28--he4a0461_0' :
        'biocontainers/minimap2:2.28--he4a0461_0' }"

    input:
    tuple val(meta), path(genome_fasta)

    output:
    tuple val(meta), path("*.mmi"), emit: index
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // add $args to script?
    """
    minimap2 \\
        -t $task.cpus \\
        -d ${genome_fasta.baseName}.mmi \\
        $genome_fasta \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    """
    touch ${genome_fasta.baseName}.mmi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
