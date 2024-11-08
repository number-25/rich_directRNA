process MINIMAP2_INDEX {
    tag "$fasta"
    label 'process_high'

    // Note: the versions here need to match the versions used in minimap2/align
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.28--he4a0461_0' :
        'biocontainers/minimap2:2.28--he4a0461_0' }"

    input:
    path(genome_fasta)

    output:
    //tuple path?
    path ("*.mmi"), emit: index
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def dRNA_preset = task.ext.dRNA_preset ?: "-ax splice -uf"
    def kmer = task.ext.kmer ?: "-k 14"
    """
    minimap2 \\
        ${dRNA_preset} \\
        ${kmer} \\
        -t $task.cpus \\
        -d ${genome_fasta.baseName}_k14.mmi \\
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
