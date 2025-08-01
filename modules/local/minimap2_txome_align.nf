process MINIMAP2_TXOME_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

    input:
    tuple val(meta), path(fasta)
    path transcriptome_index
    val reconstruction_program

    output:
    tuple val(meta), path("*.sam")                       , optional: true, emit: sam
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // This can be expanded eventually to allow cDNA mapping, etc.
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_${reconstruction_program}_minimap2"
    def dRNA_preset = task.ext.dRNA_preset ?: "-ax splice -uf"
    def kmer = task.ext.kmer ?: "-k 14"

    """
    minimap2 \\
        ${dRNA_preset} \\
        ${kmer} \\
        -N 180 \\
        -t ${task.cpus} \\
        ${transcriptome_index} \\
        $fasta \\
        > ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_${reconstruction_program}_minimap2"

    """
    touch ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
