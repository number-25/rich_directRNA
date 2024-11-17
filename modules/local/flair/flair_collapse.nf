process FLAIR_COLLAPSE {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://brookslab/flair:2.0.0' :
        'biocontainers/flair:2.0.0--pyhdfd78af_1' }"
        //'https://depot.galaxyproject.org/singularity/flair:2.0.0--pyhdfd78af_1':
    input:
    tuple val(meta), path(mapped_bed)
    tuple val(meta), path(fastq)
    path(annotation_gtf)
    path(genome_fasta)

    output:
    tuple val(meta), path("*.isoforms.bed"), emit: collapsed_isoforms_bed
    tuple val(meta), path("*.isoforms.gtf"), emit: collapsed_isoforms_gtf
    tuple val(meta), path("*.isoforms.fa"), emit: collapsed_isoforms_fa
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--check_splice'
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_collapsed"
    def reliant = task.ext.dRNA_preset ?: "generate"

    """
    flair \\
        collapse \\
        --annotation_reliant ${reliant} \\
        ${args} \\
        -q ${mapped_bed} \\
        -r ${fastq} \\
        -g ${genome_fasta} \\
        --gtf ${annotation_gtf} \\
        --threads $task.cpus \\
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flaircollapse: \$(flair --version |& sed 's/flair //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_collapsed"
    """
    touch ${prefix}.isoforms.bed
    touch ${prefix}.isoforms.gtf
    touch ${prefix}.isoforms.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flaircollapse: \$(flair --version |& sed 's/flair //')
    END_VERSIONS
    """
}
