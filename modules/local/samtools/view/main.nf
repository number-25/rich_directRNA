process SAMTOOLS_VIEW {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.20-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam") , emit: unmapped_bam
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-h -f 4 -b'
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}"

    """
    samtools \\
        view \\
        $args \\
        -@ ${task.cpus-1} \\
        ${bam} \\
        -o "$prefix"_unmapped.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: '-h -f 4 -b'
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}"

    """
    touch "${prefix}"_ummapped.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
