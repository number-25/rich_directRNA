process STRINGTIE {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/stringtie:2.2.3--h43eeafb_0' :
        'biocontainers/stringtie:2.2.3--h43eeafb_0' }"

    input:
    tuple val(meta), path(bam)
    path annotation_gtf

    output:
    tuple val(meta), path("*.transcripts.gtf"), emit: stringtie_gtf
    tuple val(meta), path("*.coverage.gtf"),    emit: stringtie_coverage
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}.${meta.replicate}.stringtie"
    def reference   = annotation_gtf ? "-G $annotation_gtf" : ""
    def coverage    = annotation_gtf ? "-C ${prefix}.coverage.gtf" : ""
    """
    stringtie \\
        $bam \\
        $reference \\
        $coverage \\
        -L \\
        -p $task.cpus \\
        -o ${prefix}.transcripts.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version 2>&1)
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}.${meta.replicate}.stringtie"
    def reference   = annotation_gtf ? "-G $annotation_gtf" : ""
    def coverage    = annotation_gtf ? "-C ${prefix}.coverage.gtf" : ""
    """
    touch ${prefix}.transcripts.gtf
    touch ${prefix}.coverage.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: \$(stringtie --version 2>&1)
    END_VERSIONS
    """
}
