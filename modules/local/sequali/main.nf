process SEQUALI {
    tag "$meta.id"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequali:0.12.0--py310h1fe012e_1':
        'biocontainers/sequali:0.12.0--py310h1fe012e_1' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.json"), emit: sequali_json
    tuple val(meta), path("*.html"), emit: sequali_html
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_sequali"
    def html_output = task.ext.html_output ?: "--html ${prefix}.html"
    def json_output = task.ext.json_output ?: "--json ${prefix}.json"

    """
    sequali \\
        $args \\
        $html_output \\
        $json_output \\
        ${fastq}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequali: \$(sequali --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_sequali"
    def html_output = task.ext.html_output ?: "--html ${prefix}.html"
    def json_output = task.ext.json_output ?: "--json ${prefix}.json"
    """
    touch ${prefix}.json
    touch ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequali: \$(sequali --version)
    END_VERSIONS
    """
}
