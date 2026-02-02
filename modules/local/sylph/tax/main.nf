process SYLPH_TAX {
    tag "SYLPH_TAX"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sylph-tax:1.8.0--pyhdfd78af_0' :
        'biocontainers/sylph-tax:1.8.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(sylph_profile)
    val database_name

    output:
    tuple val(meta), path("*.sylphmpa")     , emit: sylph_tax
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}"

    """
    sylph-tax download --download-to .

    sylph-tax taxprof \\
        $sylph_profile \\
        --taxonomy-dir !PWD \\
        -t $database_name \\
        -o ${prefix}_sylph_tax.sylphmpa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph_profile: \$(sylph-tax --version | cut -d" " -f2')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}"

    """
    touch ${prefix}_sylph_tax.sylphmpa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph_profile: \$(sylph-tax --version | cut -d" " -f2')
    END_VERSIONS
    """
}
