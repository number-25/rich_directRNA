process CHECK_SAMPLESHEET {
    tag '$samplesheet'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/julia:1.10':
        'biocontainers/julia:1.10' }"

    input:
    path samplesheet

    output:
    path samplesheet, emit: csv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    check_samplecsv.jl $samplesheet
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        julia: \$(julia --version | sed 's/julia verison //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''

    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        julia: \$(julia --version | sed 's/julia verison //g')
    END_VERSIONS
    """
}
