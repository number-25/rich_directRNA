process GET_SYLPH_TAX_DB {
    tag "SYLPH_TAX"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sylph-tax:1.8.0--pyhdfd78af_0' :
        'biocontainers/sylph-tax:1.8.0--pyhdfd78af_0' }"

    //input:
    //val outpath

    output:
    path "tax_db"       , emit: sylph_tax_db
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}"

    """
    mkdir tax_db
    sylph-tax download --download-to tax_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph_tax: \$(sylph-tax --version | cut -d" " -f2)
    END_VERSIONS
    """

/*    stub:
    //def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}"

    """
    touch *.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sylph_tax: \$(sylph-tax --version | cut -d" " -f2)
    END_VERSIONS
    """
*/
}
