process SYLPH_PREPARE_REFERENCE {
    tag "SYLPH_REFERENCE"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ricsanfre/docker-curl-jq:latest':
        'docker://ricsanfre/docker-curl-jq:latest' }"

    input:
    val url

    output:
    path("*.syldb"),        emit: database
    path "versions.yml",    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    //def retry           = task.ext.retry ?: '--retry 5'
    def mimic_browser   = task.ext.mimic_browser ?: '-A "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"'
    //def extension   = ( archive.toString() - '.gz' ).tokenize('.')[-1]
    //def name        = archive.toString() - '.gz' - ".$extension"
    //def prefix      = task.ext.prefix ?: name
    //gunzip          = prefix + ".$extension"
    """
    curl \\
        $mimic_browser \\
        -L $url

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | cut -d" " -f2 | head -n1)
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args ?: ''
    //def retry       = task.ext.args ?: '--retry 5'
    def mimic_browser   = task.ext.mimic_browser ?: '-A "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"'
    //def extension   = ( archive.toString() - '.gz' ).tokenize('.')[-1]
    //def name        = archive.toString() - '.gz' - ".$extension"
    //def prefix      = task.ext.prefix ?: name
    //gunzip          = prefix + ".$extension"
    """
    touch "*.syldb"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | cut -d" " -f2 | head -n1)
    END_VERSIONS
    """
}
