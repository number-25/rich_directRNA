process CURL {
    tag "$download"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ricsanfre/docker-curl-jq:latest':
         'docker://ricsanfre/docker-curl-jq:latest' }"
   //     'curl/curl:8.14.1':
   //     'alpine/curl:8.14.1'}"
        //'docker://quay.io/curl/curl:8.14.1' }"

    input:
    val(prefix)
    val(suffix)
    val(url)
    //tuple val(meta), path(archive)

    output:
    path("*.zip"), emit: curl
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def retry           = task.ext.retry ?: '--retry 5'
    def mimic_browser   = task.ext.mimic_browser ?: '-A "Mozilla/5.0"'
    //def extension   = ( archive.toString() - '.gz' ).tokenize('.')[-1]
    //def name        = archive.toString() - '.gz' - ".$extension"
    //def prefix      = task.ext.prefix ?: name
    //gunzip          = prefix + ".$extension"
    """
    curl \\
        $mimic_browser \\
        $retry \\
        -L $url \\
        -o ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | cut -d" " -f2 | head -n1)
    END_VERSIONS
    """

    stub:
    def args        = task.ext.args ?: ''
    def retry       = task.ext.args ?: '--retry 5'
    //def extension   = ( archive.toString() - '.gz' ).tokenize('.')[-1]
    //def name        = archive.toString() - '.gz' - ".$extension"
    //def prefix      = task.ext.prefix ?: name
    //gunzip          = prefix + ".$extension"
    """
    touch ${prefix}.${suffix}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version | cut -d" " -f2 | head -n1)
    END_VERSIONS
    """
}
