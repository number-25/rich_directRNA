process NANOQ {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoq:0.10.0--h031d066_2' :
        'biocontainers/nanoq:0.10.0--h031d066_2'}"

    input:
    tuple val(meta), path(fastq)
    //val(output_format) //One of the following: fastq, fastq.gz, fastq.bz2, fastq.lzma, fasta, fasta.gz, fasta.bz2, fasta.lzma.

    output:
    tuple val(meta), path("*.{stats,json}")                                           , emit: stats
    tuple val(meta), path("*_verbose.stats")                                          , emit: verbose_stats
    tuple val(meta), path("*.json")                                                   ,emit: json_stats
    //tuple val(meta), path("*_filtered.${output_format}")                              , emit: reads
    path "versions.yml"                                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    //if ( "${meta.replicate}" ?: '' )
    //    def prefix = task.ext.prefix ?: "${meta.id}_nanoq" // get the sample ID from the meta mapping
    //else
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_nanoq" // get the sample ID from the meta mapping
    """
    nanoq -i ${fastq} \\
        -s -H \\ 
        > ${prefix}.stats 
    
    nanoq -i ${fastq} \\
        -s -vvv \\ 
        > ${prefix}_verbose.stats 

    nanoq -i ${fastq} 
        -s -j \\ 
        > ${prefix}.json 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$(nanoq --version | sed -e 's/nanoq //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_nanoq"
    """
    touch ${prefix}.stats
    touch ${prefix}_verbose.stats
    touch ${prefix}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$(nanoq --version | sed -e 's/nanoq //g')
    END_VERSIONS
    """
}
