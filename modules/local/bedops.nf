process BEDOPS {
    tag "$phylop"
    label 'process_high'
    //conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedops:2.4.41--h4ac6f70_2':
        'biocontainers/bedops:2.4.41--h4ac6f70_2' }"

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    path phylop_wig

    output:
    //  nf-core: Named file extensions MUST be emitted for ALL output channels
    path("*.bed"), emit: phylop_bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //def args = task.ext.args ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}"
    """
    wig2bed \\
        --max-mem 70G \\
        --zero-indexed \\
        < $phylop_wig \\
        > ${phylop_bed}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedops: \$(convert2bed -w | tail -n 2 | head -n 1 | cut -d: -f2)
    END_VERSIONS
    """

    stub:
    //def args = task.ext.args ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${phylop_bed}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedops: \$(convert2bed -w | tail -n 2 | head -n 1 | cut -d: -f2)
    END_VERSIONS
    """
}
