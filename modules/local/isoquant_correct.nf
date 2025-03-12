process ISOQUANT {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoquant:3.6.2--hdfd78af_0':
        'biocontainers/isoquant:3.6.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path(genome_fasta)
    path(annotation_gtf)
    val

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    isoquant.py \\
        -d nanopore \\
        --report_noveL_unspliced true \\
        --reference ${genome_fasta} \\
        --bam $bam \\
        --threads $task.cpus \\




    samtools \\
        sort \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.bam \\
        -T $prefix \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
