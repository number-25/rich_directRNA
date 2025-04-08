process QC_SQANTI {
    tag "$meta.id"
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'anaconesalab/sqanti3:5.3.6-conda-fix'"

    input:
    tuple val(meta), path(reconstructed_transcriptome)
    path annotation_gtf
    path genome_fasta
    val program

    output:
    // Substantial outputs list
    path "*.html", optional: yes
    path "*.pdf", optional: yes

    path GMST
    path

    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                    = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.replicate}_sqantiQC"
    """
    sqantiqc.py \\
        $reconstructed_transcriptome \\
        $annotation_gtf \\
        $genome_fasta \\
        $args


        $with_cage \\
        $with_polyA_motif \\
        $with_polyA_sites \\
        $with_intron_junctions \\
        $report \\
        $task.cpus \\
        $args \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
