    tag "$meta.id"
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'anaconesalab/sqanti3:5.3.6-conda-fix'"

    input:
    tuple val(meta), path(reconstructed_transcriptome)
    path annotation_gtf
    path genome_fasta
    val cage
    path cage_path
    val polyA_motif
    path polyA_motif_path
    val polyA_sites
    path polyA_sites_path
    val intron_junctions
    path introns_junctions_path


    output:
    // Substantial outputs list


    path "*.html", optional:yes
    path "*.pdf", optional:yes

    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                    = task.ext.args ?: ''
    def with_cage               = cage ? "${cage_path}" : ''
    def with_polyA_motif        = polyA_motif ? "${polyA_motif_path}" : ''
    def with_polyA_sites        = polyA_sites ? "${polyA_sites_path}" : ''
    def with_intron_junctions   = intron_junctions ? "${intron_junctions_path}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.replicate}_sqantiQC"
    """
    sqantiqc.py \\
        $reconstructed_transcriptome \\
        $annotation_gtf \\
        $genome_fasta \\
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
