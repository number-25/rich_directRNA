process FLAIR_CORRECT {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://brookslab/flair:2.0.0' :
        'biocontainers/flair:2.0.0--pyhdfd78af_1' }"

    input:
    // Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(bed)
    path (genome_fasta)
    path (annotation_gtf)

    // FLAIR correct has many outputs
    output:
    tuple val(meta), path("*_corrected.bed"), emit: flair_corrected_bed
    tuple val(meta), path("*_inconsistent.bed"), emit: flair_inconsistent_bed
    tuple val(meta), path("*_verify.bed"), emit: flair_unverified_bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_flair_correct"
    """
    flair \\
        correct \\
        -q ${bed} \\
        -f ${annotation_gtf} \\
        -g ${genome_fasta} \\
        --threads $task.cpus \\
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flaircorrect: \$(flair --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_flair"
    """
    touch ${prefix}_all_corrected.bed
    touch ${prefix}_all_inconsistent.bed
    touch ${prefix}_cannot_verify.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flaircorrect: \$(flair --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
