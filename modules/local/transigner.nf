process TRANSIGNER {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/number_25/transigner:1.2.0':
        'number25/transigner:1.2.0' }"

    input:
    tuple val(meta), path(transcripts_fasta)
    path transcriptome_index
    val reconstruction_program
    //path transcriptome_fasta
    //val sequencing_type

    output:
    tuple val(meta), path("*.bam"),             emit: bam
    tuple val(meta), path("*.assignments.tsv"), emit: assignments
    tuple val(meta), path("*.abundances.tsv"),  emit: quant
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_${reconstruction_program}"
    """
    transigner \\
        -t $transcriptome_index \\
        $transcripts_fasta \\
        -dtype ont \\
        -p ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transigner: \$(transigner --version | sed 's#transigner ##g')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_${reconstruction_program}"

    """
    touch ${prefix}.bam
    touch ${prefix}_assignments.tsv
    touch ${prefix}_abundances.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transigner: \$(transigner --version | sed 's#transigner ##g')
    END_VERSIONS
    """
}
