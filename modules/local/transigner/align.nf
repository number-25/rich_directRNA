    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/number_25/transigner-alpine:1.1.3':
        'number25/transigner-alpine:1.1.3' }"

    //TO-DO
    // Split the different chunks of the program into separate modules and wrap into a subworkflow?
    // Need a MAPPED_BAM_TO_FASTQ module to convert the mapped bam reads to
    // fastq format, to then use here. - this is assuming that we're not taking
    // the fasta reads that we're reconstructed, or even downstream recovered
    // from sqanti

    input:
    tuple val(meta), path(transcripts)
    path transcriptome_fasta

    output:
    tuple val(meta), path("*.bam", temporary: true), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    def outdir = task.ext.args ?: './'
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_transigner"
    """
    transigner \\
        align \\
        ${args} \\
        -q ${transcripts} \\
        -t ${transcriptome_fasta} \\
        -d . \\
        -n 181
        -o ${prefix}.bam \\
        -p $tasks.cpu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transigner: \$(echo 1.1.0)
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transigner: \$(echo 1.1.0)
    END_VERSIONS
    """
}
