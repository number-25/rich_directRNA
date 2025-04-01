    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'docker/number25/transigner:' }"

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

// add -dtype $data_type options

    script:
    def align_args = task.ext.args ?: ''
    def prefilter_args = task.ext.args ?: ''
    def em_args = task.ext.args ?: ''
    def outdir = task.ext.args ?: './'
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_transigner"
    """
    transigner \\
        align \\
        ${align_args} \\
        -q ${transcripts} \\
        -t ${transcriptome_fasta} \\
        -d . \\
        -o ${prefix}.bam \\
        -p $tasks.cpu

    transigner \\
        pre \\
        ${prefilter_args} \\
        -i ${prefix}.bam \\
        -d .

    transigner \\
        em \\
        ${em_args} \\
        -s ./${prefix}_scores.tsv \\
        -u ./${prefix}_unmapped.tsv \\
        -m ./${prefix}_tmap.csv \\
        -p $tasks.cpu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transigner: \$(echo 1.1.0)
    END_VERSIONS
    """
    //transigner: \$(transigner --version |& sed '1!d ; s/samtools //')
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
