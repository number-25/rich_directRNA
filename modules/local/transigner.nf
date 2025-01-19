    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"
    //TO-DO
    // Split the different chunks of the program into separate modules and wrap into a subworkflow?
    // Need a MAPPED_BAM_TO_FASTQ module to convert the mapped bam reads to fastq format, to then use here.

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam", temporary: true), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def outdir = task.ext.args ?: './'
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}"
    """
    easy.sh \\
        ${sample_fastq} \\
        ${transcriptome_fasta} \\
        ${outdir}

    transigner \
    align \
    -q ${sample_fastq} \ ? what is it called in the minimap2 module
    -t ${transcriptome} \
    -d . \
    -o ${prefix}.bam \
    -p $tasks.cpu


    transigner \
    prefilter \
    -a ${bam} \
    -t transcripts.fa \
    -o output_dir \
    --filter \
    -tp -1 | \
    transigner \
    em \
    -s output_dir/scores.tsv \
    -i output_dir/ti.pkl \
    -o output_dir \
    --drop \
    --use-score

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transigner: 1.1.0
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
        transigner: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
