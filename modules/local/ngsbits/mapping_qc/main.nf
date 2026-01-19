process NGS_BITS {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngs-bits:2025_01--py313h6fdbb3c_0':
        'benekenobi/ngs-bits:2025_01' }"

    input:
    tuple val(meta), path(bam), path(index)
    path genome_fasta
    val build
    val contamination

    output:
    tuple val(meta), path("*.qcML"),    emit: qcML
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args ?: ''
    def contamination   = task.ext.contamination ?: "-no_cont"
    def genome_build    = task.ext.build ?: "$build"
    def prefix          = task.ext.prefix ?: "${meta.id}_${meta.replicate}_ngsbits"

    // -build $genome_build
    """
    MappingQC \\
        -in $bam \\
        -out ${prefix}.qcML \\
        -rna \\
        -ref $genome_fasta \\
        -long_read \\
        $contamination

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngsbits: \$(MappingQC --version)
    END_VERSIONS
    """

    stub:
    def args          = task.ext.args ?: ''
    def contamination = task.ext.contamination ? "-no_cont" : ""
    def build         = task.ext.build ?: "$build"
    def prefix        = task.ext.prefix ?: "${meta.id}_${meta.replicate}_ngsbits"

    """
    touch ${prefix}.qcML

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngsbits: \$(MappingQC --version)
    END_VERSIONS
    """
}
