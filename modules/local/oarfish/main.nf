process OARFISH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/oarfish:0.9.0--h5ca1c30_0':
        'biocontainers/oarfish:0.9.0--h5ca1c30_0' }"

    input:
    tuple val(meta), path(transcripts_fasta)
    path transcriptome_index
    val sequencing_type

    output:
    tuple val(meta), path("*.quant.gz"),        emit: quant
    tuple val(meta), path("*.meta_info.json"),  emit: meta_info
    path "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                = task.ext.args ?: ''
    def prefix              = task.ext.prefix ?: "${meta.id}_${meta.replicate}_$sequencing_type"
    def technology          = task.ext.technology ?: "--seq-tech $sequencing_type"
    def secondary_mappings  = task.ext.secondary_mappings ?: "--best-n 180"
    def filters             = task.ext.filters ?: "--filter-group no-filters"
    """
    oarfish \\
        --reads $transcripts_fasta \\
        --index $transcriptome_index \\
        $technology \\
        $secondary_mappings \\
        --output ${prefix} \\
        --threads ${task.cpus} \\
        --model-coverage \\
        $filters \\
        ${args}

    gzip *.quant

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oarfish: \$(oarfish --version | sed 's#oarfish ##g')
    END_VERSIONS
    """

    stub:
    def args                = task.ext.args ?: ''
    def prefix              = task.ext.prefix ?: "${meta.id}_${meta.replicate}_$sequencing_type"
    def technology          = task.ext.technology ?: "--seq-tech $sequencing_type"
    def secondary_mappings  = task.ext.secondary_mappings ?: "--best-n 180"
    def filters             = task.ext.filters ?: "--filter-group no-filters"
    """
    touch ${prefix}.meta_info.json
    touch ${prefix}.quant.gz
    """
}
