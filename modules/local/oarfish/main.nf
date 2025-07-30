process OARFISH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/oarfish:0.9.0--h5ca1c30_0':
        'biocontainers/oarfish:0.9.0--h5ca1c30_0' }"

    input:
    tuple val(meta), path(transcripts_fasta)
    path transcriptome_index
    val sequencing_type

    output:
    tuple val(meta), path("*features.tsv.gz") , emit: features
    //tuple val(meta), path("*barcodes.tsv.gz") , emit: barcodes
    tuple val(meta), path("*matrix.mtx.gz")    , emit: mtx
    tuple val(meta), path("*meta_info.json")  , emit: meta_info
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                = task.ext.args ?: ''
    def prefix              = task.ext.prefix ?: "${meta.id}"
    def secondary_mappings  = task.ext.secondary_mappings ?: '--best-n 180'
    def annotation          = task.ext.annotation ?: "--annotated $transcriptome_fasta"
    def technology          = task.ext.technology ?: "--seq-tech $sequencing_type"

    """
    oarfish \\
        $annotation \\
        --output ${prefix} \\
        --alignments $transcripts_fasta \\
        --threads ${task.cpus} \\
        --index $transcriptome_index \\
        --model-coverage \\
        $secondary_mappings \\
        $technology \\
        ${args}

    mv *features.txt features.tsv
    //mv *barcodes.txt barcodes.tsv

    grep '^%' *count.mtx > matrix.mtx
    grep -v '^%' *count.mtx | awk '{print \$2" "\$1" "\$3}' >> matrix.mtx

    for tsv_file in *features.tsv *barcodes.tsv *matrix.mtx
    do
        gzip \$tsv_file
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oarfish: \$(oarfish --version | sed 's#oarfish ##g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.features.tsv.gz
    touch ${prefix}.barcodes.tsv.gz
    touch ${prefix}.matrix.mtx.gz
    touch ${prefix}.meta_info.json
    """
}
