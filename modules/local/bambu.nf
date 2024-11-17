process BAMBU {
    tag "$meta.id"
    label 'process_medium'
    conda "conda-forge::r-base=4.0.3 bioconda::bioconductor-bambu=3.4.0 bioconda::bioconductor-bsgenome=1.74.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-bambu:3.4.0--r43hf17093f_1' :
        'quay.io/biocontainers/bioconductor-bambu:3.4.0--r43hf17093f_1' }"

    input:
    path(genome_fasta)
    //path(genome_fasta_sizes)
    path(annotation_gtf)
    tuple val(meta), path(collapsed_bams)

    output:
    tuple val(meta), path("*_counts_gene.txt")         , emit: bambu_gene_counts
    tuple val(meta), path("_counts_transcript.txt")   , emit: bambu_transcript_counts
    tuple val(meta), path("_extended_annotations.gtf"), emit: bambu_extended_gtf
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when
    //def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}_bambu"

    script:
    """
    run_bambu.r \\
        --tag=. \\
        --ncore=$task.cpus \\
        --annotation=${annotation_gtf} \\
        --fasta=${genome_fasta} \\
        ${collapsed_bams}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-bambu: \$(Rscript -e "library(bambu); cat(as.character(packageVersion('bambu')))")
        bioconductor-bsgenome: \$(Rscript -e "library(BSgenome); cat(as.character(packageVersion('BSgenome')))")
    END_VERSIONS
    """
}

