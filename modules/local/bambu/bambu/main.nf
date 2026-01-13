process BAMBU {
    tag "$meta.id"
    label 'process_medium'
    conda "conda-forge::r-base=4.0.3 bioconda::bioconductor-bambu=3.0.8 bioconda::bioconductor-bsgenome=1.66.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://number25/bambu:3.12.0' :
    'docker://quay.io/number_25/bambu:3.12.0' }"

    input:
    path(genome_fasta)
    path(annotation_gtf)
    tuple val(meta), path(bam)

    output:
    path "*.gtf"
    tuple val(meta), path("supportedTranscriptModels.gtf")  , emit: bambu_supported_gtf
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_bambu.r \\
        --tag=. \\
        --ncore=$task.cpus \\
        --annotation=${annotation_gtf} \\
        --fasta=${genome_fasta} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-bambu: \$(Rscript -e "library(bambu); cat(as.character(packageVersion('bambu')))")
        bioconductor-bsgenome: \$(Rscript -e "library(BSgenome); cat(as.character(packageVersion('BSgenome')))")
    END_VERSIONS
    """
}
