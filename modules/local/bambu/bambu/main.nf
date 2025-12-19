process BAMBU {
    tag "$meta.id"
    label 'process_medium'
    conda "conda-forge::r-base=4.0.3 bioconda::bioconductor-bambu=3.0.8 bioconda::bioconductor-bsgenome=1.66.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'number25/bambu:3.12.0' :
        'number25/bambu:3.12.0' }"

       // 'docker://quay.io/biocontainers/bioconductor-bambu:3.0.8--r42hc247a5b_0' }"
       // 'docker://quay.io/biocontainers/bioconductor-bambu:3.0.8--r42hc247a5b_0' }"
   //containerOptions = '-u $(id -u):$(id -g)'
        //'number25/bambu:3.8.0':
        //'docker://quay.io/number_25/bambu:latest' }"
        // 'https://depot.galaxyproject.org/singularity/bioconductor-bambu:3.4.0--r43hf17093f_1' :
        //'quay.io/biocontainers/bioconductor-bambu:3.4.0--r43hf17093f_1' }"

    input:
    path(genome_fasta)
    //path(genome_fasta_sizes)
    path(annotation_gtf)
    tuple val(meta), path(bam)

    output:
    //path "counts_gene.txt"                                  , emit: ch_gene_counts
    //path "counts_transcript.txt"                            , emit: ch_transcript_counts
    path "extended_annotations.gtf"                         , emit: bambu_extended_gtf
    path "allTranscriptModels.gtf"                          , emit: bambu_all_gtf
    path "novelTranscripts.gtf"                             , emit: bambu_novel_only_gtf
    tuple val(meta), path("supportedTranscriptModels.gtf")  , emit: bambu_supported_gtf
    path "versions.yml"                                     , emit: versions

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
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-bambu: \$(Rscript -e "library(bambu); cat(as.character(packageVersion('bambu')))")
        bioconductor-bsgenome: \$(Rscript -e "library(BSgenome); cat(as.character(packageVersion('BSgenome')))")
    END_VERSIONS
    """
}
