process ISOQUANT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/isoquant:3.6.2--hdfd78af_0':
        'biocontainers/isoquant:3.6.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path(annotation_gtf)
    path(genome_fasta)

    output:
    tuple val(meta), path("*.read_assignments.tsv.gz"),
    tuple val(meta), path("*.corrected_reads.bed.gz"),
    tuple val(meta), path("*.transcript_tpm.tsv"),
    tuple val(meta), path("*.transcript_counts.tsv"),
    tuple val(meta), path("*.gene_tpm.tsv"),
    tuple val(meta), path("*.gene_counts.tsv"),
    tuple val(meta), path("*.transcript_models.gtf"), emit: isoquant_transcript_gtf
    tuple val(meta), path("*.transcript_model_reads.tsv.gz"), emit: isoquant_transcript_models
    tuple val(meta), path("*.transcript_model_tpm.tsv"),
    tuple val(meta), path("*.transcript_model_counts.tsv"),
    tuple val(meta), path("*.extended_annotation.gtf"), emit isoquant_new_reference_transcriptome_gtf, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.replicate}"_isoquant
    def dRNA_preset = task.ext.dRNA_preset ?: "-d nanopore"
    def strand_preset = task.ext.dRNA_preset ?: "--stranded none"
    def input_bam = task.ext.kmer ?: "--bam $bam"
    def ref_genome = task.ext.kmer ?: "--reference $genome_fasta"
    def ref_gtf = task.ext.kmer ?: "--genedb $annotation_gtf"
    def complete = task.ext.kmer ?: "--complete_genedb"
    //def output = task.ext.kmer ?: "--output $annotation_gtf"
    """
    isoquant.py \\
        $dRNA_preset \\
        $strand_preset \\
        $input_bam \\
        $ref_genome \\
        $ref_gtf \\
        $complete \\
        --prefix $prefix \\
        --threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
