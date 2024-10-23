//               Any parameters that need to be evaluated in the context of a particular sample
// nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
//  nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process BIGWIG_TO_WIG {
    tag "$phylop"
    label 'process_medium'
    //conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-bigwigtowig:469--h9b8f530_0':
        'biocontainers/ucsc-bigwigtowig' }"

    input:
    path phylop_bigwig

    output:
    path("*.wig"), emit: phylop_wig

    when:
    task.ext.when == null || task.ext.when

    script:
    //def args = task.ext.args ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bigWigToWig $phylop_bed ${phylop_wig}.wig
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${phylop_wig}.wig
    """
}
