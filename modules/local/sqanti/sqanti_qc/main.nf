process QC_SQANTI {
    tag "$meta.id"
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'anaconesalab/sqanti3:5.3.6-conda-fix'"

    input:
    tuple val(meta), path(reconstructed_transcriptome)
    path annotation_gtf
    path genome_fasta
    val program

    output:
    // Substantial outputs list
    path "refAnnotation_*.genePred"
    path "*.params.txt"
    path "*_classification.txt"
    path "*_corrected.faa"
    path "*_corrected.fasta"
    path "*_corrected.genePred
    path "*.corrected.gtf"
    path "*.corrected.gtf.cds.gff"
    path "*.html", optional: yes
    path "*.pdf", optional: yes
    path "*.junctions.txt", optional: yes
    path "unknown_strand.gtf"
    path "GMST"
    path "RTS"

    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def output_name = task.ext.out_name ?: "--output ${program}"
    def prefix      = task.ext.prefix ?: "${meta.id}.${meta.replicate}_${program}_sqantiQC"
    """
    sqantiqc.py \\
        $reconstructed_transcriptome \\
        $annotation_gtf \\
        $genome_fasta \\
        $args \\
        $output_name \\
        $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3_qc: \$(sqanti3_qc.py | cut -d" " -f2')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    mkdir GMST
    touch GMST/gms.log
    touch GMST/GMST_tmp
    touch GMST/GMST_tmp.faa
    touch GMST/GMST_tmp.fnn
    mkdir RTS
    touch RTS/sj.rts.results.tsv
    touch unknown_strand.gtf
    touch ${program}_SQANTI3_report.html
    touch ${program}_SQANTI3_report.pdf
    touch ${program}_junctions.txt
    touch ${program}_corrected.gtf.cds.gff
    touch ${program}_corrected.gtf
    touch ${program}_corrected.genePred
    touch ${program}_corrected.fasta
    touch ${program}_corrected.faa
    touch ${program}_classification.txt
    touch ${program}.params.txt
    touch refAnnotation_${program}.genePred

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3_qc: \$(sqanti3_qc.py | cut -d" " -f2')
    END_VERSIONS
    """
}
