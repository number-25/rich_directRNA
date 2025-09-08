process SQANTI_QC {
    tag "$meta.id"
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'anaconesalab/sqanti3:v5.5.1"

    input:
    tuple val(meta), path(reconstructed_transcriptome)
    path annotation_gtf
    tuple path(genome_fasta), path(genome_fasta_index)
    val reconstruction_program

    output:
    path "isoforms_classification.txt"
    path "isoforms_corrected.cds.gff3"
    path "isoforms_corrected.faa"
    path "isoforms_corrected.fasta"
    path "isoforms_corrected.genePred"
    tuple val(meta), path "isoforms_corrected.gtf", emit: sqanti_qc_isoforms_gtf
    path "isoforms_junctions.txt"
    path "isoforms.qc_params.txt"
    path "refAnnotation_isoforms.genePred"
    path "unknown_strand.gtf"
    //path "logs/final_report.log"
    //path "logs/gtf2fasta.log
    //path "logs/GTF_to_genePred.log"
    //path "logs/normalize_gtf.log"
    //path "logs/sqanti3_qc.log"
    //path "logs/TD2_LongOrfs.log"
    //path "logs/TD2_Predict.log"
    path "RTS/sj.rts.results.tsv"                , optional: yes
    path "TD2/isoforms_corrected.fasta.TD2.bed"  , optional: yes
    path "TD2/isoforms_corrected.fasta.TD2.cds"  , optional: yes
    path "TD2/isoforms_corrected.fasta.TD2.gff2" , optional: yes
    path "TD2/longest_orfs.cds"                  , optional: yes
    path "TD2/longest_orfs.gff2"                 , optional: yes
    path "TD2/longest_orfs.pep"                  , optional: yes
    path "TD2/psauron_score.csv"                 , optional: yes
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def output_name = task.ext.out_name ?: "--output ${program}"
    def prefix      = task.ext.prefix ?:
    "${meta.id}.${meta.replicate}_${reconstruction_ program}_sqantiQC"
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
    def args        = task.ext.args ?: ''
    def output_name = task.ext.out_name ?: "--output ${program}"
    def prefix      = task.ext.prefix ?: "${meta.id}.${meta.replicate}_${program}_sqantiQC"
    """
    touch isoforms_classification.txt
    touch isoforms_corrected.cds.gff3
    touch isoforms_corrected.faa
    touch isoforms_corrected.fasta
    touch isoforms_corrected.genePred
    touch isoforms_corrected.gtf
    touch isoforms_junctions.txt
    touch isoforms.qc_params.txt
    touch refAnnotation_isoforms.genePred
    touch unknown_strand.gtf
    mkdir logs
    mkdir RTS
    touch RTS/sj.rts.results.tsv
    mkdir TD2
    touch TD2/isoforms_corrected.fasta.TD2.bed
    touch TD2/isoforms_corrected.fasta.TD2.cds
    touch TD2/isoforms_corrected.fasta.TD2.gff2
    touch TD2/longest_orfs.cds
    touch TD2/longest_orfs.gff2
    touch TD2/longest_orfs.pep
    touch TD2/psauron_score.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3_qc: \$(sqanti3_qc.py | cut -d" " -f2')
    END_VERSIONS
    """
}
