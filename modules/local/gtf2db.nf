process GTF2DB {
    tag "isoquant"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/number_25/gtf2db_iso:3.7.1' :
        'number25/gtf2db_iso:3.7.1'}"

    input:
    path transcriptome_annotation

    output:
    path("*.db"), emit: isoquant_database

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python3 gtf2db.py \\
        -i $transcriptome_annotation \\
        --complete_genedb \\
        -o ${transcriptome_annotation.baseName}.db
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch ${transcriptome_annotation.baseName}.db
    """
}
