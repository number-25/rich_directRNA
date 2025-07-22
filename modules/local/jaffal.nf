process JAFFAL {
    echo true
    label 'process_high'

    conda "bioconda::jaffa=2.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jaffa:2.3--hdfd78af_0' :
        'quay.io/biocontainers/jaffa:2.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)
    path jaffal_ref_dir

    output:
    tuple val(meta), path("*.fasta"), emit: jaffal_fasta
    path "*.csv"                    , emit: jaffal_results
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    //refBase=$jaffal_ref_dir \\
    //-p genome=Masked_hg38 \\
    //-p annotation=hg38_genCode37 \\
    script:
    """
    bpipe \\
        run \\
        JAFFAL.groovy \\
        -p refBase=$jaffal_ref_dir \\
        -n $task.cpus \\
        $fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jaffa: \$( echo '2.3' )
    END_VERSIONS
    """
}

