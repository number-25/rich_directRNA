process JAFFAL {
    echo true
    label 'process_high'
    conda "bioconda::jaffa=2.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://davidsongroup/jaffa:2.4-devel' :
        'docker://davidsongroup/jaffa:2.4-devel' }"

    //singularityRunOptions = '-B $jaffal_ref_dir:/ref'
    singularityRunOptions = 'exec'

    input:
    tuple val(meta), path(fastq)
    path jaffal_ref_dir

    output:
    tuple val(meta), path("*.fasta"), emit: jaffal_fasta
    path "*.csv"                    , emit: jaffal_results
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when
    // prefix??

    //refBase=$jaffal_ref_dir \\
    //-p genome=Masked_hg38 \\
    //-p annotation=hg38_genCode37 \\
    script:
    """
    bpipe \\
        run \\
        -p refBase=$jaffal_ref_dir \\
        --threads $task.cpus \\
        /JAFFA/JAFFAL.groovy \\
        $fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jaffa: \$( echo '2.3' )
    END_VERSIONS
    """
}

