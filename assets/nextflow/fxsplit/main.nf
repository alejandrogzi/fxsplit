process FXSPLIT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' :
        'ghcr.io/alejandrogzi/fxsplit:latest' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("chunks/*fa.gz")      , optional: true, emit: fasta_gz
    tuple val(meta), path("chunks/*.fa")        , optional: true, emit: fasta
    tuple val(meta), path("chunks/*.fq.gz")     , optional: true, emit: fastq_gz
    tuple val(meta), path("chunks/*.fq")        , optional: true, emit: fastq
    tuple val(meta), path("chunks/*.2bit")      , optional: true, emit: twobit
    path  "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def chunks        = task.ext.chunks ?: 400000
    def gzip          = reads.name.endsWith('.gz') ? true : false
    """
    fxsplit \\
        $args \\
        -f $reads \\
        -c $chunks \\
        -t $task.cpus \\
        -C \\
        --suffix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fxsplit: \$( fxsplit --version | head -n 1 | sed 's/fxsplit //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gzip   = reads.name.endsWith('.gz') ? true : false
    """
    touch chunks/gz/${prefix}.gz
    touch chunks/fa/${prefix}.fa
    touch chunks/fq/${prefix}.fq
    touch chunks/2bit/${prefix}.2bit

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fxsplit: \$( fxsplit --version | head -n 1 | sed 's/fxsplit //g' | sed 's/ (.*//g' )
    END_VERSIONS
    """
}
