process SEQKIT_STATS {
    tag "${meta.alias}"
    label 'process_low'

    container "quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0"
    input:
    tuple val(meta), path(reads_list)

    output:
    tuple val(meta), path("*.tsv"), emit: stats
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.alias}"
    """
    mkdir all_fastq
    mv $reads_list all_fastq
    
    seqkit stats \\
        --tabular \\
        --threads ${task.cpus} \\
        ${args} \\
        all_fastq/*.fastq.gz > '${prefix}_stats.tsv'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.alias}"
    """
    touch ${prefix}_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}
